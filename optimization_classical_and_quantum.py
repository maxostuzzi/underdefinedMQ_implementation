r"""

Description:

    The objective of the present algorithm is to generate all possible partitions of the input m. 
    Subsequent to this, filter mechanisms will be applied in accordance with the constraints described in Sections 3 and 4 of the paper. 
    Moving forward, the subsequent step involves the evaluation of the partitions to obtain the corresponding bit complexity. 
    The result of this process is a list of partitions, with the most efficient one appearing first in the list.

Input:

    n   number of variables
    m   number of equations
    q   characteristic of finite field F_q
    p   partition size
    c   number of cpus for computation

Optional input:

    desc            description used for output files
    t               threshold for filtering partitions, is used to skip partitions for which the exhaustive search bit complexity exceeds $t$
    maxsummand      constraint s.t b_2, ... , b_p <= maxsummand
    maxb1           constraint s.t b_1 <= maxb1
    minsumguess     k > minsumguess
    lines           limits the number of results in the output file, recommended is for instance 1000, larger values increase RAM usage!
    
Example usage: 

    python .\optimization_classical_and_quantum.py -n 860 -m 78 -q 16 -p 4 -c 8 -desc mayo_r2_sec_lvl_1_ -t 157 -maxsummand 78 -maxb1 78 -minsumguess 1 -lines 1000

Note:

    Increasing the "lines" parameter also increases the RAM usage!

Output:

    Test results classical, quantum and statistic on conditions
"""
import argparse
import heapq
import itertools
import math
from multiprocessing import Pool, cpu_count, Manager
from cryptographic_estimators.MQEstimator import MQEstimator
from tqdm.contrib.concurrent import process_map

def is_valid_partition(part, m, p, MAX_b1, MAX_SUMMAND, MIN_SUM_GUESS):
    r"""
    This function verifies the validity of a partition. 
    Involves verifying that a partition is equal to m, that the first element does not exceed MAX_b1, 
    all inner elements (b2, ..., bp) must not exceed MAX_SUMMAND, and that the last elements (k = k1 + k2) must exceed MIN_SUM_GUESS.
    """
    if sum(part) != m or part[0] > MAX_b1:
        return False
    if any(x > MAX_SUMMAND for x in part[1:p - 2]):
        return False
    if part[p - 2] + part[p - 1] < MIN_SUM_GUESS:
        return False
    return True

def filter_test(n, m, partition, threshold, exponent):
    r"""
    Application of filter mechanisms in accordance with the specifications outlined in sections 3 and 4 of the paper.
    """
    k_1, k_2 = partition[-2], partition[-1]
    p = len(partition) - 2
    conditions = {"condition_1": 0, "condition_2": 0, "condition_3": 0, "condition_4": 0, "condition_5": 0}

    for i in range(1, partition[0]):
        if ((i * (sum(partition[:-2]) + k_2 - i)) + i) > n:
            conditions["condition_1"] += 1
            return False, conditions

    if len(partition) - 2 > 2:
        for i in range(partition[1]):
            if ((partition[0] + i) * (sum(partition[1:-2]) - i) + i * partition[0] + i) > n:
                conditions["condition_2"] += 1
                return False, conditions

        for l in range(2, p - 1):
            sum_1 = sum(partition[:l - 1])
            sum_2 = sum(partition[l - 1:p])
            sum_3 = sum(partition[j] * partition[h] for j in range(l - 1) for h in range(j + 1, l - 1))
            for i in range(partition[l]):
                if ((sum_1 + i) * (sum_2 - i) + i * sum_1 + sum_3 + i) > n:
                    conditions["condition_3"] += 1
                    return False, conditions

    sum_4 = sum(partition[i] * partition[j] for i in range(p - 1) for j in range(i + 1, p))
    if (sum_4 + m) > n:
        conditions["condition_4"] += 1
        return False, conditions
    
    if threshold is not None: 
        if exponent * (k_1 + k_2) > threshold:
            conditions["condition_5"] += 1
            return False, conditions

    return True, conditions

def find_min_time(algorithm_data):
    r"""
    Identifies the algorithm that requires the minimum runtime from the mq estimator evaluation. The algorithm name, time and memory are then outputted.
    """
    return min(
        ((algo, data['estimate']['time'], data['estimate']['memory'])
         for algo, data in algorithm_data.items()
         if isinstance(data['estimate']['time'], (int, float))),
        key=lambda x: x[1],
        default=(None, None, None)
    )

def mq_problem(partition, q, cache):
    r"""
    Builds a list of instances based on the current partition. 
    Concurrently, mq estimator is executed to derive estimates on the time and memory complexity. 
    Returns a list of the best-performing algorithm results, i.e. those with minimal time complexity, for each sub-instance.
    """
    try:
        results = []
        instances = [(partition[0] - partition[-1], partition[0] + partition[-2], q)] + [
            (partition[i], partition[i], q) for i in range(1, len(partition) - 2)
        ]
        for params in instances:
            if any(val <= 2 for val in params):
                results.append(("Undefined", 0, 0))
            else:
                if params not in cache:
                    cache[params] = MQEstimator(*params).estimate()
                results.append(find_min_time(cache[params]))
        return results
    except Exception:
        return None

def process_chunk(args):
    r"""
    Provides a mechanism for generating the partitions and evaluates them. 
    Iterates through possible values of x1, builds partitions using Cartesian products. 
    Applies filters i.e tests if the partition is valid and checks if constraints as defined in the paper are reflected. 
    Maintains two heaps for storing classical time and quantum time estimates. 
    Returns the best performing partition for classical and quantum estimations and statistic on failed conditions.
    """
    start_x1, end_x1, p, m, MAX_b1, MAX_SUMMAND, MIN_SUM_GUESS, n, threshold, exponent, q, cache, max_lines = args
    passed_results = []
    passed_results_quantum = []
    condition_counts = {"condition_1": 0, "condition_2": 0, "condition_3": 0, "condition_4": 0, "condition_5": 0}   
    heap = []
    heap_2 = []

    for x1 in range(start_x1, end_x1 + 1):
        if x1 > MAX_b1:
            continue
        for middle in itertools.product(range(1, MAX_SUMMAND + 1), repeat=p - 3):
            base = (x1,) + middle
            remaining = m - sum(base)
            for x_k_1 in range(1, remaining):
                x_k = remaining - x_k_1
                if x_k >= 1:
                    full = base + (x_k_1, x_k)
                    if not is_valid_partition(full, m, p, MAX_b1, MAX_SUMMAND, MIN_SUM_GUESS):
                        continue
                    passed, conds = filter_test(n, m, full, threshold, exponent)
                    if passed:
                        result = mq_problem(full, q, cache)
                        if result:
                            max_row = max(result, key=lambda row: row[1])
                            rounded_max = round(max_row[1], 1)
                            time_est = rounded_max + exponent * (x_k + x_k_1)
                            time_est_quantum = rounded_max + exponent * ((x_k + x_k_1) / 2)
                            memory = round(max_row[2], 1)
                            item = (time_est, (full, result, time_est, memory, time_est_quantum))
                            #print(item)
                            item_2 = (time_est_quantum, (full, result, time_est_quantum, memory, time_est))
                            
                            if len(heap) < max_lines:
                                heapq.heappush(heap, (-item[0], item[1]))
                            else:
                                if -heap[0][0] > item[0]:
                                    heapq.heappushpop(heap, (-item[0], item[1]))
                            
                            if len(heap_2) < max_lines:
                                heapq.heappush(heap_2, (-item_2[0], item_2[1]))
                            else:
                                if -heap_2[0][0] > item_2[0]:
                                    heapq.heappushpop(heap_2, (-item_2[0], item_2[1]))

                    else:
                        for k in conds:
                            condition_counts[k] += conds[k]
    
    passed_results = [item[1] for item in heap]
    passed_results_quantum = [item_2[1] for item_2 in heap_2]
    return passed_results, passed_results_quantum, condition_counts


def run_parallel(m, n, q, p, exponent, desc, threshold, cpus, MAX_SUMMAND, MAX_b1, MIN_SUM_GUESS, max_lines):
    r"""
    Runs the optimization algorithm in parallel, i.e. the workload is divided into ranges of x1 across CPUs, and the results are collected and merged. 
    The results should be written into three output files, with the classical and quantum and statistical results on conditions. The results are also displayed in the CLI.
    """
    print()
    print("#############################################################")
    print(f"Running optimization with m={m}, n={n}, q={q}, p={(p-2)}, exponent={exponent}, threshold={threshold}, maxsummand={MAX_SUMMAND}, maxb1={MAX_b1}, minsumguess={MIN_SUM_GUESS}, maxlines={max_lines}, cpus={cpus}")
    x1_range = list(range(1, MAX_b1 + 1))
    chunk_size = math.ceil(len(x1_range) / cpus)

    chunks = []
    for i in range(0, len(x1_range), chunk_size):
        chunk = x1_range[i:i+chunk_size]
        if not chunk:
            continue
        start_x1 = chunk[0]
        end_x1 = chunk[-1]
        chunks.append((start_x1, end_x1, p, m, MAX_b1, MAX_SUMMAND, MIN_SUM_GUESS, n, threshold, exponent, q, {}, max_lines))

    results = process_map(process_chunk, chunks, max_workers = cpus, chunksize = 1)

    heap = []
    heap_2 = []
    stats = {"condition_1": 0, "condition_2": 0, "condition_3": 0, "condition_4": 0, "condition_5": 0}

    for data, data_2, conds in results:
        for k in conds:
            stats[k] += conds[k]

        for entry in data:
            time_est = entry[2]
            if len(heap) < max_lines:
                heapq.heappush(heap, (-time_est, entry))
            else:
                if -heap[0][0] > time_est:
                    heapq.heappushpop(heap, (-time_est, entry))
        
        for entry in data_2:
            time_est_quantum = entry[2]
            if len(heap_2) < max_lines:
                heapq.heappush(heap_2, (-time_est_quantum, entry))
            else:
                if -heap_2[0][0] > time_est_quantum:
                    heapq.heappushpop(heap_2, (-time_est_quantum, entry))

    final_results = [entry for _, entry in sorted(heap, reverse=True)]

    final_results_quantum = [entry for _, entry in sorted(heap_2, reverse=True)]

    file_prefix = f"{desc}_m_{m}_n_{n}_q_{q}_p_{(p-2)}_exp_{exponent}_threshold_{threshold}_maxsummand_{MAX_SUMMAND}_maxb1_{MAX_b1}_minsumguess_{MIN_SUM_GUESS}_maxlines_{max_lines}"
    
    with open(f"{file_prefix}_C_results.txt", "w") as f:
        for entry in final_results:
            f.write(f"{entry}\n")
    
    with open(f"{file_prefix}_Q_results.txt", "w") as f:
        for entry in final_results_quantum:
            f.write(f"{entry}\n")

    with open(f"{file_prefix}_STATS.txt", "w") as f:
        for k, v in stats.items():
            f.write(f"{k}: {v} partitions failed\n")
    
    print("-------------------------------------------------------------")
    print("Classical:")
    if final_results:
        print(final_results[0])
    else:
        print("No classical solutions found!")
    print("-------------------------------------------------------------")
    print("Quantum:")
    if final_results_quantum:
        print(final_results_quantum[0])
    else:
        print("No quantum solutions found!")
    print("-------------------------------------------------------------")
    print(f"Optimization done. {len(final_results)} results written to output files.")
    print("#############################################################")
    print()

if __name__ == "__main__":
    r"""
    Responsible for parsing command-line arguments, configuring input parameters, and adjusting default values for optional parameters.
    """
    parser = argparse.ArgumentParser(description="Optimized MQ Search Pipeline")
    
    # Required parameter
    parser.add_argument("-n", type=int, required=True)
    parser.add_argument("-m", type=int, required=True)
    parser.add_argument("-q", type=int, required=True)
    parser.add_argument("-p", type=int, required=True)
    parser.add_argument("-c", "--cpus", type=int, required=True)

    # Optional parameter
    parser.add_argument("-desc", type=str, default="output")
    parser.add_argument("-t", "--threshold", type=int, default=None)  
    parser.add_argument("-maxsummand", type=int, default=None)
    parser.add_argument("-maxb1", type=int, default=None)
    parser.add_argument("-minsumguess", type=int, default=1)
    parser.add_argument("-lines", "--max_lines", type=int, default=1000)

    args = parser.parse_args()

    maxb1 = args.maxb1 if args.maxb1 is not None else args.m
    maxsummand = args.maxsummand if args.maxsummand is not None else args.m
    exp = round(math.log2(args.q), 1)
    p = args.p + 2 # (partitions + k1 + k2)

    run_parallel(args.m, args.n, args.q, p, exp, args.desc, args.threshold, args.cpus, maxsummand, maxb1, args.minsumguess, args.max_lines)
