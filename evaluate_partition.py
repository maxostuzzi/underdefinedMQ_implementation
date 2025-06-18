import argparse
from cryptographic_estimators.MQEstimator import MQEstimator
import itertools
import math
from tqdm.contrib.concurrent import process_map
import heapq

def filter_test(n, m, partition, exponent):
    """Applies a filter mechanism"""
    k_1, k_2 = partition[-2], partition[-1]
    p = len(partition) - 2
    conditions = {"condition_1": 0, "condition_2": 0, "condition_3": 0, "condition_4": 0, "condition_5": 0}

    for i in range(1, partition[0]):
        if ((i * (sum(partition[:-2]) + k_2 - i)) + i) > n:
            conditions["condition_1"] += 1
            print(f'{((i * (sum(partition[:-2]) + k_2 - i)) + i)} > {n}')
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
        print(f'{sum_4 + m} > {n}')
        return False, conditions

    return True, conditions

def find_min_time(algorithm_data):
    """Finds the minimum runtime of the mq-solver"""
    return min(
        ((algo, data['estimate']['time'], data['estimate']['memory'])
         for algo, data in algorithm_data.items()
         if isinstance(data['estimate']['time'], (int, float))),
        key=lambda x: x[1],
        default=(None, None, None)
    )

def mq_problem(partition, q):
    """Runs a mq-solver on a selected partition"""
    try:
        results = []
        instances = [(partition[0] - partition[-1], partition[0] + partition[-2], q)] + [(partition[i], partition[i], q) for i in range(1, len(partition) - 2)] 
        for params in instances:
            if any(val <= 2 for val in params):
                results.append(("Undefined", 0, 0))
            else:
                results.append(find_min_time(MQEstimator(*params).estimate()))
        return results
    except Exception:
        return None


def eval_partition(n, m, q, partition, k):
    exponent = round(math.log2(q), 1)
    passed = filter_test(n, m, partition + k, round(math.log2(q), 1))
    if passed[0]:
        result = mq_problem(partition + k, q)
    else:
        print('Not enough variables.')
        print(passed)
        exit()
    if result:
        max_row = max(result, key=lambda row: row[1])
        rounded_max = round(max_row[1], 1)
        time_est = rounded_max + exponent * sum(k)
        time_est_quantum = rounded_max + exponent * (sum(k) / 2)
        memory = round(max_row[2], 1)
        item = (partition+k, result, time_est, memory, time_est_quantum)
        item_2 = (partition+k, result, time_est_quantum, memory, time_est)
    return item, item_2

if __name__ == "__main__":
    """Parameter description"""
    parser = argparse.ArgumentParser(description="Optimized MQ Search Pipeline")
    parser.add_argument("-n", type=int, required=True)
    parser.add_argument("-m", type=int, required=True)
    parser.add_argument("-q", type=int, required=True)
    parser.add_argument("--grover", action="store_true", help="Compute quantum (Grover) complexity if set; otherwise classical.")
    parser.add_argument("-partition", nargs='+', type=int, help="Partition")
    parser.add_argument("-k", nargs='+', type=int, help="Guessed")

    args = parser.parse_args()
    exp = round(math.log2(args.q), 1)
    assert sum(args.partition + args.k) == args.m
    final_results = eval_partition(args.n, args.m, args.q, args.partition, args.k)
    if args.grover:
        print(f'Quantum: {final_results[1][-3]}')
    else:
        print(f'Classical: {final_results[0][-3]}')
