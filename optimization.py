#!/usr/bin/env python3

import complexity_withpoly

import argparse
import math
from typing import Optional, Tuple, List

def feasible_for_pair(n: int, m: int, k: int, p: int) -> bool:
    """
    Check the two inequalities for given n,m,k,p.
    Inequality A (written to avoid floats): 2*n >= numerator
      numerator = ((m-k)*(m-k+1) - (m-k - (p-1))*(m-k - p))
    Inequality B: for all j in 1..Jmax where Jmax = m-k-(p-1),
      n >= (j-1)*(2*m - 2*k - p - (j-1)) + j
    If Jmax < 1, the second inequality is vacuously true.
    """
    m_k = m - k

    # Numerator of the first inequality's right-hand side times 2 (so we can compare integers)
    # numerator = (m_k*(m_k+1) - (m_k - (p-1))*(m_k - p))
    # numerator_1 = (m_k * (m_k + 1)) - ((m_k - p) * (m_k - p - 1)) + k
    numerator = p * (2*m - 2*k - p - 1)
    # print(numerator/2 + m) 
    # print(numerator_1/2)
    # print('######################')
    # print(numerator)
    # print(2*(n-m))
    if 2 * k + p - m < 0 or m - k - p < 0:
        return False
    if 2 * (n-m) < numerator:
        return False

    # Second inequality: compute its maximum RHS across j = 1..Jmax (if any)
    Jmax = m_k - p  # equals m-k-p+1
    if Jmax >= 2:
        max_rhs = -10**30
        for j in range(2, Jmax + 1):
            # term = (j - 1) * (2 * m - 2 * k - p - (j - 1)) + j
            term = (j - 1) * (2 * m - 2 * k - p - (j - 1)) + j-1
            if term > max_rhs:
                max_rhs = term
        if n - 1 < max_rhs:
            return False

    return True

def find_best_k_p(n: int, m: int, q: int, allow_p_zero: bool = True,polyfact = False) -> Optional[Tuple[int,int,int]]:
    """
    Search all k in [0..m] and p in [0..m-k] (or [1..m-k] if allow_p_zero=False).
    Returns (best_k, best_p, q**best_k) or None if no feasible pair found.
    """
    best_c = None  # will hold tuple (q_power, k, p)
    best_q = None
    p_start = 0 if allow_p_zero else 1

    for k in range(0, m + 1):
        for p in range(p_start, m - k + 1):
            if feasible_for_pair(n, m, k, p):
                # compute q^k (use pow to avoid float rounding; may be negative if q negative and k odd)
                cost = complexity_withpoly.total_cost(q, m, p, k, polynomial_factors = polyfact)
                # cost = [k * math.log2(q), (k/2) * math.log2(q)]
                candidate_c = (cost[0], k, p)
                candidate_q = (cost[1], k, p)
                if best_c is None:
                    best_c = candidate_c
                if best_q is None:
                    best_q = candidate_q
                else:
                    # choose smaller numeric q_power; tie-breaker: smaller k
                    if candidate_c[0] < best_c[0] or (candidate_c[0] == best_c[0] and candidate_c[1] < best_c[1]):
                        best_c = candidate_c
                    if candidate_q[0] < best_q[0] or (candidate_q[0] == best_q[0] and candidate_q[1] < best_q[1]):
                        best_q = candidate_q

    if best_c is None:
        return None
    else:
        return (best_c[1], best_c[2], best_c[0]), (best_q[1], best_q[2], best_q[0])

def all_feasible_pairs(n: int, m: int, q: int, allow_p_zero: bool = True) -> List[Tuple[int,int,int]]:
    """Return list of all feasible (k,p,q**k) tuples (unsorted)."""
    res = []
    p_start = 0 if allow_p_zero else 1
    for k in range(0, m + 1):
        for p in range(p_start, m - k + 1):
            if feasible_for_pair(n, m, k, p):
                res.append((k, p, pow(q, k)))
    return res

def main():
    parser = argparse.ArgumentParser(description="Find k and p that satisfy inequalities and minimize q^k")
    parser.add_argument("n", type=int, help="integer n")
    parser.add_argument("m", type=int, help="integer m (>=0)")
    parser.add_argument("q", type=int, help="integer q")
    parser.add_argument("--no-zero-p", dest="p_zero", action="store_false",
                        help="disallow p=0 (i.e. use p in 1..m-k instead of 0..m-k)")
    args = parser.parse_args()

    n, m, q = args.n, args.m, args.q
    res = find_best_k_p(n, m, q, allow_p_zero=args.p_zero,polyfact = False)
    res_poly = find_best_k_p(n, m, q, allow_p_zero=args.p_zero, polyfact = True)

    if res_poly is None:
        print("No feasible (k,p) pair found for the given n,m,q with current p-range.")

    best_c_poly, best_q_poly = res_poly
    best_c, best_q = res
    overhead_c = complexity_withpoly.exp_dfs(best_c[1], q)
    overhead_q = complexity_withpoly.exp_dfs(best_q[1], q)
    overhead_c_poly = complexity_withpoly.exp_dfs(best_c_poly[1], q)
    overhead_q_poly = complexity_withpoly.exp_dfs(best_q_poly[1], q)
    print(f'n = {n}')
    print(f'm = {m}')
    print(f'q = {q}')
    print(f"Best pair found Classical: k = {best_c[0]}, p = {best_c[1]}")
    print(f"Best pair found Classical with poly: k = {best_c_poly[0]}, p = {best_c_poly[1]}")
    print(f'Probability that MQ tree has p levels: {overhead_c[1]}')
    print(f'Number of nodes visited by the DFS: {overhead_c[0]}')
    print(f'Final Classical Complexity: {best_c[2]}')
    print(f'Final Classical Complexity with poly: {best_c_poly[2]}')
    print('#################################################')
    print('#################################################')
    print(f"Best pair found Quantum: k = {best_q[0]}, p = {best_q[1]}")
    print(f"Best pair found Quantum with poly: k = {best_q_poly[0]}, p = {best_q_poly[1]}")
    print(f'Probability that MQ tree has p levels: {overhead_q[1]}')
    print(f'Number of nodes visited by the DFS: {overhead_q[0]}')
    print(f'Final Quantum Complexity: {best_q[2]}')
    print(f'Final Quantum Complexity with poly: {best_q_poly[2]}')

if __name__ == "__main__":
    main()
