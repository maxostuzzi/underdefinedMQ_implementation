#!/usr/bin/env python3

import compute_complexity

import argparse
import math
from prettytable import PrettyTable
from typing import Optional, Tuple, List

def feasible_for_pair(n: int, m: int, k: int, p: int) -> bool:

    dimension_deficiency = 0
    max_def = 0
    if 2 * k + p - m < 0: # structurally not allowed
        return False
    if m - k - p < 0:
        return False

    numerator = p * (2*m - 2*k - p - 1)
    if 2 * (n-m) < numerator or n - 1 < (m - k - p - 1) * (m - k + 2):
        return False

    return True

def find_best_k_p(n: int, m: int, q: int, allow_p_zero: bool = True, polyfact = False, bitops = True) -> Optional[Tuple[int,int,int]]:
    """
    Search all k in [0..m] and p in [0..m-k] (or [1..m-k] if allow_p_zero=False).
    Returns (best_k, best_p, q**best_k) or None if no feasible pair found.
    """
    best_c = None  # will hold tuple (q_power, k, p)
    best_q = None
    p_start = 0 if allow_p_zero else 1

    for k in range(0, m + 1):
        for p in range(p_start, m - k + 1):
            feasible = feasible_for_pair(n, m, k, p)
            if feasible:
                solving_phase_cost = compute_complexity.complexity(q, m, p, k, polynomial_factors = polyfact, bitops = bitops)
                if polyfact:
                    gaussian_elimination = 2/3 * n**3
                else:
                    gaussian_elimination = 1
                cost = {'classical':{'transformation': gaussian_elimination, 'solving': 2 ** solving_phase_cost[0]}, 'quantum':{'transformation': gaussian_elimination, 'solving': 2 ** solving_phase_cost[1]}}
                # cost = [math.log2(2 ** solving_phase_cost[0] + q ** dimension_deficiency), math.log2(2 ** solving_phase_cost[1] + q ** (dimension_deficiency/2))]
                candidate_c = {'cost': {'total': math.log2(sum(cost['classical'].values())), 'transformation': math.log2(cost['classical']['transformation']), 'solving': math.log2(cost['classical']['solving'])}, 'k': k, 'p': p}
                candidate_q = {'cost': {'total': math.log2(sum(cost['quantum'].values())), 'transformation': math.log2(cost['quantum']['transformation']), 'solving': math.log2(cost['quantum']['solving'])}, 'k': k, 'p':p}
                if best_c is None:
                    best_c = candidate_c
                if best_q is None:
                    best_q = candidate_q
                else:
                    if candidate_c['cost']['total'] < best_c['cost']['total'] or (candidate_c['cost']['total'] == best_c['cost']['total'] and candidate_c['k'] > best_c['k']):
                        best_c = candidate_c
                    if candidate_q['cost']['total'] < best_q['cost']['total'] or (candidate_q['cost']['total'] == best_q['cost']['total'] and candidate_q['k'] > best_q['k']):
                        best_q = candidate_q

    if best_c is None:
        return None
    else:
        return {'classical': best_c, 'quantum': best_q}

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
    parser.add_argument("--quantum", dest="quantum", action="store_true")
    parser.add_argument("--all", dest="all", action="store_true")
    args = parser.parse_args()

    n, m, q, quantum, both = args.n, args.m, args.q, args.quantum, args.all
    res = find_best_k_p(n, m, q, polyfact = False)
    res_poly = find_best_k_p(n, m, q, polyfact = True)

    if res_poly is None:
        print("No feasible (k,p) pair found for the given n,m,q with current p-range.")

    k_poly, p_poly = res_poly['classical']['k'], res_poly['classical']['p']
    k, p = res['classical']['k'], res['classical']['p']
    k_poly_q, p_poly_q = res_poly['quantum']['k'], res_poly['quantum']['p']
    k_q, p_q = res['quantum']['k'], res['quantum']['p']
    # overhead_c_poly = compute_complexity.overhead(q, m, best_c_poly[1], best_c_poly[0])
    # overhead_q_poly = compute_complexity.overhead(q, m, best_q_poly[1], best_q_poly[0])
    table_classical = PrettyTable()
    table_classical.field_names = ['\x1b[32mCLASSICAL\x1b[0m', 'Complexity', 'Complexity w/ Poly']
    table_classical.add_rows(
        [
        ['Transformation Phase', round(res['classical']['cost']['transformation'], 2), round(res_poly['classical']['cost']['transformation'], 2)],
        ['Just Guess Phase', round(res['classical']['cost']['solving'], 2), round(res_poly['classical']['cost']['solving'], 2)],
        ['Total', round(res['classical']['cost']['total'], 2), round(res_poly['classical']['cost']['total'], 2)]
        ])
    table_quantum = PrettyTable()
    table_quantum.field_names = ['\x1b[32mQUANTUM\x1b[0m', 'Complexity', 'Complexity w/ Poly']
    table_quantum.add_rows(
        [
        ['Transformation Phase', round(res['quantum']['cost']['transformation'], 2), round(res_poly['quantum']['cost']['transformation'], 2)],
        ['Just Guess Phase', round(res['quantum']['cost']['solving'], 2), round(res_poly['quantum']['cost']['solving'], 2)],
        ['Total', round(res['quantum']['cost']['total'], 2), round(res_poly['quantum']['cost']['total'], 2)]
        ])
    print(f'n = {n}')
    print(f'm = {m}')
    print(f'q = {q}')
    print('')
    if both:
        print(f"Best pair without poly: k = {k}, p = {p}")
        print(f"Best pair with poly: k = {k_poly}, p = {p_poly}")
        print(table_classical)
        print('')
        print(f"Best pair without poly: k = {k_q}, p = {p_q}")
        print(f"Best pair with poly: k = {k_poly_q}, p = {p_poly_q}")
        print(table_quantum)
    elif quantum:
        print(f"Best pair without poly: k = {k_q}, p = {p_q}")
        print(f"Best pair with poly: k = {k_poly_q}, p = {p_poly_q}")
        print(table_quantum)
    elif not quantum and not both:
        print(f"Best pair without poly: k = {k}, p = {p}")
        print(f"Best pair with poly: k = {k_poly}, p = {p_poly}")
        print(table_classical)
















    # print(f'\x1b[32mCLASSICAL\x1b[0m')
    # print(f"Best pair found: k = {k}, p = {p}")
    # print(f'Dimension deficiency: {res['classical']['dim_deficiency']}')
    # print(f'Transformation Phase Complexity: {round(res['classical']['cost']['transformation'], 2)}')
    # print(f'Just Guess Phase Complexity: {round(res['classical']['cost']['solving'], 2)}')
    # print(f'Total Complexity: {round(res['classical']['cost']['total'], 2)}')
    # print('')
    # print(f'\x1b[32mCLASSICAL WITH POLYFACTORS\x1b[0m')
    # print(f"Best pair found: k = {k_poly}, p = {p_poly}")
    # print(f'Dimension deficiency: {res_poly['classical']['dim_deficiency']}')
    # print(f'Transformation Phase Complexity: {round(res_poly['classical']['cost']['transformation'], 2)}')
    # print(f'Just Guess Phase Complexity: {round(res_poly['classical']['cost']['solving'], 2)}')
    # print(f'Total Complexity: \x1b[32m{round(res_poly['classical']['cost']['total'], 2)}\x1b[0m')
    # print('')
    # print('\x1b[32mQUANTUM\x1b[0m')
    # print(f"Best pair found: k = {k_q}, p = {p_q}")
    # print(f'Dimension deficiency: {res['quantum']['dim_deficiency']}')
    # print(f'Transformation Phase Complexity: {round(res['quantum']['cost']['transformation'], 2)}')
    # print(f'Just Guess Phase Complexity: {round(res['quantum']['cost']['solving'], 2)}')
    # print(f'Total Complexity: {round(res['quantum']['cost']['total'], 2)}')
    # print('')
    # print('\x1b[32mQUANTUM WITH POLYFACTORS\x1b[0m')
    # print(f"Best pair found: k = {k_poly_q}, p = {p_poly_q}")
    # print(f'Dimension deficiency: {res_poly['quantum']['dim_deficiency']}')
    # print(f'Transformation Phase Complexity: {round(res_poly['quantum']['cost']['transformation'], 2)}')
    # print(f'Just Guess Phase Complexity: {round(res_poly['quantum']['cost']['solving'], 2)}')
    # print(f'Total Complexity: \x1b[32m{round(res_poly['quantum']['cost']['total'], 2)}\x1b[0m')

if __name__ == "__main__":
    main()
