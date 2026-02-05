import argparse
import math

from fractions import Fraction
from math import log1p, exp
import numpy as np

def _compute_transition_probabilities(q):
    """
    Compute the transition probabilities P0, P1, P2, and Pq for the DFS process.

    Args:
        q: Branching factor (number of children per node).

    Returns:
        Tuple of probabilities (P0, P1, P2, Pq).
    """
    P0 = (q**3 - 2.0 * q**2 + 3.0 * q - 2.0) / (2.0 * q**3)
    P1 = 2.0 * (q - 1.0) / (q**2)
    P2 = (q - 1.0)**2 / (2.0 * q**2)
    Pq = 1.0 / (q**3)
    return P0, P1, P2, Pq

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
    # dimension_deficiency = 0
    # for i in range(2, m - k - p + 1):
    #     m_i = 2 * m - 2 * k - p - i + 1
    #     print(f'm_i = {m_i}')
    #     if n - i - m_i < 0:
    #         dimension_deficiency += m_i + i - n
    #         print(m_i + i - n)
    #         print(i)
    # print('End First Part')
    # for i in range(m - k - p + 1, m - k + 1):
    #     m_i = sum(j for j in range(m - k - p, i - 1)) + (i - 1) * (m - k - i + 1)
    #     print(f'm_i = {m_i}')
    #     if n - i - m_i < 0:
    #         dimension_deficiency += m_i + i - n
    #         print(m_i + i - n)
    #         print(i)
    # for i in range(m - k + 1, m + 1):
    #     m_i = p * (m - k) - sum(j for j in range(1, p+1))
    #     print(f'm_i = {m_i}')
    #     if m_i > n - i:
    #         dimension_deficiency += m_i + i - n
    #         print(m_i + i - n)
    #         print(i)
    # print(f'Dim-Def: {dimension_deficiency}')
    if 2 * k + p - m < 0 or m - k - p < 0:
        print('STRUCTURALLY NOT POSSIBLE')
        return False
    if 2 * (n-m) < numerator or n - 1 < (m - k - p - 1) * (m - k + 2):
        print(f'Required dimension: {max((p * (2*m - 2*k - p - 1))/2 + m, (m - k - p - 1) * (m - k + 2) + 1)}')
        return False

    

    return True

def prob_power_one_minus(q, k, i, s, precise = False):
    """Return 1 - (1 - q**{-k})**{q**i} computed stably as a float."""
    if k <= 0:
        return 1.0 if i > 0 else 0.0
    if precise:
        ln = math.log1p(- q**(-k))          # ln(1 - q**{-k})
    else:
        ln = math.log1p(- s*q**(-k))          # ln(1 - q**{-k})
    exponent = (q ** i) * ln           # q**i * ln(1 - q**{-k})
    # use expm1 for accuracy: exp(exponent) - 1, then take negative
    inner_minus1 = math.expm1(exponent)
    factor = -inner_minus1
    # clamp for safety
    if factor < 0.0:
        factor = 0.0
    elif factor > 1.0:
        factor = 1.0
    return factor

def cost_assign_assign(n):
    '''
    Returns [# field mult, # field add] for matrix vector multiplication u**tAu
    '''
    return [n * (n+1)/2 + n, n * (n+1)/2 + n - 1]

def cost_assign_var(n):
    '''
    Returns [# of field multiplications, # of field additions] to perform an inner-product between two n-dim vectors.
    '''
    return [n,n-1]

def cost_sub_assign(k_assigned, n_lin, check = False):
    assign_assign = [cost *n_lin for cost in cost_assign_assign(k_assigned)]
    if check:
        assign_var = 0
        return [assign_assign[0] + assign_var, assign_assign[1] + assign_var]
    else:
        assign_var = [cost * n_lin for cost in cost_assign_var(k_assigned)]
        return [assign_assign[0] + assign_var[0], assign_assign[1] + assign_var[1]]

def overhead(q,m,p,k, poly_factors = True, bitop = True):
    C_MQtree = p
    if poly_factors:
        for i in range(p):
            C_MQtree += sum(cost_sub_assign(k+i,1))
        C_linsub = (m - k - p) * sum(cost_sub_assign(k + p, m-k-p))
        C_linsolve = 2/3*(m - k - p)**3
        C_check = sum(cost_sub_assign(m, 1, check = True)) * q * (1 - q**(-k))/(q-1)
        cost = C_MQtree + C_linsub + C_linsolve + C_check
        if bitop:
            cost *= math.log2(q)
    else:
        cost = C_MQtree
    return cost

def complexity(q,m,p,k, polynomial_factors = True, bitops = True):
    iterations_c = q**k
    iterations_q = q**(k/2)
    cost = overhead(q,m,p,k, poly_factors = polynomial_factors, bitop = bitops)
    if p!=0 or polynomial_factors:
        return  math.log2(iterations_c * cost), math.log2(iterations_q * cost)
    else:
        return  math.log2(iterations_c), math.log2(iterations_q)



if __name__ == "__main__":
    """Parameter description"""
    parser = argparse.ArgumentParser(description="DFS Data")
    parser.add_argument("-n", type=int, required=True)
    parser.add_argument("-m", type=int, required=True)
    parser.add_argument("-p", type=int, required=True)
    parser.add_argument("-q", type=int, required=True)
    parser.add_argument("-k", type=int, help="Guessed", required=True)
    parser.add_argument("--polyfactors", action="store_true", help="Add polynomial factors.")

    args = parser.parse_args()
    feasible_for_pair(args.n, args.m, args.k, args.p)
    overh = overhead(args.q, args.m, args.p, args.k)
    complexity = complexity(args.q, args.m, args.p, args.k, polynomial_factors = args.polyfactors)
    # print(winning_probability(args.q, args.m, args.k, args.p))
    print(f'Final Classical Complexity: {complexity[0]}')
    print(f'Final Quantum Complexity: {complexity[1]}')
