import argparse
import math

from fractions import Fraction
from math import log1p, exp
import numpy as np
from numba import njit


@njit
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


@njit
def exp_dfs(t, q):
    """
    Evaluate the depth-first search (DFS) success probability and expected
    number of solves up to depth t in the MQ tree.

    Args:
        t: Maximum search depth (number of levels).
        q: Branching factor at each node.

    Returns:
        Tuple (p_t, exp_solves_minus_one):
            p_t: Probability of success at depth t.
            exp_solves_minus_one: Expected number of solves minus one.
    """
    s = np.zeros(t + 1, np.float64)

    s[0] = 1.0

    P0, P1, P2, Pq = _compute_transition_probabilities(q)

    for depth in range(1, t + 1):
        prev_s = s[depth - 1]

        # Update success probability
        s[depth] = 1 - (
            P0
            + P1 * (1.0 - prev_s)
            + P2 * (1.0 - prev_s)**2
            + Pq * (1.0 - prev_s)**q
        )

    E = 1 + s[t] * sum([1/s[i] for i in range(t)])

    # return E, s[t]
    return E, s[t], sum([1/s[i] for i in range(t)])


#          

def log_N_over_qsqr(s, q, r):
    """Return ln( N(s,r) / q**{s**2} )."""
    logN = 0.0
    for j in range(r):
        logN += 2.0 * math.log(q**s - q**j) - math.log(q**r - q**j)
    logQsqr = (s * s) * math.log(q)
    return logN - logQsqr

def prob_power_one_minus(q, k, i, s, precise = False):
    """Return 1 - (1 - q**{-k})**{q**i} computed stably as a float."""
    if k <= 0:
        return 1.0 if i > 0 else 0.0
    if precise:
        ln = math.log1p(- q**(-k))          # ln(1 - q**{-k})
    else:
        ln = math.log1p(- s*q**(-k))          # ln(1 - q**{-k})
    exponent = (q ** i) * ln           # q**i * ln(1 - q**{-k})
    inner_minus1 = math.expm1(exponent)
    factor = -inner_minus1
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
    '''
    Return [# of field multiplication, # of field additions] to assign k variables to n_lin equations. If check = True, this is in the check phase, so no inner product to be performed.
    '''
    assign_assign = [cost *n_lin for cost in cost_assign_assign(k_assigned)]
    if check:
        assign_var = 0
        return [assign_assign[0] + assign_var, assign_assign[1] + assign_var]
    else:
        assign_var = [cost * n_lin for cost in cost_assign_var(k_assigned)]
        return [assign_assign[0] + assign_var[0], assign_assign[1] + assign_var[1]]

def total_cost(q,m,p,k, polynomial_factors = False):
    e_p, s_p, _ = exp_dfs(p, q)
    expcost_trivialMQstep = 0
    if polynomial_factors:
        for i in range(p):
            expcost_trivialMQstep += sum(cost_sub_assign(k+i,1)) * (exp_dfs(i+1,q)[0] - exp_dfs(i,q)[0])
        expcost_trivialMQstep += e_p
        expcost_linsolve = 2/3*(m - k - p)**3
        expcost_trivialMQsteplin = (m - k - p) * sum(cost_sub_assign(k + p, m-k-p))
        cost_check = k * sum(cost_sub_assign(m, 1, check = True))
        P = s_p*q**(-k)
        exp_checks = 1
    else:
        expcost_trivialMQstep = e_p
        cost_check = 0
        expcost_trivialMQsteplin = 0
        expcost_linsolve = 0
        P = s_p*q**(-k)
        exp_checks = 1
    cost = expcost_trivialMQstep + s_p * (expcost_trivialMQsteplin + expcost_linsolve + exp_checks * cost_check)
    # iterations_q = math.sqrt(prob_power_one_minus(q, k, k, s_p)/(P * s_p))
    # iterations_c = iterations_q**2
    iterations_c = q**(k)/s_p
    iterations_q = math.sqrt(iterations_c)
    return  math.log2(iterations_c * cost), math.log2(iterations_q * cost)

def winning_probability(q,m,k,p, precise = False):
    e_p, s_p, _ = exp_dfs(p,q)
    if precise:
        pass
    else:
        return prob_power_one_minus(q, k, k, s_p), s_p, e_p



if __name__ == "__main__":
    """Parameter description"""
    parser = argparse.ArgumentParser(description="DFS Data")
    parser.add_argument("-m", type=int, required=True)
    parser.add_argument("-p", type=int, required=True)
    parser.add_argument("-q", type=int, required=True)
    parser.add_argument("-k", type=int, help="Guessed", required=True)
    parser.add_argument("--polyfactors", action="store_true", help="Add polynomial factors.")

    args = parser.parse_args()
    result = exp_dfs(args.p, args.q)
    print(f'Success Probability MQ-DFS: {result[1]}')
    print(f'Expected solves per MQ tree: {result[0]}')
    complexity = total_cost(args.q, args.m, args.p, args.k, polynomial_factors = args.polyfactors)
    # print(winning_probability(args.q, args.m, args.k, args.p))
    print(f'Final Classical Complexity: {complexity[0]}')
    print(f'Final Quantum Complexity: {complexity[1]}')
