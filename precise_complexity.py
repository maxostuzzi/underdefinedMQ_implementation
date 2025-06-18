import argparse
import math

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
    number of solves up to depth t in a q-ary tree.

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
    # assume s is already filled to length t+1
    # x = [None] + [s[r]/s[r-1] for r in range(1, t+1)]
    # E = 0.0
    # for k in range(0, t+1):
    #     prod = 1.0
    #     for r in range(k+1, t+1):
    #         prod *= x[r]
    #     E += prod

    E = 1 + s[t] * sum([1/s[i] for i in range(t)])

    # return E, s[t]
    return E, s[t], sum([1/s[i] for i in range(t)])

def main():
    parser = argparse.ArgumentParser(
        description="Compute DFS success probability and complexity metrics."
    )
    parser.add_argument(
        "-p", type=int, required=True,
        help="Maximum search depth (integer)."
    )
    parser.add_argument(
        "-q", type=int, required=True,
        help="Branching factor per node (integer)."
    )
    parser.add_argument(
        "-k", nargs='+', type=int, default=[],
        help="List of guessed solution depths for complexity calculation."
    )
    parser.add_argument(
        "--grover", action="store_true",
        help="Compute quantum (Grover) complexity if set; otherwise classical."
    )

    args = parser.parse_args()
    exp_solves, new_prob, factor = exp_dfs(args.p, args.q)

    print(f"Probability of success at depth {args.p}: {new_prob:.6f}")
    print(f"Expected number of solves: {exp_solves:.6f}")

    if args.k:
        total_k = sum(args.k)
        if args.grover:
            complexity = math.log(
                args.q**(total_k / 2) * math.sqrt(exp_solves) / new_prob,
                2
            )
            print(f"Final quantum complexity (log₂): {complexity:.6f} = {math.log(exp_solves, 2):.6f} + {math.log(args.q**total_k /new_prob,2)/2:.6f}")
        else:
            complexity = math.log(
                args.q**total_k * exp_solves / new_prob,
                2
            )
            print(f"Final classical complexity (log₂): {complexity:.6f} = {math.log(exp_solves, 2):.6f} + {math.log(args.q**total_k /new_prob,2):.6f}")


if __name__ == "__main__":
    main()
