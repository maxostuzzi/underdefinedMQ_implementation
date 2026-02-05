
#
# Conventions used here:
#   - P_list is a list of m matrices, each n x n, representing quadratic polynomials
#     with the paper's convention p(x) = sum_{i <= j} P[i,j] * x_i * x_j
#   - variables indices are 0-based in the code
#   - the first m variables are x_0 .. x_{m-1}, the "extra" variables are x_m .. x_{n-1}
#
# Returns a Sage vector y (length n) which is a solution for P(y) = t, or None if not found.

load('just_guess_transformation.sage')
from itertools import product
import argparse

def _extract_top_left(Pmat, m):
    """Return the top-left m x m submatrix of Pmat (Pmat is n x n)."""
    return Matrix(Pmat.rows()[:m]).matrix_from_rows_and_columns(range(m), range(m))

def _evaluate_quadratic_value(Pmat_m, x_values):
    """
    Evaluate polynomial p(x) = sum_{i<=j} P[i,j]*x_i*x_j where x_values is a dict index->value
    Only indices present in x_values are used; missing indices treated as 0.
    """
    FF = list(x_values.values())[0].parent() if x_values else None
    const = FF(0) if FF is not None else 0
    m = Pmat_m.nrows()
    for i in range(m):
        xi = x_values.get(i, None)
        for j in range(i, m):
            coef = Pmat_m[i, j]
            if coef == 0:
                continue
            xj = x_values.get(j, None)
            if xi is not None and xj is not None:
                const += coef * xi * xj
            elif xi is not None and xj is None:
                # xj assumed zero => contributes 0
                pass
            elif xi is None and xj is not None:
                pass
            else:
                pass
    return const

def _build_univariate_poly_in_var(Pmat_m, unknown_idx, assigned_vars, FF):
    """
    Build polynomial in z over FF for the equation p(x)=rhs when all variables except
    unknown_idx are substituted by assigned_vars (dictionary index->FF element).
    Returns a Sage polynomial in variable z and a flag indicating if other unknowns
    appear (in which case we should not proceed).
    """
    R = PolynomialRing(FF, 'z')
    z = R.gen()
    poly = R(0)
    m = Pmat_m.nrows()
    # check for dependencies on other unknown variables
    for i in range(m):
        for j in range(i, m):
            coef = Pmat_m[i, j]
            if coef == 0:
                continue
            is_i_known = (i in assigned_vars)
            is_j_known = (j in assigned_vars)
            if i == unknown_idx and j == unknown_idx:
                # coef * z^2
                poly += R(coef) * z**2
            elif i == unknown_idx and is_j_known:
                poly += R(coef) * z * R(assigned_vars[j])
            elif j == unknown_idx and is_i_known:
                poly += R(coef) * z * R(assigned_vars[i])
            elif (i == unknown_idx) and (not is_j_known):
                # other unknown present
                return None, True
            elif (j == unknown_idx) and (not is_i_known):
                return None, True
            elif (not is_i_known) and (not is_j_known):
                # both unknown and not our variable -> degree>1 in other unknowns
                return None, True
            else:
                # both known -> constant
                poly += R(coef) * R(assigned_vars[i]) * R(assigned_vars[j])
    return poly, False

def _roots_univariate(poly, FF):
    """Return list of roots in FF of polynomial poly (poly is a polynomial over FF)."""
    # If poly is zero polynomial -> every element is a root
    if poly == 0:
        return list(FF)  # all elements
    # Use Sage's roots method on the polynomial ring; convert roots to field elements
    roots_mult = poly.roots(multiplicities=False)
    # .roots returns elements of FF; ensure unique
    roots = []
    for r in roots_mult:
        if r not in roots:
            roots.append(r)
    return roots

def _solve_linear_system_all_solutions(A, b):
    """
    Solve A x = b over finite field. Return a generator of all solutions as vectors.
    A : Matrix over FF (r x n)
    b : vector over FF (length r)
    Yields each solution as a vector of length n.
    """
    FF = A.base_ring()
    r, n = A.nrows(), A.ncols()
    # check consistency using ranks:
    rankA = A.rank()
    Ab = A.augment(b.column())
    rankAb = Ab.rank()
    if rankA != rankAb:
        # inconsistent
        return []
    # find a particular solution (if A has full column rank this is unique)
    # Sage's solve_right works if system consistent:
    try:
        x_part = A.solve_right(b)
    except Exception:
        # fallback: gaussian elimination on augmented matrix
        # Solve using basic linear algebra: compute row-reduced echelon and then parametrize kernel
        # For simplicity: use right_kernel and find a particular by solving using a subset of columns
        x_part = zero_vector(FF, n)
        # attempt to pivot on independent columns
        pivcols = A.pivots()
        # set non-pivots to zero, compute pivots solution
        # Form submatrix for pivot columns:
        if len(pivcols) > 0:
            A_piv = A.matrix_from_columns(pivcols)
            try:
                sol_piv = A_piv.solve_right(b)
                for i, col in enumerate(pivcols):
                    x_part[col] = sol_piv[i]
            except Exception:
                # best-effort: use least-squares-ish approach; but this situation is rarely needed
                raise
    K = A.right_kernel()  # vector space basis for kernel as list of vectors dimension d
    basis = K.basis()
    d = len(basis)
    if d == 0:
        return [x_part]  # unique
    # enumerate all combinations of kernel basis coefficients over FF
    sols = []
    # iterate over q^d combinations
    if d >= 6:
        # guard: enumerating q^d may be huge; still we follow the algorithm's behavior.
        pass
    # produce all combinations in cartesian product
    for coeffs in product(FF, repeat=d):
        vec = vector(FF, list(x_part))
        for c, bv in zip(coeffs, basis):
            vec += c * bv
        sols.append(vec)
    return sols

def just_guess_solver(FF, P_list, t_vec, k, p, transform):
    """
    Main Just-Guess solver.
    - FF: finite field (Sage GF(q) or GF(p^r))
    - P_list: list of m matrices, each n x n (upper-triangular representation)
    - t_vec: vector of length m over FF (target)
    - k: number of guessed coordinates (0 <= k <= m)
    - p: number of MQ(1,1) equations to solve (0 <= p <= m-k)
    - transformation: function(FF, P_list, k, p) -> [T (m x m), S (n x n), Ptilde_list]
    Returns solution y (length n vector over FF) or None if not found.
    """
    m = len(P_list)
    n = P_list[0].nrows()
    # get transformation
    T, S, Ptilde_full = transform
    # transformed target
    t_tilde = T * t_vec

    # restrict to top-left m x m submatrices because variables x_{m}..x_{n-1} are fixed to zero
    Ptilde = [_extract_top_left(M, m) for M in Ptilde_full]

    # Predefine last n-m variables = 0
    assigned_fixed_rest = {i: FF(0) for i in range(m, n)}

    # guessed variable indices (we follow the paper: the guessed k coordinates are
    # the last k among the first m variables: indices m-k .. m-1)
    guess_indices = list(range(m - k, m))

    # indices of MQ(1,1) equations to solve: we start at Ptilde[m-1], then Ptilde[m-2], ...
    # unknowns solved by these are indices: m-k-1, m-k-2, ...
    print('Guessing...')
    def process_guess_assignment(guess_values):
        # guess_values: tuple of length k with elements in FF, corresponds to guess_indices
        assigned = dict(assigned_fixed_rest)  # start with fixed zeros for x_m..x_{n-1}
        for idx, val in zip(guess_indices, guess_values):
            assigned[idx] = FF(val)
        # Now do MQ tree step solving p single-variable equations sequentially
        solved_order = []  # list of (eq_idx, solved_var_idx, value)
        # recursion backtracking: we need to manage branching when MQ(1,1) has multiple roots
        solutions_after_mqtree = []

        def dfs_mqtree(j, assigned_local):
            """
            j: how many of the p MQ(1,1) we've already solved (0..p)
            assigned_local: dict of current assignments (indices -> FF elements)
            """
            if j == p:
                # Completed MQ tree step; record assigned_local as candidate
                solutions_after_mqtree.append(dict(assigned_local))
                return

            eq_idx = m - 1 - j           # equation index to solve (0-based)
            unknown_idx = m - k - 1 - j  # variable index that should be solved by this equation

            # Build univariate polynomial in z
            poly, has_other_unknowns = _build_univariate_poly_in_var(Ptilde[eq_idx], unknown_idx, assigned_local, FF)
            if has_other_unknowns:
                # transformation was supposed to prevent this; abandon this branch
                return

            # But polynomial includes constant terms; we need polynomial - t_tilde[eq_idx] = 0
            R = PolynomialRing(FF, 'z')
            z = R.gen()
            if poly is None:
                return
            poly_minus_rhs = poly - R(t_tilde[eq_idx])

            # get roots
            roots = _roots_univariate(poly_minus_rhs, FF)
            # if no roots, dead branch
            if len(roots) == 0:
                return
            # branch on all roots
            for r in roots:
                new_assigned = dict(assigned_local)
                new_assigned[unknown_idx] = FF(r)
                dfs_mqtree(j + 1, new_assigned)

        dfs_mqtree(0, assigned)

        # For each candidate after MQ-tree, do linearization step and check step
        for cand_assigned in solutions_after_mqtree:
            # unknown variables left for linear system: indices 0 .. m-k-p-1
            free_vars = list(range(0, m - k - p))
            num_eqs = (m - p) - (k)  # equations k .. m-p-1 inclusive
            if num_eqs < 0:
                num_eqs = 0
            # build linear system A x = b
            if len(free_vars) == 0:
                # nothing to solve; just build candidate full vector
                u = [None] * m
                for i in range(m):
                    if i in cand_assigned:
                        u[i] = cand_assigned[i]
                    else:
                        u[i] = FF(0)
                # Verify check step: first k equations must hold
                good = True
                # extend to full n vector with zeros for m..n-1
                full_u = [u_i for u_i in u] + [FF(0) for _ in range(n - m)]
                # verify Ptilde[0..k-1]
                for i_eq in range(0, k):
                    # evaluate p(x) = x^T Ptilde[i_eq] x
                    val = FF(0)
                    for a in range(m):
                        va = full_u[a]
                        for b in range(a, m):
                            coef = Ptilde[i_eq][a, b]
                            if coef != 0:
                                val += coef * va * full_u[b]
                    if val != t_tilde[i_eq]:
                        good = False
                        break
                if good:
                    # Map back: u_full_n (length n vector)
                    u_vec = vector(FF, full_u)
                    y = S * u_vec
                    return y
                else:
                    continue

            # build matrices
            FF_mat = Ptilde[0].base_ring()
            Arows = []
            brows = []
            consistent = True
            for eq_i in range(k, m - p):
                const = FF(0)
                coeffs = [FF(0)] * len(free_vars)
                # iterate over all terms a<=b
                for a in range(m):
                    for b in range(a, m):
                        coef = Ptilde[eq_i][a, b]
                        if coef == 0:
                            continue
                        a_known = a in cand_assigned
                        b_known = b in cand_assigned
                        if (a in free_vars) and (b in free_vars):
                            # degree 2 in free vars -> fails (structure should avoid this)
                            # If coef != 0 then not linear -> abandon candidate
                            consistent = False
                            break
                        elif (a in free_vars) and b_known:
                            idx = free_vars.index(a)
                            coeffs[idx] += coef * cand_assigned[b]
                        elif (b in free_vars) and a_known:
                            idx = free_vars.index(b)
                            coeffs[idx] += coef * cand_assigned[a]
                        elif a_known and b_known:
                            const += coef * cand_assigned[a] * cand_assigned[b]
                        else:
                            # term involves unknown not in free_vars and not assigned -> abort
                            consistent = False
                            break
                    if not consistent:
                        break
                if not consistent:
                    break
                # equation: sum coeffs * free_vars = t_tilde[eq_i] - const
                Arows.append(coeffs)
                brows.append(t_tilde[eq_i] - const)
            if not consistent:
                continue

            if len(Arows) == 0:
                # no linear equations to solve; go directly to check step
                # assemble u and test first k equations
                u = [None] * m
                for i in range(m):
                    if i in cand_assigned:
                        u[i] = cand_assigned[i]
                    else:
                        u[i] = FF(0)
                full_u = [u_i for u_i in u] + [FF(0) for _ in range(n - m)]
                good = True
                for i_eq in range(0, k):
                    val = FF(0)
                    for a in range(m):
                        va = full_u[a]
                        for b in range(a, m):
                            coef = Ptilde[i_eq][a, b]
                            if coef != 0:
                                val += coef * va * full_u[b]
                    if val != t_tilde[i_eq]:
                        good = False
                        break
                if good:
                    u_vec = vector(FF, full_u)
                    y = S * u_vec
                    return y
                else:
                    continue

            A = Matrix(FF_mat, Arows)
            b = vector(FF_mat, brows)

            # Solve linear system and enumerate all solutions
            sols = _solve_linear_system_all_solutions(A, b)
            for solv in sols:
                # build full m-vector u
                u = [None] * m
                # free vars (sanity check)
                for i, idx in enumerate(free_vars):
                    if not (0 <= idx < m):
                        raise ValueError(f"free_vars contains invalid index {idx} (m={m})")
                    u[idx] = solv[i]
                # assigned vars from cand_assigned (only indices < m)
                for idx, val in cand_assigned.items():
                    if 0 <= idx < m:
                        u[idx] = val
                # fill remaining with zero
                for idx in range(m):
                    if u[idx] is None:
                        u[idx] = FF(0)
                # Check step: first k equations must match
                full_u = [u_i for u_i in u] + [FF(0) for _ in range(n - m)]
                ok = True
                for i_eq in range(0, k):
                    val = FF(0)
                    for a in range(m):
                        va = full_u[a]
                        for b in range(a, m):
                            coef = Ptilde[i_eq][a, b]
                            if coef != 0:
                                val += coef * va * full_u[b]
                    if val != t_tilde[i_eq]:
                        ok = False
                        break
                if ok:
                    u_vec = vector(FF, full_u)
                    y = S * u_vec
                    return y
        # none of the MQ-tree candidates led to full solution
        return None

    # iterate all assignments for guessed variables (k variables)
    # iterate lexicographically over product(FF, repeat=k)
    if k == 0:
        # single empty tuple
        guess_iter = [tuple()]
    else:
        guess_iter = product(FF, repeat=k)

    for guess_vals in guess_iter:
        y = process_guess_assignment(guess_vals)
        if y is not None:
            return y
    return None


def main(argv=None):
    parser = argparse.ArgumentParser(description='Run just_guess_solver from terminal (sage).')
    parser.add_argument('-n', type=int, required=True, help='ambient dimension n')
    parser.add_argument('-m', type=int, required=True, help='number of polynomials m')
    parser.add_argument('-q', type=int, required=True, help='field size q')
    parser.add_argument('-k', type=int, required=True, help='number of guessed coordinates k')
    parser.add_argument('-p', type=int, required=True, help='number of MQ(1,1) equations p')
    parser.add_argument('--seed', type=int, default=None, help='optional RNG seed')
    args = parser.parse_args(argv)

    n = args.n
    m = args.m
    k = args.k
    p = args.p
    q = args.q

    if args.seed is not None:
        random.seed(args.seed)
        set_random_seed(args.seed)

    print(f'Running just_guess_solver with n={n}, m={m}, q={q}, k={k}, p={p}')
    print('-------------------------------------------')

    # Finite field
    FF = GF(q, 'w')

    average_rerand = 0
    iterations = 1
    for _ in range(iterations):
        solved = False
        attempts = 0
        # original random eqs
        eqs = [uptriag(random_matrix(FF, n)) for i in range(m)]
        target = vector([FF.random_element() for i in range(m)])
        current_eqs = eqs
        M = None
        while not solved:
            print(f'Attempt n{attempts + 1}:')
            print(f'Computing S and T...')
            rerandomisations = 0
            while True:
                try:
                    transform, rank_loss = linear_transformations(FF, current_eqs, k, p)    
                    break
                except:
                    print(f'Randomising the Transformation')
                    rerandomisations += 1
                    M = random_matrix(FF, n)
                    while M.rank() != n:
                        M = random_matrix(FF, n)
                    current_eqs = [M.transpose() * eq * M for eq in eqs]    
            print(f'Rerandomised {rerandomisations} times to compute S, T.')
            print(f'Rank Loss: {rank_loss}')
            y = just_guess_solver(FF, current_eqs, target, k, p, transform)
            if y is None:
                print('Not found... Randomizing!')
                print('-------------------------------------------')
                M = random_matrix(FF, n)
                while M.rank() != n:
                    M = random_matrix(FF, n)
                current_eqs = [M.transpose() * eq * M for eq in eqs]
                attempts += 1
                continue
            print('\x1b[32mFound!\x1b[0m')
            attempts += 1
            if attempts > 1 and M is not None:
                y = M * y
            solved = True
        print('###########################################')
        sol = vector(y)
        print(f'y = {y}')
        print('###########################################')
        print('Verifying the solution:')
        print(f'P(y) = {vector([sol * eq * sol for eq in eqs])}')
        print(f't    = {target}')
        print(f'Number of attempts: {attempts}')
        average_rerand += rerandomisations/attempts
    # print(average_rerand/iterations.n())


if __name__ == '__main__':
    main()
