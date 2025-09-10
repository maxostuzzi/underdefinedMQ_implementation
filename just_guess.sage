load('transformation.sage')

import argparse

def compute_x(FF, P, t, p, k, max_initial_attempts=None, verbose=False):
    F = FF

    m = len(P)
    if m == 0:
        raise ValueError("P must be non-empty")
    n = P[0].nrows()

    if not (0 <= k <= m):
        raise ValueError("k must satisfy 0 <= k <= m")
    if not (0 <= p <= m - k):
        raise ValueError("p must satisfy 0 <= p <= m - k")

    transform = linear_transformations(F, P, [m - k - p] + [1 for _ in range(p)], [2 * k + p - m, m - k - p])
    T = transform[0]
    S = transform[1]
    P_tilde_raw = transform[2]
    P_tilde = [ uptriag(M) for M in P_tilde_raw ]
    t_vec = vector(F, t)
    t_tilde = T * t_vec

    Pnames = ['x{}'.format(i) for i in range(n)]
    R = PolynomialRing(F, Pnames)
    gens = list(R.gens())

    variables_vector = vector(list(gens)[0:m] + [R(0) for _ in range(n - m)])
    eqs = [ (variables_vector * P_tilde[i] * variables_vector) - t_tilde[i] for i in range(m) ]

    def poly_roots(poly):
        if poly == 0:
            return list(F)
        raw = poly.roots()
        roots = []
        for item in raw:
            if isinstance(item, tuple):
                roots.append(item[0])
            else:
                roots.append(item)
        return roots

    def univariate_poly_for_var(f_poly, var_index, assigned_map):
        PR = PolynomialRing(F, 't')
        t = PR.gen()
        gens_in_f = set(f_poly.variables())
        images = []
        for j, g in enumerate(gens):
            if j == var_index:
                images.append(t)
            else:
                if g in gens_in_f:
                    if j not in assigned_map:
                        raise ValueError("Variable x{} (appears in f) not assigned".format(j))
                    images.append(PR(F(assigned_map[j])))
                else:
                    images.append(PR(0))
        phi = R.hom(images, PR)
        return phi(f_poly)

    def eval_residual(f_poly, full_assignment):
        images_eval = [ full_assignment[j] if j in full_assignment else F(0) for j,_ in enumerate(gens) ]
        phi = R.hom(images_eval, F)
        return phi(f_poly)

    last_k_indices = list(range(m - k, m)) if k > 0 else []

    if max_initial_attempts is None:
        max_initial_attempts = float('inf')
    field_elems = list(F)
    total_possible = (len(field_elems)**k) if k > 0 else 1
    initial_iter = product(field_elems, repeat=k) if k > 0 else [()]

    tried = 0

    # DFS that returns full_assign dict or None
    def dfs_solve_eq_index(i, assigned_map):
        if i < m - p:
            # remaining vars among first m
            remaining_vars = [j for j in range(0, m) if j not in assigned_map]
            num_remaining = len(remaining_vars)

            if num_remaining == 0:
                full_assign = { idx: assigned_map[idx] for idx in range(0, m) }
                for idx in range(m, n):
                    full_assign[idx] = F(0)
                for j in range(m):
                    if eval_residual(eqs[j], full_assign) != 0:
                        return None
                return full_assign

            # we'll collect linear rows (j, coeffs_dict[F], const_F) and nonlinear indices
            linear_rows = []
            nonlinear_rows = []

            # Precompute two kinds of homomorphisms per row:
            # - phi_const: assigned_map -> values, remaining vars -> 0
            # - phi_u1: assigned_map -> values, var u -> 1, other remaining -> 0
            for j in range(0, m - p):
                # images_const maps assigned -> R(constants), remaining -> 0, last n-m -> 0
                images_const = []
                for ind, g in enumerate(gens):
                    if ind < m:
                        if ind in assigned_map:
                            images_const.append(R(F(assigned_map[ind])))
                        else:
                            images_const.append(R(0))
                    else:
                        images_const.append(R(0))
                # homomorphism to R but since images_const are constants, we can map directly to F
                phi_const_toF = R.hom([F(assigned_map[ind]) if (ind < m and ind in assigned_map) else F(0) for ind,_ in enumerate(gens)], F)
                const_F = phi_const_toF(eqs[j])

                # now get linear coefficients by setting each remaining var to 1 in turn
                linear_coeffs = {var_idx: F(0) for var_idx in remaining_vars}
                is_linear = True
                # To compute coeff_u: evaluate with u=1 and other remaining = 0, then subtract const_F
                for u in remaining_vars:
                    # build images_u1 mapping to field F (constants): assigned -> value, u -> 1, other remaining -> 0, last n-m -> 0
                    images_toF = []
                    for ind, _ in enumerate(gens):
                        if ind < m:
                            if ind in assigned_map:
                                images_toF.append(F(assigned_map[ind]))
                            else:
                                if ind == u:
                                    images_toF.append(F(1))
                                else:
                                    images_toF.append(F(0))
                        else:
                            images_toF.append(F(0))
                    phi_u1 = R.hom(images_toF, F)
                    eval_u1 = phi_u1(eqs[j])
                    coeff_u = eval_u1 - const_F
                    linear_coeffs[u] = F(coeff_u)

                                # Now reconstruct linear polynomial in R: linear_R = R(const_F) + sum coeff_u * gens[u]
                linear_R = R(F(const_F))
                for u in remaining_vars:
                    linear_R = linear_R + R(F(linear_coeffs[u])) * gens[u]

                # Build images_R for full substitution except remaining vars remain symbolic
                images_R = []
                for ind, g in enumerate(gens):
                    if ind < m:
                        if ind in assigned_map:
                            images_R.append(R(F(assigned_map[ind])))
                        else:
                            images_R.append(g)
                    else:
                        images_R.append(R(0))
                phi_R = R.hom(images_R, R)
                poly_sub = phi_R(eqs[j])

                diff = poly_sub - linear_R

                if diff.is_zero():
                    linear_rows.append((j, linear_coeffs, const_F))
                else:
                    nonlinear_rows.append(j)


            # If there are linear rows, build linear system from them (using const_F and coeffs in F)
            if len(linear_rows) > 0:
                rows = len(linear_rows)
                cols = num_remaining
                A = Matrix(F, rows, cols, lambda rr, cc: F(0))
                b = vector(F, [F(0) for _ in range(rows)])
                for row_idx, (j, lin_coeffs, const_F) in enumerate(linear_rows):
                    for col_idx, var_idx in enumerate(remaining_vars):
                        A[row_idx, col_idx] = lin_coeffs[var_idx]
                    b[row_idx] = -F(const_F)

                aug = A.augment(b)
                if A.rank() != aug.rank():
                    return None

                try:
                    y0 = A.solve_right(b)
                except Exception:
                    return None
                y0 = vector([F(c) for c in list(y0)])

                K = A.right_kernel()
                basis = K.basis()
                d = len(basis)

                # enumerate all affine solutions
                if d == 0:
                    candidates = [y0]
                else:
                    candidates = []
                    for coeffs in product(field_elems, repeat=d):
                        y = vector(list(y0))
                        for idx_b, c in enumerate(coeffs):
                            if c == 0:
                                continue
                            y = y + F(c) * vector(basis[idx_b])
                        candidates.append(y)

                for y in candidates:
                    full_assign = dict(assigned_map)
                    for col_idx, var_idx in enumerate(remaining_vars):
                        full_assign[var_idx] = F(y[col_idx])
                    for idx in range(m, n):
                        full_assign[idx] = F(0)
                    ok = True
                    for j in range(m):
                        if eval_residual(eqs[j], full_assign) != 0:
                            ok = False
                            break
                    if ok:
                        return full_assign
                return None
            else:
                # No linear rows: enumerate all assignments on remaining_vars
                if verbose:
                    print("No linear eqs among first m-p; enumerating all {}^{} assignments.".format(F.order(), num_remaining))
                for tup2 in product(field_elems, repeat=num_remaining):
                    full_assign = dict(assigned_map)
                    for idx_col, var_idx in enumerate(remaining_vars):
                        full_assign[var_idx] = tup2[idx_col]
                    for idx in range(m, n):
                        full_assign[idx] = F(0)
                    ok = True
                    for j in range(m):
                        if eval_residual(eqs[j], full_assign) != 0:
                            ok = False
                            break
                    if ok:
                        return full_assign
                return None

        # DFS case i >= m-p
        start = max(0, i - k)
        vars_indices = list(range(start, m))
        free_vars = [j for j in vars_indices if j not in assigned_map]

        if len(free_vars) == 0:
            residual = eval_residual(eqs[i], assigned_map)
            if residual != 0:
                return None
            else:
                return dfs_solve_eq_index(i-1, assigned_map)

        if len(free_vars) > 1:
            return None

        v = free_vars[0]
        try:
            poly_univ = univariate_poly_for_var(eqs[i], v, assigned_map)
        except ValueError:
            return None
        roots = poly_roots(poly_univ)
        if not roots:
            return None
        for r in roots:
            assigned_map[v] = F(r)
            res = dfs_solve_eq_index(i-1, assigned_map)
            if res is not None:
                return res
            del assigned_map[v]
        return None

    # main loop
    for tup in initial_iter:
        if k > 0:
            tried += 1
            if tried > max_initial_attempts:
                break
            assigned_map = { idx: tup[j] for j, idx in enumerate(last_k_indices) }
        else:
            assigned_map = {}
        if verbose:
            if k > 0:
                printable = { 'x{}'.format(idx): assigned_map[idx] for idx in assigned_map }
                print("\n---\nTrying initial assignment #{} / {} : {}".format(min(tried, total_possible), total_possible, printable))
        full_assign = dfs_solve_eq_index(m-1, assigned_map)
        if full_assign is not None:
            return vector([ full_assign[i] for i in range(n) ]), S*vector([ full_assign[i] for i in range(n) ])
    return None


def main():
    parser = argparse.ArgumentParser(description="Find k and p that satisfy inequalities and minimize q^k")
    parser.add_argument("n", type=int, help="integer n")
    parser.add_argument("m", type=int, help="integer m (>=0)")
    parser.add_argument("k", type=int, help="integer k")
    parser.add_argument("p", type=int, help="integer p")
    args = parser.parse_args()
    print(f'Solving MQ problem with Just-Guess-Solver with:')
    print(f'n = {args.n}')
    print(f'm = {args.m}')
    print(f'q = {7}')
    print(f'k = {args.k}')
    print(f'p = {args.p}')

    FF = GF(7)
    eqs = [uptriag(random_matrix(FF, args.n)) for i in range(args.m)]
    target = [FF.random_element() for i in range(args.m)]
    try:
        x, Sx = compute_x(FF, eqs, target, args.p, args.k, verbose = False)
        print(f'The solution is')
        print(f'y = {Sx}')
        print(f'Check whether Py = t...')
        sol = vector(Sx)
        if [sol * eq * sol for eq in eqs]:
            print(f'True!')
        else:
            print('Wrong...')
    except:
        print("No solution")


if __name__ == "__main__":
    main()