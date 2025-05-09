from itertools import product

load('mq_solver.sage')
load('transformation.sage')
#load('underdefined_solver2.sage')


# linear_transformations(FF, matrix_list, partition, guessed)

# MQ_square_solver(FF, matrix_list, t, guess, hybrid_guess):

# last_variables(FF, matrix_list, partition, guessed):

def last_variables(FF, matrix_list, partition, guessed):
    n_vars = matrix_list[0].nrows()
    m_eqs = len(matrix_list)
    k = sum(guessed)
    constraints = Matrix(FF, zero_vector(FF, n_vars - m_eqs))

    #print("constraints shape:", constraints.nrows(), constraints.ncols())

    for p in range(1, len(partition)):
        for M in matrix_list[k + sum(partition[:p]) : k + sum(partition[: p + 1])]:

            #print(M[:sum(partition[:p]),sum(partition) + k:])

            constraints = constraints.stack(M[:sum(partition[:p]), m_eqs:])
    
    #print("constraints shape:", constraints.nrows(), constraints.ncols())

    constraints = constraints[1:]
    while constraints.nrows() < n_vars - m_eqs:
        constraints = constraints.stack(zero_vector(FF, n_vars - m_eqs))
    print(constraints.ncols())
    print(constraints.nrows())
    last_vars = constraints.right_kernel()
    return last_vars


def compute_x(FF, P, t, partition, k, q, omega = 3):
    """Computes x = (x_1, ... ,x_n)"""
    k_total = sum(k)

    # computes (T,S,\tilde{P}) and assigns them to corresponding variables
    transform = linear_transformations(FF, P, partition, k)
    T = transform[0]
    S = transform[1]
    P_tilde = transform[2]
    t_tilde = T * vector(FF, t)
    print(t_tilde)

    # Fixing last variables
    last = last_variables(FF, P_tilde, partition, k) # columns are not equal
    vectorize_last = last.basis()[0]

    sol = None

    # Try to find a solution with some u otherwise select different u
    for u in product([element for element in list(FF)], repeat = sum(k)): # integer rep or field rep? ._integer_representation()
        print(u)
        vecorize_u = vector(FF, u)

        sol = vecorize_u.concatenate(vectorize_last)
        
        solved = True

        for i in reversed(range(1, len(partition))):
            try:
                print(f'Solving {i+1}-th block')
                print(f'The Hybrid Guess is {find_guess_hybrid(partition[i], q, omega)}')
                try:
                    hybrid_guess = find_guess_hybrid(partition[i], q, omega)
                except:
                    hybrid_guess = partition[i]

                sol = MQ_square_solver(FF, [P[sum(partition[:i]): , sum(partition[:i]):] for P in P_tilde[k_total + sum(partition[:i]): k_total + sum(partition[:i+1])]], t_tilde[k_total + sum(partition[:i]):], tuple(sol), 1)
                if sol is None:
                    solved = False
                    print(f'Failed {i+1}-th block')
                    break
                else:
                    print(f'Solved {i+1}-th block')
                    #sol = vector(FF, sol_).concatenate(sol)
            except:
                solved = False
                raise
                break
        if not solved:
            continue
        sol = MQ_overdefined_solver(FF, P_tilde[:(partition[0] + k_total)], t_tilde[:(partition[0] + k_total)], tuple(sol), 1)
        if sol is None: 
            print(f'Failed at last block')
            continue
        else:
            return S * vector(FF, sol)
    return None
        

        

    """
            try:
                sol_ = MQ_square_solver(FF, [P[partition[0]:, partition[0]:] for P in P_tilde[partition[0]:(partition[0] + partition[1])]], t[partition[0]:(partition[0] + partition[1])], sol, find_guess_hybrid(partition[1], q, omega))
            if sol_ is None:
                continue
            else:
                sol = sol_ + sol
            except:
                continue

        try:
            sol_ = MQ_square_solver(FF, [P[partition[0]:, partition[0]:] for P in P_tilde[partition[0]:(partition[0] + partition[1])]], t[partition[0]:(partition[0] + partition[1])], sol, find_guess_hybrid(partition[1], q, omega))
            if sol_ is None:
                continue
            else:
                sol = sol_ + sol
        except:
            continue
        
        try:
            sol_ = MQ_square_solver(FF, [P for P in P_tilde[:partition[0]]], t[:partition[0]], sol, find_guess_hybrid(partition[1], q, omega))
            if sol_ is None:
                continue
            else: 
                sol = sol_ + sol
        except:
            continue
    
    return sol
    """
    





q = 2^3
FF.<w> = GF(q)
B = [3,3,3]
k1 = 1
k2 = 1
k = [k1, k2]
n = 73
m = sum(B) + sum(k) 
omega = 3

assert sum(B)+sum(k) == m, f'Error: the partition does not sum to {m}'

P = [uptriag(random_matrix(FF, n)) for i in range(m)]
t = [FF.random_element() for i in range(m)]

"""res = compute_x(FF, P, t, B, k, q, omega)
print("------------------------")
print(res)"""


tmp = compute_x(FF, P, t, B, k, q, omega)
print("-----------------------------------------------")
if tmp is None:
    print("No sol")
else:
    sol = vector(tmp)
    print(sol)
    print([sol * eq * sol for eq in P])
    print(t)