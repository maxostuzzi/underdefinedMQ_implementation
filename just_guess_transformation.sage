from itertools import product
import subprocess
import argparse

def uptriag(matrix):
    '''
    Input: a matrix matrix
    Output: Upper(matrix)
    '''
    for i in range(matrix.nrows()):
        for j in range(i):
            matrix[j, i] = matrix[i, j] + matrix[j, i]
            matrix[i, j] = 0
    return matrix

def elim_square(FF, matrix_list, eliminator_index, index):
    '''
    Generalised Gaussian elimination for MQ systems of equations. It makes 0 the coefficient of the square monomial of index "index", in all matrices starting from "eliminator_index + 1"
    '''
    # eliminator_index -= 1
    # index -= 1
    T = identity_matrix(FF, len(matrix_list))
    counter = 0
    pivot_found = True
    if matrix_list[eliminator_index][index, index] == 0:
        pivot_found = False
        if eliminator_index == len(matrix_list) - 1:
            return T
        while counter <= len(matrix_list) - eliminator_index - 1:
            if matrix_list[eliminator_index + counter][index, index] != 0:
                pivot_found = True
                break
            counter += 1
        if not pivot_found:
            return T
    for i in list(range(eliminator_index, eliminator_index + counter)) + list(range(eliminator_index + counter + 1, len(matrix_list))):
        factor = matrix_list[eliminator_index + counter][index, index].inverse() * matrix_list[i][index, index]
        T[i, eliminator_index + counter] = - factor 
    T[eliminator_index], T[eliminator_index + counter] = T[eliminator_index + counter], T[eliminator_index] 
    return T


def fetch_constraints_first(FF, matrix_list, k, p, col_index):
    '''
    Fetches the linear constraints in the case col_index is in A_1(-1). Col_index here starts from 0.
    '''
    m_eqs = len(matrix_list)
    assert 0 < col_index and col_index < m_eqs - k - p
    n_vars = matrix_list[0].nrows()
    constraints = Matrix(FF, zero_vector(FF, n_vars))
    for M in matrix_list[2 * k + p - m_eqs + col_index:]:
        constraints = constraints.stack(M[:col_index])
    return constraints[1:]

# def fetch_constraints_second(FF, matrix_list, k, p, col_index):
#     '''
#     Fetches the constraints in the case that col_index is in A_2(-1).
#     '''
#     m_eqs = len(matrix_list)
#     assert 0 <= col_index and col_index < p
#     n_vars = matrix_list[0].nrows()
#     constraints = Matrix(FF, zero_vector(FF, n_vars))

#     for M in matrix_list[m_eqs - p + col_index:]:
#         constraints = constraints.stack(M[:m_eqs - k - p + col_index])

#     if col_index > 0:
#         for i in range(col_index):
#             constraints = constraints.stack(matrix_list[m_eqs - p +i][:m_eqs - p + i])

#     return constraints[1:]

def fetch_constraints_second(FF, matrix_list, k, p, col_index):
    '''
    Fetches the constraints in the case that col_index is in A_2(-1).
    '''
    m_eqs = len(matrix_list)
    assert m_eqs - k - p <= col_index and col_index < m_eqs - k
    n_vars = matrix_list[0].nrows()
    constraints = Matrix(FF, zero_vector(FF, n_vars))

    for M in matrix_list[k + col_index:]:
        constraints = constraints.stack(M[: col_index])

    if col_index > m_eqs - k - p:
        for i in range(m_eqs - k - p, col_index):
            constraints = constraints.stack(matrix_list[k + i][: i])

    return constraints[1:]

def fetch_constraints_third(FF, matrix_list, k, p, col_index):
    '''
    Fetches the constraints in the case col_index is in A_3(-1).
    '''
    n_vars = matrix_list[0].nrows()
    m_eqs = len(matrix_list)
    assert m_eqs - k <= col_index and col_index < m_eqs
    constraints = Matrix(FF, zero_vector(FF, n_vars))

    for i in range(m_eqs - p, m_eqs):
        constraints = constraints.stack(matrix_list[i][:i - k])

    return constraints[1:]

def elim_above(FF, matrix_list, k, p, col_index):
    '''
    Eliminates the coefficients of mixed terms above the diagonal of the column col_index.
    It returns the (transposed) desired change of variables in the domain of the MQ map.
    '''
    m_eqs = len(matrix_list)
    if col_index < m_eqs - k - p:
        constraints = fetch_constraints_first(FF, matrix_list, k, p, col_index)
    elif m_eqs - k - p <= col_index and col_index < m_eqs - k:
        constraints = fetch_constraints_second(FF, matrix_list, k, p, col_index)
    elif m_eqs - k <= col_index and col_index < m_eqs:
        constraints = fetch_constraints_third(FF, matrix_list, k, p, col_index)
    # Number of variables
    n_vars = matrix_list[0].nrows()
    # Pad with rows of zeros
    while constraints.nrows() < n_vars - 1:
        constraints = constraints.stack(zero_vector(FF, n_vars))
    # Fix the diagonal coefficient to 1, for invertibility
    inv_constraint = zero_vector(FF, n_vars)
    # Double check, there might be a +something
    inv_constraint[col_index] = 1
    constraints = constraints.stack(inv_constraint)
    # assert constraints.nrows() > n_vars, 'Error: Not enough variables'
    # Solve the system
    try:
        S_col = constraints.solve_right(vector(FF, [0]*(n_vars-1) + [1]))
    except ValueError as e:
        print("Failed to solve constraints matrix:")
        print(constraints.str())
        raise
    # Build transformation matrix
    # Double check, there might be a +something
    S = identity_matrix(FF, n_vars)
    S[col_index] = S_col
    return S

def elim_above_guessing(FF, matrix_list, k, p, col_index):
    '''
    Eliminates the coefficients of mixed terms above the diagonal of the column col_index.
    It returns the (transposed) desired change of variables in the domain of the MQ map.
    '''
    m_eqs = len(matrix_list)
    assert m_eqs - k <= col_index and col_index < m_eqs
    constraints = fetch_constraints_third(FF, matrix_list, k, p, col_index - (m_eqs - k))
    # Number of variables
    n_vars = matrix_list[0].nrows()
    # Pad with rows of zeros
    while constraints.nrows() < n_vars - 1:
        constraints = constraints.stack(zero_vector(FF, n_vars))
    #switch to 0-based
    col_index -= 1
    # Fix the diagonal coefficient to 1, for invertibility
    inv_constraint = zero_vector(FF, n_vars)
    # Double check, there might be a +something
    inv_constraint[col_index] = 1
    constraints = constraints.stack(inv_constraint)
    # assert constraints.nrows() > n_vars, 'Error: Not enough variables'
    # Solve the system
    try:
        S_col = constraints.solve_right(vector(FF, [0]*(n_vars-1) + [1]))
    except ValueError as e:
        print("Failed to solve constraints matrix:")
        print(constraints.str())
        raise
    # Build transformation matrix
    # Double check, there might be a +something
    S = identity_matrix(FF, n_vars)
    S[col_index] = S_col
    return S

def one_step(FF, matrix_list, k, p, col_index):
    '''
    It returns a list of 2 elements: the 1st is the change of variables in the codomain after handling the column col_index, and the 2nd is the change of variables in the domain after handling the column col_index
    '''
    m_eqs = len(matrix_list)
    if col_index == 0:
        return [elim_square(FF, matrix_list, 2 * k + p - m_eqs, 0), identity_matrix(FF, matrix_list[0].nrows())]
    elif 0 < col_index and col_index < m_eqs - k - p:
        S = elim_above(FF, matrix_list, k, p, col_index)
        copy = [uptriag(S * M * S.transpose()) for M in matrix_list]
        return [elim_square(FF, copy, 2 * k + p - m_eqs + col_index, col_index), S]
    elif m_eqs - k - p <= col_index and col_index < m_eqs - k:
        S = elim_above(FF, matrix_list, k, p, col_index)
        copy = [uptriag(S * M * S.transpose()) for M in matrix_list]
        return [elim_square(FF, copy, k + col_index, col_index), S]
    else:
        S = elim_above(FF, matrix_list, k, p, col_index)
        return [identity_matrix(FF, len(matrix_list)), S]

def linear_transformations(FF, matrix_list, k, p):
    '''
    It returns a list of 3 elements: 1st is the full change of variables in the codomain, 2nd is the full change of variables in the domain and 3rd is the MQ system after the change of variables
    '''
    n = matrix_list[0].nrows()
    m = len(matrix_list)
    S, T = identity_matrix(FF, n), identity_matrix(FF, m)
    after_S = [M for M in matrix_list]
    after_TS = [M for M in matrix_list]
    for column_index in range(m):
        try:
            T1, S1 = one_step(FF, after_TS, k, p, column_index)
        except ValueError as e:
            print(f'\x1b[31mFailed at {column_index+1}-th column\x1b[0m')
            raise
        T = T1 * T
        S = S1 * S
        after_S = [uptriag(S * M * S.transpose()) for M in matrix_list]
        after_TS = [sum(T[h, k] * after_S[k] for k in range(m)) for h in range(m)]
    return [T, S.transpose(), after_TS], n - S.rank()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tranformation Algorithm")
    parser.add_argument("-n", type=int, help="Number of variables")
    parser.add_argument("-m", type=int, help="Number of equations")
    parser.add_argument("-q", type=int, help="Field Size")
    parser.add_argument("-k", type=int, help="Field Size")
    parser.add_argument("-p", type=int, help="Field Size")
    arg = parser.parse_args()
    FF.<w> = GF(arg.q)
    eqs = [uptriag(random_matrix(FF, arg.n)) for i in range(arg.m)]
    current_eqs = [Mat for Mat in eqs]
    print(f'Parameters:')
    print(f'q = {arg.q}')
    print(f'n = {arg.n}')
    print(f'm = {arg.m}')
    print(f'k = {arg.k}')
    print(f'p = {arg.p}')
    solved = False
    rerand_counter = 0
    while solved == False:
        try:
            transform, rank_loss = linear_transformations(FF, current_eqs, arg.k, arg.p)
            solved = True
        except:
            rerand_counter += 1
            M = random_matrix(FF, arg.n)
            while not M.is_invertible():
                M = random_matrix(FF, arg.n)
            current_eqs = [uptriag(M.transpose() * eq * M) for eq in eqs]
    eqs_transformed = transform[2]
    print('Success! S and T have successfully been found. ')
    for M in eqs_transformed[:arg.k]:
        print(M)
        print('######################')
    print('######################')
    print('######################')
    for M in eqs_transformed[arg.k: arg.m - arg.p]:
        print(M)
        print('######################')
    print('######################')
    print('######################')
    for M in eqs_transformed[arg.m - arg.p:]:
        print(M)
        print('######################')
    print(f'Necessary rerands: {rerand_counter}')


