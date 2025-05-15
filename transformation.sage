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
    Generalised Gaussian elimination for MQ systems of equations. It makes 0 the coefficient of the square monomial of index index, in all matrices starting from eliminator_index + 1
    '''
    eliminator_index -= 1
    index -= 1
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

def fetch_constraints_first(FF, matrix_list, partition, guessed, col_index):
    '''
    Fetches the linear equations in the case of the first part b1 of the partition
    '''
    col_index -= 1
    assert 0 < col_index and col_index < partition[0]
    n_vars = matrix_list[0].nrows()
    constraints = Matrix(FF, zero_vector(FF, n_vars))
    for M in matrix_list[guessed[0] + col_index:]:
        constraints = constraints.stack(M[:col_index])
    return constraints[1:]

def fetch_constraints(FF, matrix_list, partition, guessed, p_index, col_index):
    '''
    Fetches the constraints in the case of intermediate parts bp_index of the partition (i.e. not b1 and bp)
    '''
    p_index -= 1
    col_index -= 1
    assert 0 < p_index and p_index < len(partition) - 1
    assert 0 <= col_index and col_index < partition[p_index]
    n_vars = matrix_list[0].nrows()
    constraints = Matrix(FF, zero_vector(FF, n_vars))

    for p in range(1, p_index):
        for M in matrix_list[sum(guessed) + sum(partition[:p]) : sum(guessed) + sum(partition[:p + 1])]:
            constraints = constraints.stack(M[:sum(partition[:p])])

    for M in matrix_list[sum(guessed) + sum(partition[:p_index]) : sum(guessed) + sum(partition[:p_index]) + col_index]:
        constraints = constraints.stack(M[:sum(partition[:p_index])])

    for M in matrix_list[sum(guessed) + sum(partition[:p_index]) + col_index :]:
        constraints = constraints.stack(M[: sum(partition[:p_index]) + col_index])
    return constraints[1:]

def fetch_constraints_last(FF, matrix_list, partition, guessed, col_index):
    '''
    Fetches the constraints for the last part bp of the partition
    '''
    col_index -= 1
    n_vars = matrix_list[0].nrows()
    constraints = Matrix(FF, zero_vector(FF, n_vars))
    for p in range(1, len(partition) - 1):
        # Double check here what happens if it overshoots
        for M in matrix_list[sum(guessed) + sum(partition[:p]) : sum(guessed) + sum(partition[:p + 1])]:
            constraints = constraints.stack(M[: sum(partition[:p])])
    for M in matrix_list[sum(guessed) + sum(partition[:len(partition) - 1]) : ]:
        constraints = constraints.stack(M[:sum(partition[:len(partition) - 1])])
    return constraints[1:]

def elim_above(FF, matrix_list, partition, guessed, p_index, col_index):
    '''
    Eliminates the coefficients of mixed terms above the diagonal of the column col_index.
    It returns the (transposed) desired change of variables in the domain of the MQ map.
    '''
    if p_index == 1:
        constraints = fetch_constraints_first(FF, matrix_list, partition, guessed, col_index)
    elif p_index == len(partition):
        constraints = fetch_constraints_last(FF, matrix_list, partition, guessed, col_index)
    else:
        constraints = fetch_constraints(FF, matrix_list, partition, guessed, p_index, col_index)
    # Number of variables
    n_vars = matrix_list[0].nrows()
    # Pad with rows of zeros
    while constraints.nrows() < n_vars - 1:
        constraints = constraints.stack(zero_vector(FF, n_vars))
    #switch to 0-based
    p_index -= 1
    col_index -= 1
    # Fix the diagonal coefficient to 1, for invertibility
    inv_constraint = zero_vector(FF, n_vars)
    # Double check, there might be a +something
    inv_constraint[sum(partition[:p_index]) + col_index] = 1
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
    S[sum(partition[:p_index]) + col_index] = S_col
    return S

def elim_above_guessing(FF, matrix_list, partition, guessed, col_index):
    '''
    Eliminates the coefficients of mixed terms above the diagonal of the column col_index.
    It returns the (transposed) desired change of variables in the domain of the MQ map.
    '''
    constraints = fetch_constraints_last(FF, matrix_list, partition, guessed, col_index)
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
    inv_constraint[sum(partition) + col_index] = 1
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
    S[sum(partition) + col_index] = S_col
    return S

def one_step(FF, matrix_list, partition, guessed, p_index, col_index):
    '''
    It returns a list of 2 elements: the 1st is the change of variables in the codomain after handling the column col_index, and the 2nd is the change of variables in the domain after handling the column col_index
    '''
    if p_index == 1 and col_index == 1:
        return [elim_square(FF, matrix_list, guessed[0] + col_index, col_index), identity_matrix(FF, matrix_list[0].nrows())]
    elif p_index == 1:
        S = elim_above(FF, matrix_list, partition, guessed, p_index, col_index)
        copy = [uptriag(S * M * S.transpose()) for M in matrix_list]
        return [elim_square(FF, copy, guessed[0] + col_index, col_index), S]
    elif p_index == len(partition):
        S = elim_above(FF, matrix_list, partition, guessed, p_index, col_index)
        return [identity_matrix(FF, len(matrix_list)), S]
    else:
        S = elim_above(FF, matrix_list, partition, guessed, p_index, col_index)
        copy = [uptriag(S * M * S.transpose()) for M in matrix_list]
        return [elim_square(FF, copy, sum(guessed) + sum(partition[: p_index - 1]) + col_index, sum(partition[: p_index - 1]) + col_index), S]

def linear_transformations(FF, matrix_list, partition, guessed):
    '''
    It returns a list of 3 elements: 1st is the full change of variables in the codomain, 2nd is the full change of variables in the domain and 3rd is the MQ system after the change of variables
    '''
    assert sum(partition) + sum(guessed) == len(matrix_list)
    n = matrix_list[0].nrows()
    m = len(matrix_list)
    S, T = identity_matrix(FF, n), identity_matrix(FF, m)
    after_S = [M for M in matrix_list]
    after_TS = [M for M in matrix_list]
    for i in range(len(partition)):
        for j in range(partition[i]):
            try:
                T1, S1 = one_step(FF, after_TS, partition, guessed, i + 1, j + 1)
            except ValueError as e:
                print(f'\x1b[31mFailed at {i+1}-th block, {j+1}-th column\x1b[0m')
                raise
            T = T1 * T
            S = S1 * S
            after_S = [uptriag(S * M * S.transpose()) for M in matrix_list]
            after_TS = [sum(T[h, k] * after_S[k] for k in range(m)) for h in range(m)]
    for i in range(sum(guessed)):
        S1 = elim_above_guessing(FF, after_TS, partition, guessed, i + 1)
        S = S1 * S
        after_S = [uptriag(S * M * S.transpose()) for M in matrix_list]
        after_TS = [sum(T[h, k] * after_S[k] for k in range(m)) for h in range(m)]
    return [T, S.transpose(), after_TS]



########################################################################################################################################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tranformation Algorithm")
    parser.add_argument("-n", type=int, help="Number of variables")
    parser.add_argument("-m", type=int, help="Number of equations")
    parser.add_argument("-q", type=int, help="Field Size")
    parser.add_argument("-partition", nargs='+', type=int, help="Partition")
    parser.add_argument("-guessed", nargs='+', type=int, help="Guessed")    
    arg = parser.parse_args()
    assert sum(arg.partition)+sum(arg.guessed) == arg.m, f'Error: the partition does not sum to {arg.m}'
    FF.<w> = GF(arg.q)
    eqs = [uptriag(random_matrix(FF, arg.n)) for i in range(arg.m)]
    print(f'Parameters:')
    print(f'q = {arg.q}')
    print(f'n = {arg.n}')
    print(f'm = {arg.m}')
    eqs_transformed = linear_transformations(FF, eqs, arg.partition, arg.guessed)[2]
    print('Success! S and T have successfully been found. ')
    print('Printing the final form of the matrices...')
    print(f'{arg.guessed[0]} unused Matrices')
    for j in range(arg.guessed[0]):
        print(eqs_transformed[j])
        print('--------------------------------------------------------------------------------------')
    print('--------------------------------------------------------------------------------------')
    print(f'First Block of size {arg.partition[0]}')
    for j in range(arg.guessed[0], arg.guessed[0]+arg.partition[0]):
        print(eqs_transformed[j])
        print('--------------------------------------------------------------------------------------')
    print('--------------------------------------------------------------------------------------')
    print(f'{arg.guessed[1]} to linearize')
    for j in range(arg.guessed[0]+arg.partition[0], arg.guessed[0]+arg.partition[0]+arg.guessed[1]):
        print(eqs_transformed[j])
        print('--------------------------------------------------------------------------------------')
    print('--------------------------------------------------------------------------------------')
    for i in range(1, len(arg.partition)):
        print(f'{i+1}-th Block of size {arg.partition[i]}')
        for j in range(sum(arg.guessed)+sum(arg.partition[:i]), sum(arg.guessed)+sum(arg.partition[:i+1])):
            print(eqs_transformed[j])
            print('--------------------------------------------------------------------------------------')
        print('--------------------------------------------------------------------------------------')