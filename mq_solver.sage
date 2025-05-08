from itertools import product
import subprocess

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

def to_FF(n, F):
    '''
    Input: a integer n encoding a finite field element in F
    Output: the encoded finite field element
    '''
    p = F.characteristic()
    e = F.degree()
    V, v2e, _ = F.vector_space()
    # Convert integer to base-p representation (least significant digit first)
    coeffs = []
    for _ in range(e):
        coeffs.append(n % p)
        n //= p
    v = V(coeffs)  # Create a vector in the vector space
    return v2e(v)  # Convert back to field element

def to_int(x):
    '''
    Input: a finite field element
    Output: an integer encoding the finite field element
    '''
    F = x.parent()
    p = F.characteristic()
    e2v = F.vector_space()[2]
    coeffs = e2v(x)  # This gives the [b_0, b_1, ..., b_{e-1}]
    return sum(int(c) * (p^i) for i, c in enumerate(coeffs))
    
def field_poly_to_integer_poly(f):
    '''
    Input: a multivariate polynomial over a finite field
    Output: the same polynomial over the integers
    '''
    R = f.parent()             # Polynomial ring
    F = R.base_ring()          # Finite field GF(p^e)
    p = F.characteristic()
    e = F.degree()

    # Define new polynomial ring over ZZ with same variables
    R_ZZ = PolynomialRing(ZZ, R.variable_names())

    # Map coefficients
    new_terms = []
    for mon, coeff in f.dict().items():
        int_coeff = to_int(coeff)
        new_terms.append((mon, int_coeff))
    
    return R_ZZ(dict(new_terms))

def polyring(FF, n_vars):
    '''
    Input: a finite field FF and an integer n
    Output: a multivariate polynomial ring over FF with n variables
    '''
    var_names = ['X' + str(i) for i in range(n_vars)]
    return PolynomialRing(FF, var_names)

def matrix_to_input(R, matrix_list, guess, t):
    '''
    Input: a multivariate polyring R over a finite field F, a list of matrices, a vector with coordinates in FF, a target vector over FF
    Output: the input for M4GB, consisting of the fieldsize, the number of variables, a list of equations in these variables (but the last variables are substituted with the guess) plus the target
    '''
    n_vars = matrix_list[0].nrows()
    m_eqs = len(matrix_list)
    generators = R.gens()
    print(generators)
    variables = vector(generators + guess)
    # assert len(R.gens()) == m_eqs           # n equations == n variables
    eqs = [variables * M * variables for M in matrix_list]
    # for i in range(len(eqs)):
    #     print(eqs[i])
    #     print(matrix_list[i])
    #     print('-------------------------------------')
    MQ_input = '$fieldsize ' + str(R.base_ring().cardinality()) + '\n$vars ' + str(m_eqs) + ' X\n'
    for i in range(len(eqs)):
        MQ_input += str(field_poly_to_integer_poly(eqs[i])).strip() + ' + ' + str(to_int(t[i])) + '\n'
        # MQ_input += str(field_poly_to_integer_poly(eqs[i])).strip() + '\n'
    return MQ_input

def find_guess_hybrid(q, n, omega):
    '''
    Input: fieldsize q, number of variables n and the constant for solving linear systems omega
    Output: the optimal number of guesses for the Hybrid Approach 
    '''
    R.<x> = QQ[]
    return int(round(find_root(log(q) + omega * (log(n - x - 1) + 1 / (2 * (n - x - 1))) \
        - (omega / 2) * (1 + sqrt(n / x)) * (
            log((3 * n - x) / 2 - 1 - sqrt(n * x)) +
            1 / (2 * ((3 * n - x) / 2 - 1 - sqrt(n * x)))
        ) \
        - (omega / 2) * (1 - sqrt(n / x)) * (
            log((n + x) / 2 - sqrt(n * x)) +
            1 / (2 * ((n + x) / 2 - sqrt(n * x)))
        ), 0.5, n-2)))


def process_output(FF, R, file: str):
    '''
    Input: a finite field FF, a polynomial ring R over FF, output from M4GB file
    Output: a solution of the system as a tuple
    '''
    with open(file, "r", encoding="utf-8") as f:
        eqs_str = []
        for line in f.readlines():
            # print(f"Parsing line: {line.strip()}")  # Debugging: print each line
            if 'X' in line and '[m4gb]' not in line:
                eqs_str.append(line.replace("\n", ""))
    eqs = []
    for poly in eqs_str:
        eq = R(0)
        s = poly.split(" ")
        # print(f'poly = {poly}')
        # print(f's={s}')
        # print(f"Line: '{poly}', Split: {s}, Length: {len(s)}")  # Debugging
        if len(s) == 3 and poly.count('X') == 1:
            eq += R(s[0]) + R(to_FF(int(s[2]), FF))     #only in char2
            # Debugging: print the solution being added
            # print(f"Equation: {solution}")
        else:
            for p in s:
                try:
                    eq += R(to_FF(int(p)))
                except:
                    if p == '+' or p == '-':
                        pass
                    try:
                        p[0]
                    except:
                        print(poly)
                        print(s)
                        print(p)
                        raise
                    if p.strip() == "":
                        continue
                    if 'X' == p[0]:
                        eq += R(p)
                    elif 'X' in p and 'MAX' not in p:
                        parts = p.split('*')
                        parts[0] = to_FF(int(parts[0]), FF)
                        eq += prod(R(part) for part in parts)
        eqs.append(eq)
    # print(eqs)
    I = ideal(R, eqs)
    # Compute the variety; we use Sage's variety() which returns a list of dictionaries.
    V = I.variety()
    gens = R.gens()
    return tuple(V[0][i] for i in gens)


def MQ_square_solver(FF, matrix_list, t, guess, hybrid_guess):
    '''
    Input: a finite field FF, a list of matrices, the target over FF, the current guess, the number of coordinates to be guessed in the Hybrid Approach
    Output: the solution of the system as a tuple
    '''
    n_vars = matrix_list[0].nrows()
    m_eqs = len(matrix_list)
    field_size = FF.cardinality()
    R = polyring(FF, m_eqs)
    if n_vars - len(guess) < 3:
        generators = R.gens()
        variables = vector(generators + guess)

        # print(variables)
        # assert len(R.gens()) == m_eqs           # n equations == n variables
        eqs = [variables * matrix_list[i] * variables - t[i] for i in range(len(matrix_list))]
        # print(eqs)
        I = ideal(R, eqs)
        # print(I)
        V = I.variety()
        if V == []:
            return None
        else:
            return tuple(V[0][i] for i in generators) + guess

    MQ_input = matrix_to_input(R, matrix_list, guess, t)
    input_filename = 'input.in'
    for L in product([to_int(element) for element in list(FF)], repeat = hybrid_guess):
        print(L)
        MQ_input_overdef = MQ_input
        for g in range(hybrid_guess):
            MQ_input_overdef = MQ_input_overdef + 'X' + str(g) + ' - ' + str(L[g])
        
        with open(input_filename, 'w') as file:
            file.write(MQ_input_overdef)
            # print(f"Trying hybrid guess: {L}, writing input to {input_filename}")

        command = ['./solver.sh', '-f', str(field_size), '-n', str(n_vars), '-s', 'm4gb', input_filename]
        output_filename = 'output.txt'
        try:
            with open(output_filename, 'w') as output_file:
                result = subprocess.run(command, stdout=output_file, stderr=subprocess.PIPE, text=True)
                # result = subprocess.run(command, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print("Solver returned non-zero exit status (might be no solution for this guess).")
            print("Details:\n", e.stderr)
            continue  # Try next guess
        with open(output_filename, 'r') as f:
            out = f.read()
        # print(out)
        if "contradiction" in out:
            print("No solution for this guess.")
            continue
        if 'X' not in out:
            continue
        print("Solution found!")
        return process_output(FF, R, output_filename) + guess
    return None




# char = 2
# exp = 4
# FF.<w> = GF(2^4)
# n = 10
# m = 7

# eqs = [uptriag(random_matrix(FF, n)) for i in range(m)]
# target = [FF.random_element() for i in range(m)]
# # target = [FF(w^3 + w), FF(w^3 + w), FF(w^2 + 1), FF(0), FF(w^3 + w^2 + w + 1), FF(w^2 + 1), FF(w^3), FF(w^3 + w^2), FF(w^3 + w + 1), FF(w^3 + w)]

# tmp = MQ_square_solver(FF, eqs, target, (w^3 + w,w^2 + 1,w^3 + w + 1), 1)
# if tmp is None:
#     print("No sol")
# else:
#     sol = vector(tmp)
#     print(sol)
#     print([sol * eq * sol for eq in eqs])
#     print(target)