

# This file was *autogenerated* from the file prob_exp_test.sage
from sage.all_cmdline import *   # import sage library

_sage_const_0 = Integer(0); _sage_const_1 = Integer(1); _sage_const_4 = Integer(4); _sage_const_2 = Integer(2); _sage_const_3 = Integer(3)
import argparse 

def generate_instance(FF, n):
	matrix_list = [Matrix(FF, [FF(_sage_const_0 ) for i in range(n)]) for j in range(n)]
	for i in range(len(matrix_list)):
		for j in range(n - i - _sage_const_1 ):
			matrix_list[i] = matrix_list[i].stack(vector([FF(_sage_const_0 ) for k in range(n)]))
		for j in range(n - i - _sage_const_1 , n):
			matrix_list[i] = matrix_list[i].stack(vector([FF(_sage_const_0 ) for k in range(n - i - _sage_const_1 )] + [FF.random_element() for k in range(i + _sage_const_1 )]))
	return [M[_sage_const_1 :] for M in matrix_list]

def polyring(FF, n_vars):
	var_names = ['X' + str(i) for i in range(n_vars)]
	return PolynomialRing(FF, var_names)

def matrix_to_eq(R, matrix, given):
	generators = R.gens()
	variables = vector(generators + given)
	return variables * matrix * variables


from typing import List, Dict, Generator, Any, Tuple


def solve_triangular_matrices(
	FF: Any,
	n: int,
	mats,
	target,
	show: bool = False,
	first_only: bool = True
) -> Tuple[Tuple[Any, Ellipsis] | None, int]:
	"""
	Solves a triangular system of quadratic equations incrementally and counts the number of attempted root lookups.
	Returns the first found solution (if any) and the number of attempts made.
	"""

	SR = PolynomialRing(FF, names=('T',)); (T,) = SR._first_ngens(1)
	m = n
	found_solution = None
	total_attempts = _sage_const_0 

	def backtrack(level: int, partial: List[Any]) -> Generator[Tuple[Any, Ellipsis], None, None]:
		nonlocal total_attempts, found_solution

		if level == m:
			found_solution = tuple(partial)
			yield found_solution
			return

		R = polyring(FF, n - level)
		eq = matrix_to_eq(R, mats[level], tuple(partial))
		variables = [_sage_const_0  for _ in range(n - level - _sage_const_1 )] + [T]
		P = eq(variables) - target[n - level - _sage_const_1 ]
		if show:
			print(P)
		if not show:
			if level == _sage_const_0 :
				P += FF.random_element() * T
		total_attempts += _sage_const_1 

		if P != _sage_const_0 :
			roots = P.roots(multiplicities=False)
		else:
			roots = list(FF)

		for r in roots:
			partial = [FF(r)] + partial
			if show:
				print(partial)
			yield from backtrack(level + _sage_const_1 , partial)
			partial.pop()

	try:
		next(backtrack(_sage_const_0 , []))
	except StopIteration:
		pass

	return found_solution, total_attempts

# TO TEST TIMINGS
# n = 25
# FF = GF(16)
# matrices = generate_instance(FF, n)
# target = [FF.random_element() for i in range(n)]
# start = cputime()
# out = solve_triangular_matrices(FF, n, matrices, target)
# end = cputime()
# print(end-start)
# print(out[0])
# print(timeit('solve_triangular_matrices(FF, n, matrices, target)'))

def main():
	parser = argparse.ArgumentParser(
	)
	parser.add_argument(
		"-p", type=int, required=True,
		help="Number of trivial MQ problems."
	)
	parser.add_argument(
		"-q", type=int, required=True,
		help="Field size."
	)
	parser.add_argument(
		"-attempts", type=int, required=True,
		help="Number of experiments."
	)

	args = parser.parse_args()

	n = args.p
	FF = GF(args.q)
	solved = _sage_const_0 
	repeat = args.attempts
	total_attempts = _sage_const_0 
	for i in range(repeat):
		matrices = generate_instance(FF, n)
		target = [FF.random_element() for i in range(n)]
		out = solve_triangular_matrices(FF, n, matrices, target)
		if i == repeat//_sage_const_4 :
			print(f'1/4 is done!')
			print(f'Partial solution ratio:{(solved/(i+_sage_const_1 )).n()}')
			print(f'Partial visited nodes on average:{(total_attempts/(i+_sage_const_1 )).n()}')
			print('---------------------------------------------------------------')
		if i == repeat//_sage_const_2 :
			print(f'Half way there!')
			print(f'Partial solution ratio:{(solved/(i+_sage_const_1 )).n()}')
			print(f'Partial visited nodes on average:{(total_attempts/(i+_sage_const_1 )).n()}')
			print('---------------------------------------------------------------')
		if i == _sage_const_3 *repeat//_sage_const_4 :
			print(f'Almost there!')
			print(f'Partial solution ratio:{(solved/(i+_sage_const_1 )).n()}')
			print(f'Partial visited nodes on average:{(total_attempts/(i+_sage_const_1 )).n()}')
			print('---------------------------------------------------------------')
		if out[_sage_const_0 ] != None:
			solved += _sage_const_1 
		total_attempts += out[_sage_const_1 ]
	print(f'Solution ratio:{(solved/repeat).n()}')
	print(f'Visited nodes on average:{(total_attempts/repeat).n()}')



if __name__ == "__main__":
	main()

