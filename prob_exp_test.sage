import argparse 

def generate_instance(FF, n):
	matrix_list = [Matrix(FF, [FF(0) for i in range(n)]) for j in range(n)]
	for i in range(len(matrix_list)):
		for j in range(n - i - 1):
			matrix_list[i] = matrix_list[i].stack(vector([FF(0) for k in range(n)]))
		for j in range(n - i - 1, n):
			matrix_list[i] = matrix_list[i].stack(vector([FF(0) for k in range(n - i - 1)] + [FF.random_element() for k in range(i + 1)]))
	return [M[1:] for M in matrix_list]

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
) -> Tuple[Tuple[Any, ...] | None, int]:
	"""
	Solves a triangular system of quadratic equations incrementally and counts the number of attempted root lookups.
	Returns the first found solution (if any) and the number of attempts made.
	"""

	SR.<T> = PolynomialRing(FF)
	m = n
	found_solution = None
	total_attempts = 0

	def backtrack(level: int, partial: List[Any]) -> Generator[Tuple[Any, ...], None, None]:
		nonlocal total_attempts, found_solution

		if level == m:
			found_solution = tuple(partial)
			yield found_solution
			return

		R = polyring(FF, n - level)
		eq = matrix_to_eq(R, mats[level], tuple(partial))
		variables = [0 for _ in range(n - level - 1)] + [T]
		P = eq(variables) - target[n - level - 1]
		if show:
			print(P)
		if not show:
			if level == 0:
				P += FF.random_element() * T
		total_attempts += 1

		if P != 0:
			roots = P.roots(multiplicities=False)
		else:
			roots = list(FF)

		for r in roots:
			partial = [FF(r)] + partial
			if show:
				print(partial)
			yield from backtrack(level + 1, partial)
			partial.pop()

	try:
		next(backtrack(0, []))
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
	solved = 0
	repeat = args.attempts
	total_attempts = 0
	for i in range(repeat):
		matrices = generate_instance(FF, n)
		target = [FF.random_element() for i in range(n)]
		out = solve_triangular_matrices(FF, n, matrices, target)
		if i == repeat//4:
			print(f'1/4 is done!')
			print(f'Partial solution ratio:{(solved/(i+1)).n()}')
			print(f'Partial visited nodes on average:{(total_attempts/(i+1)).n()}')
			print('---------------------------------------------------------------')
		if i == repeat//2:
			print(f'Half way there!')
			print(f'Partial solution ratio:{(solved/(i+1)).n()}')
			print(f'Partial visited nodes on average:{(total_attempts/(i+1)).n()}')
			print('---------------------------------------------------------------')
		if i == 3*repeat//4:
			print(f'Almost there!')
			print(f'Partial solution ratio:{(solved/(i+1)).n()}')
			print(f'Partial visited nodes on average:{(total_attempts/(i+1)).n()}')
			print('---------------------------------------------------------------')
		if out[0] != None:
			solved += 1
		total_attempts += out[1]
	print(f'Solution ratio:{(solved/repeat).n()}')
	print(f'Visited nodes on average:{(total_attempts/repeat).n()}')



if __name__ == "__main__":
	main()
