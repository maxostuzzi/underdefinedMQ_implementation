import argparse

def filter_test(n, m, partition, treshold):
	k_1 = partition[-2] # aka r
	k_2 = partition[-1] # aka l
	p = len(partition) - 2
	good = True
	assert sum(partition) == m, f'The partition sums up to {sum(partition)}, while m = {m}'
	print(f'It is a {p}-partition and {k_1+k_2} coordinates are guessed.')
	
	# (i-1) * (sum_{1 <= j <= p} b_j + k_2 - (i-1)) + i <= n    for 1 <= i <= b_1
	for i in range(1, partition[0]):
		good = good and (((i * (sum(partition[:-2]) + k_2 - i)) + i) <= n)
		if not (((i * (sum(partition[:-2]) + k_2 - i)) + i) <= n):
			print(f'\x1b[31mFirst block, column {i+1} CANNOT be transformed\x1b[0m')
			print(f'{i * (sum(partition[:-2]) + k_2 - i) + i} > {n}')
		elif (((i * (sum(partition[:-2]) + k_2 - i)) + i) <= n):
			print(f'\x1b[32mFirst block, column {i+1} can be transformed\x1b[0m')
	if p > 2:
		# (b_1 + (i-1)) * (sum_{2 <= j  <= p} b_j - (i-1)) + i * b_1 + i <= n   for 1 <= i <= b_2
		for i in range(0, partition[1]):
			good = good and ((((partition[0] + i ) * (sum(partition[1:-2]) - i)) + (i * partition[0]) + i) <= n)
			if not ((((partition[0] + i ) * (sum(partition[1:-2]) - i)) + (i * partition[0]) + i) <= n):
				print(f'\x1b[31mSecond block, column {i+1} CANNOT be transformed\x1b[0m')
			elif ((((partition[0] + i ) * (sum(partition[1:-2]) - i)) + (i * partition[0]) + i) <= n):
				print(f'\x1b[32mSecond block, column {i+1} can be transformed\x1b[0m')

		for l in range(2,p-1):              # 3 <= l <= p-1
		
			sum_1 = sum(partition[:l-1])    # sum_{1 <= j = l-1} b_j
		
			sum_2 = sum(partition[l-1:p])   # sum_{l <= j <= p} b_j
			
			sum_3 = 0        
			for j in range(l-1):            # sum_{1 <= j < h <= l - 1}
				for h in range(j+1, l-1):
					sum_3 = sum_3 + partition[j] * partition[h]
		
			for i in range(0, partition[l]):
				good = good and (((( sum_1 + i ) * ( sum_2 - i )) + (i * sum_1) +  sum_3 + i) <= n)
				if not (((( sum_1 + i ) * ( sum_2 - i )) + (i * sum_1) +  sum_3 + i) <= n):
					print(f'\x1b[31m{l+1}-th block, column {i} CANNOT be transformed\x1b[0m')
				elif (((( sum_1 + i ) * ( sum_2 - i )) + (i * sum_1) +  sum_3 + i) <= n):
					print(f'\x1b[32m{l+1}-th block, column {i} can be transformed\x1b[0m')

	sum_4 = 0
	for i in range(p):
		for j in range(i+1,p):
			# print(f'b{i+1} = {partition[i]}')
			# print(j)
			# print(f'b{j+1} = {partition[j]}')
			# print('===========================================')
			sum_4 = sum_4 + partition[i] * partition[j] 
	good = good and ((sum_4 + m - partition[-1] - partition[-2])  <= n)
	# print(sum_4 + m - partition[-1] - partition[2])
	# print(n)
	if not ((sum_4 + m - partition[-1] - partition[-2])  <= n):
		print(f'\x1b[31mLast block CANNOT be transformed\x1b[0m')
	elif ((sum_4 + m - partition[-1] - partition[-2])  <= n):
		print(f'\x1b[32mLast block can be transformed\x1b[0m')
	good = ((sum_4 + m)  <= n)
	if not good:
		print(f'\x1b[31mMQ System CANNOT be solved\x1b[0m')
		print(f'{sum_4 + m} >= {n}')
	elif good:
		print(f'\x1b[32mMQ System can be solved\x1b[0m')



if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Optimization Algorithm")
	parser.add_argument("-n", type=int, help="Number of variables")
	parser.add_argument("-m", type=int, help="Number of equations")
	parser.add_argument("-partition", nargs='+', type=int, help="Partition")
	parser.add_argument("-k", nargs='+', type=int, help="Guessed")
	parser.add_argument("-t", type=int, default=999, help="Defines a treshold to reduce the number of combinations")
	
	args = parser.parse_args()
	filter_test(args.n, args.m, args.partition + args.k, args.t)
