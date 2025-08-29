# Just Guess!
This is the implementation of the algorithms described in the paper Just Guess: Improved Algorithm for the Underdetermined MQ Problem.
***

## Transformation Algorithm

Our algorithm `transformation.sage` takes as input:
- `n` number of variables;
- `m` number of equations;
- `q` the characteristic of the field $\mathbb{F}_q$;
- `k` number of guessed coordinates;
- `p` the number of MQ $(1,1)$. 

Note that the `k` and the `p` must satisfy the constraints in the paper.

`transformation.sage` then generates a random MQ map with $P:\mathbb{F}^n_q\rightarrow\mathbb{F}^m_q$ described by a list of matrices and then computes the change of variables $S:\mathbb{F}^n_q\rightarrow\mathbb{F}^n_q$ and the system transformation $T:\mathbb{F}^m_q\rightarrow\mathbb{F}^m_q$ as described in Section 3 of the paper.
It will print all the matrices of the new MQ map $\tilde P = T\circ P\circ S$.

An example:

- `./sage ./transformation.sage -n 11 -m 6 -q 16 -k 2 -p 2`
  
***

## Precise Complexity Algorithm

Our algorithm `precise_complexity.py` takes as input the following parameters:
- `q` the characteristic of the field $\mathbb{F}_q$;
- `p` a partition number, which ultimately corresponds to the number of trivial MQ problems;
- `k` two integers, which correspond to $k_1$ and $k_2$ in the paper;
- `grover` an optional toggle input, which triggers quantum complexity instead of the classical one.

The algorithm computes the precise complexity relative to the input in the case we enter the polynomial regime, as described in the paper. It outputs the success probability $s_p$, the expected number of trivial MQ problem solved $e_p$ (both computed using the recursive formulas) and the precise complexity.
To run the algorithm from the terminal, run for example

`python precise_complexity.py -q 16 -p 25 -k 15 19 --grover`.

***

## Test Probability and Expectation

This algorithm, called `prob_exp_test.py` takes as input the following parameters:
- `q` the characteristic of the field $\mathbb{F}_q$;
- `p` the number of MQ $(1,1)$;
- `attempts` the number of experiments we want to conduct,
and returns the experimental probability $s_p$ and the experimental expectation $e_p$.

To run the algorithm from the terminal, run for example

`python prob_exp_test.py -q 16 -p 25 -attempts 10000`.

***

## Optimization Algorithm

Our optimization algorithm `optimizationm.py` takes as input the following parameters:
- `n` number of variables;
- `m` number of equations;
- `q` the characteristic of the field $\mathbb{F}_q$;

The algorithm exhaustively searches for the optimal choice of the parameters $(k,p)$. Then, it prints
- the optimal parameters found excluding polynomial factors;
- the optimal parameters found including polynomial factors;
- the probability $s_p$;
- the number of expected nodes in the MQ tree $e_p$ visited by the DFS;
- the complexity without polynomial factors;
- the complexity with polynomial factors;
both classically and quantumly.

![](example_output_optimization.png)

### Prerequisites

For our algorithm the library `njit` is additionally required.

- `python .\optimization.py 860 78 16`
  ***
