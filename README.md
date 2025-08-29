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

Note that the `guessed` and the `trivial` must satisfy the constraints in the paper.

`transformation.sage` then generates a random MQ map with $P:\mathbb{F}^n_q\rightarrow\mathbb{F}^m_q$ described by a list of matrices and then computes the change of variables $S:\mathbb{F}^n_q\rightarrow\mathbb{F}^n_q$ and the system transformation $T:\mathbb{F}^m_q\rightarrow\mathbb{F}^m_q$ as described in Section 3 of the paper.
It will print all the matrices of the new MQ map $\tilde P = T\circ P\circ S$.

An example:

- `./sage ./transformation.sage -n 100 -m 14 -q 16 -k 10 -p 2`
  
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
- `p` the number of MQ$(1,1)$;
- `attempts` the number of experiments we want to conduct,
and returns the experimental probability $s_p$ and the experimental expectation $e_p$.

To run the algorithm from the terminal, run for example

`python prob_exp_test.py -q 16 -p 25 -attempts 10000`.

***

## Optimization Algorithm

Our optimization algorithm `optimization_classical_and_quantum.py` takes as input the following parameters:
- `n` number of variables;
- `m` number of equations;
- `q` the characteristic of the field $\mathbb{F}_q$;
- `p` a partition number;
- `c` number of cpus for running in parallel.

The algorithm generates all possible partition of the input `m`, filters them according to the constraints described in Section 3 and 4 of the paper, evaluates the corresponding bit complexity (without accounting for the failure probability) and ouputs a list of partitions sorted starting from the most efficient one. The results are displayed on the terminal window, as shown below.

![](example_output_optimization.png)

### Prerequisites

Our algorithm `optimization_classical_and_quantum.py` runs with **Python 3.13.2**. For our algorithm the libraries `tqdm` and `cryptographic_estimators` are additionally required. These can be installed for example using pip via the following two commands:
- `pip install tqdm`
- `pip install cryptographic_estimators`. 

For a more detailed description of the libraries and customized installations please follow the corresponding references:
- `tqdm` documentation [TQDM](https://tqdm.github.io/)
- `cryptographic_estimators` documentation [Crypto-TII](https://github.com/Crypto-TII/CryptographicEstimators)

An **example** is provided below:

- `python .\optimization_classical_and_quantum.py -n 1680 -m 190 -q 7 -p 8 -c 14`

### Optional Parameters

The algorithm can take as input the following optional parameters:
- `t` threshold that filters out partitions for which the exhaustive search bit complexity superseeds $t$;
- `maxsummand` i.e $b_{2},...,b_{p} \leq \text{maxsummand}$
- `maxb1` i.e $b_1 \leq \text{maxb1}$
- `minsumguess` $k \gt \text{minsumguess}$
- `lines` number of results in output file

  ***
