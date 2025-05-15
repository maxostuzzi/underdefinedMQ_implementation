# Just Guess!
This is the implementation of the algorithms described in the paper Just Guess: Improved Algorithm for the Underdetermined MQ Problem.

***

## Filtering Algorithm

Our filtering algorithm `filter_test.py` takes as input:
- `n` number of variables;
- `m` number of equations;
- `guessed` a partition of $k$ as described in the paper;
- `partition` a proper partition of $m-k$.

It will check whether `guessed` and `partition` satisfy the constraints from the paper and it will print a detailed overview of what steps can be done and what steps cannot.

An example:

`python filter_test.py -n 1848 -m 142 -guessed 1 1 -partition 9 1 1 1 1 127`

***

## Transformation Algorithm

Our algorithm `transformation.sage` takes as input:
- `n` number of variables;
- `m` number of equations;
- `q` the characteristic of the field $\mathbb{F}_q$;
- `guessed` a partition of $k$ as described in the paper;
- `partition` a proper partition of $m-k$. 

Note that the `guessed` and the `partition` must satisfy the constraints in the paper.
To check whether a `guessed` and `partition` are valid, one can use the filtering algorithm `full_filter_test.py`, see the dedicated description above.

`transformation.sage` then generates a random MQ map with $P:\mathbb{F}^n_q\rightarrow\mathbb{F}^m_q$ described by a list of matrices and then computes the change of variables $S:\mathbb{F}^n_q\rightarrow\mathbb{F}^n_q$ and the system transformation $T:\mathbb{F}^m_q\rightarrow\mathbb{F}^m_q$ as described in Section 3 of the paper.
It will print all the matrices of the new MQ map $\tilde P = T\circ P\circ S$.

An example:

- `python ./transformation.sage -n 100 -m 14 -q 16 -guessed 1 1 -partiton 4 4 4`
  
***

## Optimization Algorithm

Our optimization algorithm `optimization_classical_and_quantum.py` takes as input the following parameters:
- `n` number of variables;
- `m` number of equations;
- `q` the characteristic of the field $\mathbb{F}_q$;
- `p` a partition number;
- `c` number of cpus for running in parallel.

The algorithm generates all possible partition of the input `m`, filters them according to the constraints described in Section 3 and 4 of the paper, evaluates the corresponding bit complexity and ouputs a list of partitions sorted starting from the most efficient one. The results are displayed on the terminal window, as shown below.

![](https://anonymous.4open.science/r/underdefinedMQ_implementation-6A0E/example_output_optimization.png)

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
