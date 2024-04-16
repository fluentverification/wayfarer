# Wayfarer - A proof-of-concept seed-path (counterexample) heuristic

If successful, may be integrated into [STAMINA](https://github.com/fluentverification/stamina-storm) and/or [RAGTIMER](https://github.com/fluentverification/ragtimer)

Wayfarer treats a CRN as a VASS and computes a vectoral distance to the "boundary" of a desired counterexample.

## Methods to find probability bounds

### Lower probability boundary

#### Single-order priority

Single-order Priority assumes all solution states reside in a closed subspace, $S_s$, which in the RAGTIMER format, is always true. It then constructs a shortest distance to that closed subspace. In the RAGTIMER format, since the solution space is orthogonal to one or more species dimensions, there is no need to construct a projection matrix, just measure the distance on each non- "don't care" species. However, this method can be expanded to subspaces where a projection matrix ${P}_s$ can be constructed and prioritize states with low $\epsilon_s = |{P}_s{s}_{adj} - {s}_{adj}|_2$.

#### Iterative Subspace Reduction

Iterative subspace reduction takes the reactions in the dependency graph, and constructs a set of nested subspaces $S_0 \subseteq S_1 \subseteq S_2 \cdots \subseteq S_n$ where $S_n \cap S_s \neq \emptyset$. Note that $S_i = \{\ M_i{x} + {s}_0 + {f} | {x} \in {R}^m\ \}$ with ${f}$ ensuring that $S_n \cap S_s \neq \emptyset$. It then prioritizes states in smaller subspaces to approach $S_s$. These are calculated using a set of distances $\epsilon_i = |{P}_i{s}_{adj} - {s}_{adj}|_2$. ${s}_{adj} = {s} - ({s}_0 + {f})$ and ${P}_i = M_i(M_i^T M_i)^+ M_i^T$. Since $j < i \wedge \epsilon_i = 0 \implies \epsilon_j = 0$ we start calculating at $n$ and short circuit when zero distance is encountered.

From there, it's a priority first search, prioritizing first:

1. Higher indexes $i$ in which $\epsilon_i = 0$
2. Lower values of $\epsilon_{i + 1}$

It creates a partial state graph and seeks $K$ satisfying states.

## Installation

It is recommended to use a virtual environment to locally install the dependencies (since StormPy is a bit of a beast). In the project directory:

```sh
# Create a virtual environment
python3 -m venv .venv
# Activate the virtual environment
source .venv/bin/activate
# Install the dependencies
pip -r install dependencies.txt
```

When done, the `deactivate` command can be used.

## Usage

```bash
# For single order priority
./main.py -r $RAGTIMER_FILE -V -n $NUMDER_DESIRED_SATISFYING_STATES
# For iterative subspace reduction
./main.py -r $RAGTIMER_FILE -S -n $NUMDER_DESIRED_SATISFYING_STATES
```
