# Wayfarer - A proof-of-concept seed-path (counterexample) heuristic

If successful, may be integrated into [STAMINA](https://github.com/fluentverification/stamina-storm) and/or [RAGTIMER](https://github.com/fluentverification/ragtimer)

Wayfarer treats a CRN as a VASS and computes a vectoral distance to the "boundary" of a desired counterexample.

## Methods to find probability bounds

### Lower probability boundary

#### Single-order priority

Single-order Priority assumes all solution states reside in a closed subspace, $S_s$, which in the RAGTIMER format, is always true. It then constructs a shortest distance to that closed subspace. In the RAGTIMER format, since the solution space is orthogonal to one or more species dimensions, there is no need to construct a projection matrix, just measure the distance on each non- "don't care" species. However, this method can be expanded to subspaces where a projection matrix $\mathcal{P}_s$ can be constructed.

#### Iterative Subspace Reduction

Iterative subspace reduction takes the reactions in the dependency graph, and constructs a set of nested subspaces $S_0 \subseteq S_1 \subseteq S_2 \cdots \subseteq S_n$ where $S_n \cap S_s \neq \emptyset$. It then prioritizes states in smaller subspaces to approach $S_s$.
