# Wayfarer - A proof-of-concept seed-path (counterexample) heuristic

If successful, may be integrated into [STAMINA](https://github.com/fluentverification/stamina-storm) and/or [RAGTIMER](https://github.com/fluentverification/ragtimer)

Wayfarer treats a CRN as a VASS and computes a vectoral distance to the "boundary" of a desired counterexample.

## Methods to find probability bounds

### Lower probability boundary

#### Single-order priority

These methods explore a partial state space using a priority queue to explore states with a certain priority metric. When a satisfying state is found, a "traceback" method is used to find a certain  number of traces back to the init state. The main difference is the priority metric used. The following are things which can be incormporated into the priority.

- **Vectoral Distance to Boundary**: takes the shortest vector to the boundary of satisfying states and computes the L2 norm of that. This distance, which we want to minimize, is the priority.
- **Angles**: There are a number of angle elements that can be included:
	+ **Angle between flow and transition vector:** We create a flow vector based on the available transitions, weighted by their rates. We wish to minimize the angle between the flow and transition vector, however, we can just maximize the transition rate.
	+ **Angle between transition vector and vector to boundary:** Are we going in the right direction?

#### Multiple-order priority

*Not yet implemented*

This method makes use of Landon's dependency graph, but takes the reactions from the root node and assigns them an "order", which just corresponds to the level in the graph, from the bottom. The root node is order 1, and every reaction beyond that is order 2, 3, etc. This gives us an order of reaction vectors:

R1, R2, R3, ..., RN

This ordering is crucial, as we then take the reaction vectors in this order and create subspaces with these vectors as basis vectors. If a particular state, shifted by the init state vector, is not in this space, it cannot reach the satisfying state, and is not explored. Obviously, we can restrict reactions to those in the graph, and then this is not necessary. However, we then create a subspace:

span(R2, R3, ..., RN)

and prioritize low distances to that subspace. Once in that subspace, we create another subspace

span(R3, R4, ..., RN)

And prioritize on low distances to that subspace. We do this until we are within a subspace of rank 1, which we then just prioritize on distance to the satisfying state.

### Upper bound probability

-
