# Wayfarer - A proof-of-concept seed-path (counterexample) heuristic

If successful, may be integrated into [STAMINA](https://github.com/fluentverification/stamina-storm) and/or [RAGTIMER](https://github.com/fluentverification/ragtimer)

Wayfarer treats a CRN as a VASS and computes a vectoral distance to the "boundary" of a desired counterexample. It can then prioritize on a number of things:

- Just distance (implemented, works well)
- Angle between distance and probability flow (attempted, doesn't seem to work well)
- Some combination of the two
