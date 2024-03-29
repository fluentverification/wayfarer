# Notes

- Using just distance as priority, on our toy model (`test1.py`), we find estimated lower bound of 3.3845545545021705e-15, whereas random exploration gets values such as 3.618909437343446e-75, 4.570745425741801e-47, 2.140080580527486e-59, orders of magnitude closer to 0.
- Again using just distance as priority, we get about 0.007-0.008 s on guided method, but 0.001-0.008 s on random exploration. Although longer runtime, the guided method produces a much more accurate bound.
- We really only find one satisfying state, and then traceback gives us a number of high-probability counterexamples
- For the second (more challenging test) wayfarer finds a lower bound of 6.852600255650947e-75 in 0.016871929168701172 s, whereas random exploration did not immediately terminate.

## Ideas to prioritize on:

- Angle to boundary (are we moving in the right direction?)
- Weighted average of angle to boundary and angle between flow (are we moving in the most probable direction?)

## Things to make our lives easier:

- Use PyPy for speed
