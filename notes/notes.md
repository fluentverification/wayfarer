# Notes

- Using just distance as priority, on our toy model (`test1.py`), we find estimated lower bound of 3.3845545545021705e-15, whereas random exploration gets values such as 3.618909437343446e-75, 4.570745425741801e-47, 2.140080580527486e-59, orders of magnitude closer to 0.
- Again using just distance as priority, we get about 0.007-0.008 s on guided method, but 0.001-0.008 s on random exploration. Although longer runtime, the guided method produces a much more accurate bound.
