import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

import numpy as np
from fem.gauss import gaussPoints

# 10.15(a)
I = 0

p, w = gaussPoints(n = 2)

for xi, wi in zip(p, w):
    I += np.cos(xi/2)*wi

print(I)

# 10.15(b)
I = 0

p, w = gaussPoints(n = 2)

for xi, wi in zip(p, w):
    I += xi**2*wi

print(I)

# 10.15(c)
I = 0

p, w = gaussPoints(n = 3)

for xi, wi in zip(p, w):
    I += xi**4*wi

print(I)

# 10.15(d)
I = 0

p, w = gaussPoints(n = 2)

for xi, wi in zip(p, w):
    I += np.cos(xi)/(1 - xi**2)*wi

print(I)

# 10.15(e)
I = 0

p, w = gaussPoints(n = 2)

for xi, wi in zip(p, w):
    I += xi**3*wi

print(I)

# 10.15(f)
I = 0

p, w = gaussPoints(n = 3)

for xi, wi in zip(p, w):
    I += xi*np.cos(xi)*wi

print(I)

# 10.15(g)
I = 0

p, w = gaussPoints(n = 3)

for xi, wi in zip(p, w):
    I += (4**xi - 2*xi)*wi

print(I)