import math
from uppg10 import *

def gaussQuad(f, a, b, l):
  abscissae = rootsFinder(lambda x: legendre(x, l), -1, 1, 1e-15)
  weights = [2 / ((1 - x*x) * legendreDeriv(x, l)**2) for x in abscissae]
  print("absc:", abscissae)
  print("weights:", weights)
  new_f = lambda t: (b - a) / 2 * f((t + 1)*(b-a)/2 + a)
  res = 0
  for x, w in zip(abscissae, weights):
    res += w*new_f(x)

  return res

def f(x):
  print(x)
  return math.sqrt(1 + x)

print(gaussQuad(f, 0, 3, 3))
