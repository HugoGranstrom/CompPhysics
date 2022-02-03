import math
import numpy as np
import matplotlib.pyplot as plt

def legendre(x, l):
  p = [0*x + 1, x]
  for i in range(1, l):
    p.append(((2*i + 1)*x*p[i] - i*p[i-1]) / (i + 1))
  return p[l]

def legendreDeriv(x, l):
  if l == 0: return x*0
  return (-l*x*legendre(x, l) + l*legendre(x, l-1)) / (1 - x*x)

def bisection(f, x0, x1, tol):
  a = x0
  b = x1
  fa = f(a)
  fb = f(b)
  assert fa * fb < 0
  while abs(b - a) > tol:
    c = (a + b) / 2
    fc = f(c)
    if fc*fa < 0:
      b = c
      fb = fc
    elif fc*fb < 0:
      a = c
      fa = fc
    elif fc == 0:
      return fc
    else:
      print("Something went wrong in bisection!")
  return (a + b) / 2

def rootsFinder(f, a, b, tol):
  intervals = []
  currentValue = f(a)
  step = 0.01
  for x in np.arange(a, b, step):
    temp = f(x)
    if temp*currentValue < 0:
      intervals.append((x-step, x))
    currentValue = temp
  roots = [bisection(f, interval[0], interval[1], tol) for interval in intervals]
  return roots

if __name__ == "__main__":
  for l in range(3, 11):
    f = lambda x: legendre(x, l)
    roots = rootsFinder(f, -1, 1, 1e-8)
    assert len(roots) == l
    print(f"Roots of P_{l} are: {roots}")

  x = np.linspace(-1, 1, 50)
  plt.figure()
  for l in range(0, 6):
    plt.plot(x, legendre(x, l), label=f"l = {l}")
  plt.title("Legendre polynomials of degree l")
  plt.legend()
  plt.show()

  plt.figure()
  for l in range(0, 6):
    plt.plot(x, legendreDeriv(x, l), label=f"l = {l}")
  plt.title("Legendre polynomial derivatives of degree l")
  plt.legend()
  plt.show()

