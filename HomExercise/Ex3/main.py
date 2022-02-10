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

def gaussQuad(f, a, b, l):
  abscissae = rootsFinder(lambda x: legendre(x, l), -1, 1, 1e-15)
  weights = [2 / ((1 - x*x) * legendreDeriv(x, l)**2) for x in abscissae]
  new_f = lambda t: (b - a) / 2 * f((t + 1)*(b-a)/2 + a)
  res = 0
  for x, w in zip(abscissae, weights):
    res += w*new_f(x)

  return res

def simps(a,b,f,N):
    h = (b-a)/N
    res = 0
    x = a+h
    for i in range(N//2):
        res += (h/3)*(f(x-h)+4*f(x)+f(x+h))
        x += 2*h
    return res

def bode(a,b,f,N):
    h = (b-a)/N
    res = 0
    x = a
    for i in range(N//4):
        res += ((2*h)/45)*(7*f(x)+32*f(x+h)+12*f(x+2*h)+32*f(x+3*h)+7*f(x+4*h))
        x += 4*h
    return res

def f(x):
  return (1 - x*x) ** 0.5

correct = math.pi / 2
a = -1
b = 1

ls = list(range(3, 16))
gaussErrors = []

for l in ls:
  I = gaussQuad(f, a, b, l)
  err = abs(correct - I)
  gaussErrors.append(err)
  print(f"l = {l}: {err}")

fevalsSimp = []
simpErrors = []
for i in [1, 2, 3, 4, 5]:
  simpErrors.append(abs(correct - simps(a, b, f, i*2)))
  fevalsSimp.append(i*3)

fevalsBode = []
bodeErrors = []
for i in [1, 2, 3]:
  bodeErrors.append(abs(correct - bode(a, b, f, 4*i)))
  fevalsBode.append(i*5)

Isimp = simps(a, b, f, 2)
Ibode = bode(a, b, f, 4)
print("simpson:", abs(correct - Isimp))
print("bode:", abs(correct - Ibode))

plt.figure()
plt.semilogy(ls, gaussErrors, label="Gauss")
plt.semilogy(fevalsSimp, simpErrors, label="Simpson")
plt.semilogy(fevalsBode, bodeErrors, label="Bode")
plt.legend()
plt.xlabel("Function evaluations")
plt.ylabel("Error")
plt.title("Comparisons between Gauss, Simpson and Bode")
plt.savefig("comparison.png")

plt.figure()
x = np.linspace(a, b, 100)
y = f(x)
plt.plot(x, y)
plt.title("Plot of f(x)")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.savefig("f(x).png")





