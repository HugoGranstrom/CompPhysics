import math
import numpy as np
import matplotlib.pyplot as plt

def bode(a,b,f,N):
    h = (b-a)/N
    res = 0
    x = a
    for i in range(N//4):
        res += ((2*h)/45)*(7*f(x)+32*f(x+h)+12*f(x+2*h)+32*f(x+3*h)+7*f(x+4*h))
        x += 4*h
    return res

a = 4
r_max = 30
N = 500

def phi_low(r):
  return np.sqrt(1/(2*a)) * (np.exp(a*r) - np.exp(-a*r))

def phi_high(r):
  return -np.sqrt(1/(2*a)) * np.exp(-a*r)

def phi_exact(r):
  return (1 / (1 - a*a))**2 * (np.exp(-a*r) - np.exp(-r) * (1 + 0.5 * (1 - a*a) * r))

S = lambda r: -0.5 * r * np.exp(-r)

def f_low(r):
  return phi_low(r) * S(r)

def f_high(r):
  return phi_high(r) * S(r)

def phi(r):
  return phi_high(r) * bode(0, r, f_low, N) + phi_low(r) * bode(r, r_max, f_high, N)

t = np.linspace(0, r_max, 200)

y = np.array([phi(t[i]) for i in range(len(t))])
yExact = phi_exact(t)

plt.figure()
plt.plot(t, y, label="Numerical")
plt.plot(t, yExact, label="Exact", linestyle="dashed")
plt.legend()
plt.title(f"Comparison of numerical and exact solution (N = {N})")
plt.xlabel("r")
plt.ylabel("φ(r)")
plt.savefig("phiSol.png")

plt.figure()
upper_limit = int(200 / r_max * 3)
plt.plot(t[:upper_limit], y[:upper_limit], label="Numerical")
plt.plot(t[:upper_limit], yExact[:upper_limit], label="Exact", linestyle="dashed")
plt.legend()
plt.title(f"Comparison of numerical and exact solution (N = {N})")
plt.xlabel("r")
plt.ylabel("φ(r)")
plt.savefig("phisolZoom.png")

plt.figure()
plt.semilogy(t, np.abs(y - yExact))
plt.xlabel("r")
plt.ylabel("log(|φ - φ_exact|)")
plt.title(f"Error of numerical solution (N = {N})")
plt.savefig("phiError.png")