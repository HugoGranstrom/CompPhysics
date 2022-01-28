import numpy as np
import matplotlib.pyplot as plt

def numerov(a,b, y0, S, N):
  h = (b - a) / N
  r = a + h
  y = [y0[0], y0[1]]
  r = a + h

  for i in range(N-1):
    y.append(2*y[-1] - y[-2] + h*h / 12 * (S(r+h) + 10*S(r) + S(r-h)))
    r += h

  return y

N = 1000
a = 0
b = 30
h = (b - a) / N
fExact = lambda r: 1 - 0.5 * (r + 2) * np.exp(-r)
S = lambda r: -0.5 * r * np.exp(-r)

y0_exact = [0, fExact(a + h)]
y0_approx = [0, fExact(a + h)*1.05]
ySol = np.array(numerov(a, b, y0_exact, S, N))
ySolApprox = np.array(numerov(a, b, y0_approx, S, N))

x = np.linspace(a, b, N+1)
yExact = fExact(x)

k = (ySolApprox[-1] - ySolApprox[-2]) / (x[-1] - x[-2])

corrected = ySolApprox - k*x

plt.plot(x, ySolApprox, label="Approx y(h)")
plt.plot(x, ySol, label="Exact y(h)")
plt.plot(x, corrected, label="Corrected")
plt.legend()
plt.show()

plt.semilogy(x, np.abs(ySol - yExact), label="Exact y(h)")
plt.semilogy(x, np.abs(ySolApprox - yExact), label="Approx y(h)")
plt.semilogy(x, np.abs(corrected - yExact), label="Corrected")
plt.legend()
plt.show()

