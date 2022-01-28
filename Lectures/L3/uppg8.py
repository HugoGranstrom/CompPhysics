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
a = 100
b = 0
h = (b - a) / N
fExact = lambda r: 1 - 0.5 * (r + 2) * np.exp(-r)
S = lambda r: -0.5 * r * np.exp(-r)

y0_exact = [fExact(a), fExact(a + h)]
y0_approx = [fExact(a), fExact(a + h)*0.95]
ySol = np.array(numerov(a, b, y0_exact, S, N))
ySolApprox = np.array(numerov(a, b, y0_approx, S, N))

x = np.linspace(a, b, N+1)

k = (ySolApprox[3] - ySolApprox[4]) / (x[3] - x[4])

corrected = ySolApprox - k * x - ySolApprox[-1]

yExact = fExact(x)

print(f"Error at end: {abs(ySol[-1])} and {abs(ySolApprox[-1])} and {abs(corrected[-1])}")

plt.plot(x, ySol)
#plt.plot(x, ySolApprox)
plt.plot(x, corrected)
plt.legend(["Exact", "Approx", "Corrected"])
plt.show()

plt.semilogy(x, np.abs(ySol - yExact))
plt.semilogy(x, np.abs(ySolApprox - yExact))
plt.semilogy(x, np.abs(corrected - yExact))
plt.legend(["Exact", "Approx", "Corrected"])
plt.show()