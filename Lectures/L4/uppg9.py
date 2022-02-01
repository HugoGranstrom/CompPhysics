import matplotlib.pyplot as plt
import numpy as np

def search(f,x0,tol,stepsize):
    x = x0
    xlist = [x]
    ylist = [f(x)]
    val1 = f(x)
    step = stepsize
    while abs(step) > tol:
        x = x + step
        val2 = f(x)
        if val1*val2 < 0:
            step /= -2
        val1 = val2
        xlist.append(x)
        ylist.append(f(x))
    return xlist,ylist

def secant(f,x0, x1 ,tol):
    x = x1
    xprev = x0
    diff = abs(x-xprev)
    y = f(x1)
    yprev = f(x0)
    xlist = [x0, x1]
    ylist = [f(x0),f(x1)]
    while diff > tol:
        xtemp = x - y*(x-xprev)/(y-yprev)
        xprev = x
        x = xtemp
        yprev = y
        y = f(x)
        diff = abs(x-xprev)
        xlist.append(x)
        ylist.append(y)
    return xlist,ylist

def numerov(a,b, y0, k, N):
  h = (b - a) / N
  r = a + h
  y = [y0[0], y0[1]]
  r = a + h

  for i in range(N-1):
    y.append(2 * (1 - 5*h*h/12*k*k)/(1 + h*h/12*k*k) * y[-1] - y[-2])
    r += h

  return y

def f(k):
    return numerov(0, 1, [0, 1e-4], k, 1000)[-1]

k = 1
N = 1000
xs, ys = search(f, k, 1e-8, 1)
print(len(xs))
print(xs[-1])

sol = numerov(0, 1, [0, 1e-4], xs[-1], N)
t = np.linspace(0, 1, N+1)
plt.plot(t, sol, label="Numeric")
plt.legend()
plt.show()



