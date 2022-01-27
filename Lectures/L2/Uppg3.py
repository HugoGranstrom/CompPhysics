import math
import matplotlib.pyplot as plt

def euler(f,a,b,y0,N):
    h = (b-a)/N
    x = a
    y = y0
    for i in range(N):
        y += h*f(x,y)
        x += h
    return y

def taylor(f,dfx,dfy,a,b,y0,N):
    h = (b-a)/N
    x = a
    y = y0
    for i in range(N):
        y += h*f(x,y) + (1/2)*h**2*(dfx(x,y)+f(x,y)*dfy(x,y))
        x += h
    return y

def imp(f,g,a,b,y0,N):
    h = (b-a)/N
    x = a
    y = y0
    for i in range(N):
        y = (1+0.5*g(x)*h)/(1-0.5*g(x+h)*h)*y
        x += h
    return y

f = lambda x,y:-x*y
dfx = lambda x,y:-y
dfy = lambda x,y:-x
g = lambda x:-x
a = 0
b = 3
N = [2**k for k in range(10)]
y0 = 1
fexact = lambda x:math.exp(-x**2/2)
erreuler = [abs(euler(f,a,b,y0,i)-fexact(b)) for i in N]
errtay = [abs(taylor(f,dfx,dfy,a,b,y0,i)-fexact(b)) for i in N]
errimp = [abs(imp(f,g,a,b,y0,i)-fexact(b)) for i in N]

""" plt.loglog(N,erreuler)
plt.loglog(N,errtay)
plt.loglog(N,errimp)
plt.legend(["Euler","Taylor","Implicit"])
plt.show() """

euy0 = [euler(f,a,b,y0,i) for i in N]
tayy0 = [taylor(f,dfx,dfy,a,b,y0,i) for i in N]
impy0 = [imp(f,g,a,b,y0,i) for i in N]
erreu2 = [abs(euler(f,b,a,euy0[i],N[i])-y0) for i in range(len(N))]
errtay2 = [abs(taylor(f,dfx,dfy,b,a,euy0[i],N[i])-y0) for i in range(len(N))]
errimp2 = [abs(imp(f,g,b,a,euy0[i],N[i])-y0) for i in range(len(N))]

plt.loglog(N,erreu2)
plt.loglog(N,errtay2)
plt.loglog(N,errimp2)
plt.legend(["Euler","Taylor","Implicit"])
plt.show()