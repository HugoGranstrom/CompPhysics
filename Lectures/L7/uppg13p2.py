import random as rd
import matplotlib.pyplot as plt
import math

def variance(f):
    N = len(f)
    return 1/N*(1/N*sum([x*x for x in f])-(1/N*sum(f))**2)

def f(x):
    return 1/(1+x*x)

def w1(x):
    return 1

def w2(x):
    return (4-2*x)/3

def x1(y):
    return y

def x2(y):
    return 2-math.sqrt(4-3*y)

def mcintegrate(f,w,xfunc,n):
    flist = []
    integlist = []
    varlist = []
    intsum = 0
    for i in range(n):
        y = rd.uniform(0,1)
        x = xfunc(y)
        fval = f(x)
        wval = w(x)
        flist.append(fval)
        intsum += fval/wval
        varlist.append(variance(flist))
        integlist.append(1/(i+1)*intsum)
    return flist,integlist,varlist

n = 10000
f1,int1,var1 = mcintegrate(f,w1,x1,n)
f2,int2,var2 = mcintegrate(f,w2,x2,n)
pi = [math.pi/4 for _ in range(n)]

plt.plot(pi[1000:],label="exact")
plt.plot(int1[1000:],label="1")
plt.plot(int2[1000:],label="2")
plt.legend()
plt.show()

plt.semilogy(var1)
plt.semilogy(var2)
plt.show()