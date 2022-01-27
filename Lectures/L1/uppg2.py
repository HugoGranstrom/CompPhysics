import math
import matplotlib.pyplot as plt

def trap(a,b,f,N):
    h = (b-a)/N
    res = 0
    x = a+h
    for i in range(N//2):
        res += (h/2)*(f(x-h)+2*f(x)+f(x+h))
        x += 2*h
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

a = 0
b = 1
f = math.exp
exact = 1.718282

N = [4, 8, 12, 16, 20, 24, 28]
restrap = []
ressimps = []
resbode = []
for i in N:
    restrap.append(abs(trap(a,b,f,i)-exact))
    ressimps.append(abs(simps(a,b,f,i)-exact))
    resbode.append(abs(bode(a,b,f,i)-exact))

print(restrap)
print(ressimps)

plt.loglog(N,restrap)
plt.loglog(N,ressimps)
plt.loglog(N,resbode)
plt.legend(["Trapz", "Simpson","Bode"])
plt.show()