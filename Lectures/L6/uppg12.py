import timeit
import numpy as np
import matplotlib.pyplot as plt

def variational(S,w,BC,phi0,N,tol,Efunc):

    h = 1/N
    phi = phi0
    phi[0] = BC[0]
    phi[-1] = BC[1]
    Ediff = tol+1
    E = Efunc(phi,S,N)
    Elist = [E]
    while tol < Ediff:
        for i in range(1, N-1):
            phi[i] = (1-w)*phi[i] + w/2*(phi[i+1]+phi[i-1]+h*h*S(i*h))
        newE = Efunc(phi,S,N)
        Elist.append(newE)
        Ediff = abs(E-newE)
        E = newE
    return Elist,phi

def Efunc(phi,S,N):
    h = 1/N
    res = 0
    for i in range(1,N):
        res += 1/(2*h)*(phi[i]-phi[i-1])**2-h*S(i*h)*phi[i]
    return res

def S(x):
    return 12*x*x

def interp(phi):
    output = []
    for i in range(len(phi)-1):
        output.append(phi[i])
        output.append((phi[i]+phi[i+1])/2)
    output.append(phi[-1])
    return np.array(output)

def single(N):
    BC = [0, 0]
    w = 1.9
    tol = 1e-8
    phi = np.zeros((N))
    Energy,Sol = variational(S,w,BC,phi,N,tol,Efunc)
    return Energy

def interpolated(N):
    BC = [0, 0]
    w = 1.9
    tol = 1e-8
    phi = np.zeros((N))
    Energy,Sol = variational(S,w,BC,phi,N,tol,Efunc)
    solinterp = interp(Sol)
    Energyrefined, Solrefined = variational(S,w,BC,solinterp,N*2-1,tol,Efunc)
    return Energyrefined

N = 100

t1 = timeit.timeit("single(2*N)",setup = "N = 100",number=1,globals=globals())
t2 = timeit.timeit("interpolated(N)",setup ="N = 100",number=1,globals=globals())
E = single(2*N)
Eref = interpolated(N)

print(t1)
print(t2)
print(E[-1])
print(Eref[-1])

""" x = np.linspace(0,1,N)
Solinterpt = interp(Sol)
newx = interp(x)
plt.plot(newx,Solinterpt,label="Interpolated")
plt.plot(x,Sol,Label="Solution")
plt.legend()
plt.show()
plt.plot(Energy)
plt.plot(np.ones(len(Energy))*(-9/14))
plt.show()
EnergyRefined,SolRefined = variational(S,w,BC,Solinterpt,N*2-1,tol,Efunc)

plt.plot(newx,SolRefined)
plt.plot(x,Sol)
plt.show()
print(f"Energy = {Energy[-1]}, Energy Refined = {EnergyRefined[-1]}") """