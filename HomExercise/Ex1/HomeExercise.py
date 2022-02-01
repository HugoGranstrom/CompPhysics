import matplotlib.pyplot as plt
import numpy as np

def rk2(f,a,b,y0,N):
    y = np.zeros((len(y0),N))
    y[:,0] = y0
    h = (b-a)/N
    x = a

    for i in range(N-1):
        yi = y[:,i]
        k = h*f(x,yi)
        y[:,i+1] = yi + h*f(x+h/2,yi+k/2)
        x += h
    return y

def rk3(f,a,b,y0,N):
    y = np.zeros((len(y0),N))
    y[:,0] = y0
    h = (b-a)/N
    x = a

    for i in range(N-1):
        yi = y[:,i]
        k1 = h*f(x,yi)
        k2 = h*f(x+h/2,yi+k1/2)
        k3 = h*f(x+h,yi-k1+2*k2)
        y[:,i+1] = yi + (1/6)*(k1+4*k2+k3)
        x += h
    return y


def rk4(f,a,b,y0,N):
    y = np.zeros((len(y0),N))
    y[:,0] = y0
    h = (b-a)/N
    x = a

    for i in range(N-1):
        yi = y[:,i]
        k1 = h*f(x,yi)
        k2 = h*f(x+h/2,yi+k1/2)
        k3 = h*f(x+h/2,yi+k2/2)
        k4 =h*f(x+h,yi+k3)
        y[:,i+1] = yi + (1/6)*(k1+2*k2+2*k3+k4)
        x += h
    return y

def f(x,y):
    result = np.ones((2))
    result[0] = y[1]
    result[1] = -4*np.pi*np.pi*y[0]
    return result

a = 0
b = 40
y0 = np.array([1,0])
N = 10000
yexact = lambda x:np.cos(2*np.pi*x)
pexact = lambda x:-np.sin(2*np.pi*x)*2*np.pi
""" solrk2 = rk2(f,a,b,y0,N)
solrk3 = rk3(f,a,b,y0,N) """
""" solrk4 = rk4(f,a,b,y0,N)
x = np.linspace(a,b,N)
soly = solrk4[0,:]
solp = solrk4[1,:]
erry = abs(soly-yexact(x))
errp = abs(solp-pexact(x))
plt.subplot(2,1,1)
plt.plot(x,solrk2[0,:],label="RK2")
plt.plot(x,solrk3[0,:],label="RK3")
plt.plot(x[-1000:],soly[-1000:],label="RK4")
plt.plot(x[-1000:],yexact(x[-1000:]))
plt.xlabel("t")
plt.ylabel("y")
plt.title("Solutions for t in range 0-40 of y")
plt.subplot(2,1,2)
plt.plot(x,solrk2[1,:],label="RK2")
plt.plot(x,solrk3[1,:],label="RK3")
plt.plot(x[-1000:],solp[-1000:],label="RK4")
plt.plot(x[-1000:],pexact(x[-1000:]))
plt.xlabel("t")
plt.ylabel("p")
plt.title("Solutions at t in range 0-40 of p")
plt.suptitle("y and p calculated with Runge-Kutta 4 in the timerange 0-40")
plt.tight_layout()
plt.show()
plt.plot(x,erry)
plt.plot(x,errp)
plt.show()
 """

N = np.array([500*2**i for i in range(4)])
h = (b-a)/N
tLast = -100
allsoly = []
allsolp = []
totery = []
toterp = []
xs = []
for i in range(4):
    sol = rk4(f,a,b,y0,N[i])
    soly = sol[0,:]
    solp = sol[1,:]
    x = np.linspace(a,b,N[i])
    erry = abs(soly-yexact(x))
    errp = abs(solp-pexact(x))
    xs.append(x)
    allsoly.append(soly)
    allsolp.append(solp)
    totery.append(erry)
    toterp.append(errp)

plt.figure(1)
for i in range(len(N)):
    plt.plot(xs[i],totery[i],label=f"h = {h[i]}")
plt.legend()
plt.xlabel("t")
plt.ylabel("|y-y_exact|")
plt.title("Error of y for different h")
plt.savefig("errory.png")

plt.figure(2)
for i in range(len(N)):
    plt.plot(xs[i],toterp[i],label=f"h = {h[i]}")
plt.legend()
plt.xlabel("t")
plt.ylabel("|p-p_exact|")
plt.title("Error of p for different h")
plt.savefig("errorp.png")

plt.figure(3)
plt.plot(xs[-1],yexact(xs[-1]),label="y (Exact)")
plt.plot(xs[-1],allsoly[-1],label =f"y (RK4 with h = {h[-1]})")
plt.legend()
plt.xlabel("t")
plt.ylabel("y")
plt.title("Exact and approximated solution of y")
plt.savefig("ysol.png")

plt.figure(4)
plt.plot(xs[-1],pexact(xs[-1]),label="p (Exact)")
plt.plot(xs[-1],allsolp[-1],label =f"p (RK4 with h = {h[-1]})")
plt.legend()
plt.xlabel("t")
plt.ylabel("p")
plt.title("Exact and approximated solution of p")
plt.savefig("psol.png")

plt.figure(5)
plt.plot(xs[-1][tLast:],yexact(xs[-1][tLast:]),label="y (Exact)")
plt.plot(xs[-1][tLast:],allsoly[-1][tLast:],label =f"y (RK4 with h = {h[-1]})")
plt.legend()
plt.xlabel("t")
plt.ylabel("y")
plt.title("Exact and approximated solution of y")
plt.savefig("ysolzoom.png")

plt.figure(6)
plt.plot(xs[-1][tLast:],pexact(xs[-1][tLast:]),label="p (Exact)")
plt.plot(xs[-1][tLast:],allsolp[-1][tLast:],label =f"p (RK4 with h = {h[-1]})")
plt.legend()
plt.xlabel("t")
plt.ylabel("p")
plt.title("Exact and approximated solution of p")
plt.savefig("psolzoom.png")