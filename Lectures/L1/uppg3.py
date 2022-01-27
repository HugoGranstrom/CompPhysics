import matplotlib.pyplot as plt

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

def newton(f,derf,x0,tol):
    x = x0
    xlist = [x]
    ylist = [f(x)]
    diff = tol+1
    while diff > tol:
        x2 = x-(f(x)/derf(x))
        diff = abs(x2-x)
        x = x2
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


f = lambda x:x**2-5
derf = lambda x:2*x
x0 = 0.1
x1 = 1
tol = 1e-6
stepsize = 1
searchx,searchy = search(f,x0,tol,stepsize)
newx,newy = newton(f,derf,x0,tol)
secx,secy = secant(f,x0,x1,tol)
exact = 5**0.5
print(searchx)
print(searchy)
print(newx)
print(newy)
print(secx)
print(secy)
errsearch = [abs(i-exact) for i in searchx]
errnew = [abs(i-exact) for i in newx]
errsec = [abs(i-exact) for i in secx]
plt.semilogy(errsearch)
plt.semilogy(errnew)
plt.semilogy(errsec)
plt.legend(["Search","Newton","Secant"])
plt.show()