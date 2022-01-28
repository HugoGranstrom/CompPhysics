import matplotlib.pyplot as plt
import numpy as np

def rk2(f,a,b,y0,N):
    y = np.zeros(len(y0),N)
    y[:,0] = y0
    h = (b-a)/N
    x = a

    for i in range(N):
        yi = y[:,i]
        k = h*f(x,yi)
        y[:,i+1] = yi + h*f(x+h/2,yi+k/2)
        x += h
    return y

def rk3(f,a,b,y0,N):
    y = np.zeros(len(y0),N)
    y[:,0] = y0
    h = (b-a)/N
    x = a

    for i in range(N):
        yi = y[:,i]
        k1 = h*f(x,yi)
        k2 = h*f(x+h/2,yi+k1/2)
        k3 = h*f(x+h,yi-k1+2*k2)
        y[:,i+1] = yi + (1/6)*(k1+4*k2+k3)
        x += h
    return y


def rk4(f,a,b,y0,N):
    y = np.zeros(len(y0),N)
    y[:,0] = y0
    h = (b-a)/N
    x = a

    for i in range(N):
        yi = y[:,i]
        k1 = h*f(x,yi)
        k2 = h*f(x+h/2,yi+k1/2)
        k3 = h*f(x+h/2,yi+k2/2)
        k4 =h*f(x+h,yi+k3)
        y[:,i+1] = yi + (1/6)*(k1+2*k2+2*k3+k4)
        x += h
    return y

def f(x,y):
    result = np.ones(2,1)
    result[0] = y[1]
    result[1] = -4*np.pi*np.pi*y[0]
    return result