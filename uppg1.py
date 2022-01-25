import math

def sym3(f,x,h):
    return (f(x+h)-f(x-h))/(2*h)

def forward(f,x,h):
    return (f(x+h)-f(x))/h

def backward(f,x,h):
    return (f(x)-f(x-h))/h

def sym5(f,x,h):
    (1/(12*h))*(f(x-2*h)-8*f(x-h)+8*f(x+h)-f(x+2*h))

f = math.sin
h = [5e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5]
x = 1

ressym3 = []
resforward = []
resbackward = []
ressym5 = []
exact = math.cos(1)

for i in h:
    ressym3.append(abs(sym3(f,x,i)-exact))
    resforward.append(abs(forward(f,x,i)-exact))
    resbackward.append(abs(backward(f,x,i)-exact))
    ressym5.append(abs(sym5(f,x,i)-exact))

print(ressym3)
print(resforward)
print(resbackward)
print(ressym5)