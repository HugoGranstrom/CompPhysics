import math
import random as rd
import matplotlib.pyplot as plt

def piapprox(n):
    totpoint = [0]
    circpoint = [0]
    for i in range(n):
        x = rd.uniform(-1,1)
        y = rd.uniform(-1,1)
        if x**2+y**2 < 1:
            circpoint.append(circpoint[-1]+1)
        else:
            circpoint.append(circpoint[-1])     
        totpoint.append(totpoint[-1]+1)
    piapprox = [4*circpoint[i]/totpoint[i] for i in range(1,len(totpoint))]
    err = [abs(piapprox[i] - math.pi) for i in range(len(piapprox))]

    plt.plot(piapprox)
    plt.show()
    plt.semilogy(err)
    plt.show()
    print(piapprox[-1])
    print(math.pi)

piapprox(10000000)


