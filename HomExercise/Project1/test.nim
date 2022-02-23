import utils
import math

proc exactSolution*(b, r_max, u0, E: float): float =
  let r_min = b / sqrt(1 - u0 / E) # calculate this by finding root of 1 - b^2 / r^2 - V / E
  if E < u0:
    result = arctan(b / r_max / sqrt(1 - b*b / (r_max*r_max) - u0 / E)) - arctan(b / r_min / sqrt(1 - b*b / (r_min*r_min) - u0 / E)) + arcsin(1.0) - arcsin(b / r_max)
    result *= 2



proc thetafunc(u0,r_max,E,b: float) =

    let r_min = b / sqrt(1 - u0 / E)

    proc f1(x: float): float =
        return 1/sqrt(x-b)*1/sqrt(1-abs(x-b)-u0/E)

    proc f2(x: float): float =
        return 1/sqrt(x-b)*1/sqrt(1-abs(x-b))

    let a1 = b^2/r_min^2 + b
    let a2 = 1+b
    let b1 = b^2/r_max^2 + b
    let b2 = b^2/r_max^2 + b
    let int1 = gaussQuad(f1,a1,b1,5)
    let int2 = gaussQuad(f2,a2,b2,5)
    let exact = exactSolution(b,r_max,u0,E)
    let res = int1+int2
    echo res
    echo exact
    echo r_min
thetafunc(-1,1,-2,0.1)