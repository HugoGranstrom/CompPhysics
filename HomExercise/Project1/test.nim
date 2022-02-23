import utils
import math

proc exactSolution*(theta: float, b, r_max, u0, E: float): float =
  let r_min = b / sqrt(1 - u0 / E) # calculate this by finding root of 1 - b^2 / r^2 - V / E
  if E < u0:
    result = arctan(b / r_max / sqrt(1 - b*b / (r_max*r_max) - u0 / E)) - arctan(b / r_min / sqrt(1 - b*b / (r_min*r_min) - u0 / E)) + arcsin(1.0) - arcsin(b / r_max)
    result *= 2

proc f(x: float): float =
    return x^2

echo gaussQuad(f, 0, 3, 3)
