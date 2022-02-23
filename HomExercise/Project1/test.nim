import utils
import math

proc f(x: float): float =
    return x^2

echo gaussQuad(f, 0, 3, 3)