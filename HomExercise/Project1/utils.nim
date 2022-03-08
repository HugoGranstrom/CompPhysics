# Gaussian quadrature - Legendre Polynomials
# Root finder
import std/[sequtils, math]
#import numericalnim, ggplotnim


proc legendre*(x: float, l: int): float =
  var p = @[1.0, x]
  for i in 1 ..< l:
    let iFloat = i.toFloat
    p.add ((2*iFloat + 1)*x*p[i] - iFloat*p[i-1]) / (iFloat + 1)
  return p[l]

proc legendreDeriv*(x: float, l: int): float =
  if l == 0: return 0.0
  let lFloat = l.toFloat
  return (-lFloat*x*legendre(x, l) + lFloat*legendre(x, l-1)) / (1 - x*x)

#[ let x = linspace(-1, 1, 100)
let y = x.mapIt(legendre(it, 5))

let df = seqsToDf(x, y)
ggplot(df, aes("x", "y")) +
  geom_line() +
  xlim(-1, 1) +
  ggsave("legendre.png") ]#

proc bisection*(f: proc (x: float): float,x0: float,x1: float, tol: float): float =
    var a = x0
    var b = x1
    var c: float
    var fa = f(a)
    var fb = f(b)
    var fc: float
    assert fa*fb < 0
    while abs(b-a) > tol:
        c = (a+b)/2
        fc = f(c)
        if fc*fa < 0:
            b = c
            fb = fc
        elif fc*fb < 0:
            a = c
            fa = fc
        elif fc == 0:
            return c
        else:
            echo "Something went wrong in bisection!"
    return (a+b)/2

proc rootsFinder*(f: proc (x: float): float,a: float,b: float, tol: float): seq[float] =
    var intervals: seq[(float, float)]
    var currentValue = f(a)
    const step = 0.01
    var x = a
    while x < b:
      #echo x
      var temp = f(x)
      if temp*currentValue < 0:
        #echo "hah"
        intervals.add((x-step,x))
      currentValue = temp
      x += step
    var roots: seq[float]
    for interval in intervals:
        roots.add(bisection(f, interval[0],interval[1],tol))
    return roots

proc gaussQuad*(f: proc (x: float): float, a, b: float, l: int): float =
  let abscissae = rootsFinder(proc (x: float): float = legendre(x, l), -1, 1, 1e-15) 
  #echo abscissae
  let weights = abscissae.mapIt(2 / ((1 - it*it) * legendreDeriv(it, l)^2))
  #echo "weights: ", weights
  let new_f =
    proc (t: float): float =
      return (b - a) / 2 * f((t+1)*(b-a)/2 + a)
  for (x, w) in zip(abscissae, weights):
    result += w * new_f(x)

