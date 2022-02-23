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
            return fc
        else:
            echo "Something went wrong in bisection!"
    return (a+b)/2

proc rootsFinder*(f: proc (x: float): float,a: float,b: float, tol: float): seq[float] =
    var intervals: seq[(float, float)]
    var currentValue = f(a)
    const step = 0.01
    var x = a
    while x < b:
        var temp = f(x)
        if temp*currentValue < 0:
            intervals.add((x-step,x))
        currentValue = temp
        x += step
    var roots: seq[float]
    for interval in intervals:
        roots.add(bisection(f, interval[0],interval[1],tol))
    return roots

#[ def legendre(x, l):
  p = [0*x + 1, x]
  for i in range(1, l):
    p.append(((2*i + 1)*x*p[i] - i*p[i-1]) / (i + 1))
  return p[l]

def legendreDeriv(x, l):
  if l == 0: return x*0
  return (-l*x*legendre(x, l) + l*legendre(x, l-1)) / (1 - x*x)

def bisection(f, x0, x1, tol):
  a = x0
  b = x1
  fa = f(a)
  fb = f(b)
  assert fa * fb < 0
  while abs(b - a) > tol:
    c = (a + b) / 2
    fc = f(c)
    if fc*fa < 0:
      b = c
      fb = fc
    elif fc*fb < 0:
      a = c
      fa = fc
    elif fc == 0:
      return fc
    else:
      print("Something went wrong in bisection!")
  return (a + b) / 2

def rootsFinder(f, a, b, tol):
  intervals = []
  currentValue = f(a)
  step = 0.01
  for x in np.arange(a, b, step):
    temp = f(x)
    if temp*currentValue < 0:
      intervals.append((x-step, x))
    currentValue = temp
  roots = [bisection(f, interval[0], interval[1], tol) for interval in intervals]
  return roots ]#