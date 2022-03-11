import ../Project1/utils
import math, strformat
import numericalnim
import ggplotnim




proc f1(b:float,p:float):float =
  return 2*b*(2*p/(p^2+b)^2*(1/sqrt(1-(b^2)/(p^2+b)^2)))

proc f2(b:float,p:float,u:proc(r: float): float,e:float,rmin:float):float =
  return -2*b*(2*p/(p^2+rmin)^2*(1/sqrt(1-(b^2)/(p^2+rmin)^2-u(p*p + rmin)/e)))

proc main() =
    let a = 1.0
    let rmax = 3 * a
    let v0 = 1.0
    let b = 0.5
    let E = linspace(0.1*v0, 100*v0,100)
    let min = 0.0001
    let U0 = proc (r: float): float = 4*v0*((a / r)^12 - (a / r)^6)
    let l = 20
    var ftot:seq[float]
    var ftot2:seq[float]
    for e in E:
      let e = e
      proc fhelp(x:float):float =
        f1(b,x)
      
      proc rminRoot(r: float): float =
        1 - b*b / (r*r) - U0(r) / e

      let rmin = bisection(rminRoot, 0.75, 2.25, 1e-8)
      echo "rmin: ", rmin

      proc fhelp2(x:float):float =
        f2(b,x,U0,e,rmin)
      ftot.add(gaussQuad(fhelp,min,sqrt(rmax-b),l))
      ftot2.add(gaussQuad(fhelp,min,sqrt(rmax-b),l)+gaussQuad(fhelp2,min,sqrt(rmax-rmin),l))
    echo ftot


    #[ var df1 = seqsToDf({"b":E, "Numerical":ftot})
    df1 = df1.gather(["Numerical"], key="type", value="f value")
    echo df1
    ggplot(df1,aes("b","f value", color = "type")) +
      geom_line() +
      ggtitle("Plot over f") +
      ggsave("fplot.png")
 ]#
    var df2 = seqsToDf({"b":E, "Numerical":ftot2})
    df2 = df2.gather(["Numerical"], key="type", value="f value")
    ggplot(df2,aes("b","f value", color = "type")) +
      geom_line() +
      ggtitle("Plot over f") +
      ggsave("fplot2.png")

# xs = b / sin(θ) * |(dθ/db)^-1|

proc crossSection(allTheta: seq[seq[float]], B: seq[float]): seq[seq[float]] =
  for thetaB in allTheta:
    var resTheta: seq[float]
    for i in 0 .. thetaB.high:
      let deriv =
        if i == 0:
          (thetaB[i+1] - thetaB[i]) / (B[i+1] - B[i])
        elif i == thetaB.high:
          (thetaB[i] - thetaB[i-1]) / (B[i] - B[i-1])
        else:
          (thetaB[i+1] - thetaB[i-1]) / (B[i+1] - B[i-1])
      resTheta.add B[i] / sin(thetaB[i]) * abs(1 / deriv)
    result.add resTheta

proc plotData(y: seq[seq[float]], x: seq[float], E: seq[float], filename, title, ylabel: string) =
  var df = seqsToDf({"b": x, fmt"E = {E[0]:6.3f}": y[0]})
  var columnNames: seq[string] = @[fmt"E = {E[0]:6.3f}"]
  for i in 1 .. y.high:
    df = innerJoin(df, seqsToDf({"b": x, fmt"E = {E[i]:6.3f}": y[i]}), by = "b")
    columnNames.add fmt"E = {E[i]:6.3f}"

  df = df.gather(columnNames, key="Energy", value=yLabel)

  ggplot(df,aes("b",yLabel, color = "Energy")) +
    geom_line() +
    ggtitle(title) +
    ggsave(filename)

proc main2() =
  let a = 1.0
  let rmax = 3 * a
  let v0 = 1.0
  let B = linspace(0.1, rmax, 100)
  let E = geomspace(0.1*v0, 100*v0,10)
  let min = 0.0001
  let U0 = proc (r: float): float = 4*v0*((a / r)^12 - (a / r)^6)
  let l = 20
  var allTheta: seq[seq[float]]
  for e in E:
    let e = e
    var theta: seq[float]
    for b in B:
      let b = b
      proc fhelp(x:float):float =
        f1(b,x)
      
      proc rminRoot(r: float): float =
        1 - b*b / (r*r) - U0(r) / e

      let rmin = max(rootsFinder(rminRoot, 0.75, 3.25, 1e-8))
      echo "rmin: ", rmin

      proc fhelp2(x:float):float =
        f2(b,x,U0,e,rmin)
      theta.add(gaussQuad(fhelp,min,sqrt(rmax-b),l)+gaussQuad(fhelp2,min,sqrt(rmax-rmin),l))
    allTheta.add theta

  let xs = crossSection(allTheta, B)
  plotData(allTheta, B, E, "thetaPlot.png", "θ as a function of b for different energies", "θ")
  plotData(xs, B, E, "xsPlot.png", "Cross section as a function of b for different energies", "θ")
main()
main2()