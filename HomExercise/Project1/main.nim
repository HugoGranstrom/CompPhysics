import utils
import math
import numericalnim
import ggplotnim

proc f1(b:float,p:float):float =
  return 2*b*(2*p/(p^2+b)^2*(1/sqrt(1-(b^2)/(p^2+b)^2)))

proc f2(b:float,p:float,u:float,e:float,rmin:float):float =
  return -2*b*(2*p/(p^2+rmin)^2*(1/sqrt(1-(b^2)/(p^2+rmin)^2-u/e)))

proc exact(b:float,rmax:float):float =
  return PI-2*arcsin(b/rmax)

proc exact2(b:float, rmax:float, u:float, e:float):float =
  return 2*arcsin(b/(rmax*sqrt(1-u/e)))-2*arcsin(b/rmax)

proc main() =
    let rmax = 100.0
    let min = 0.01
    let b = linspace(min,rmax,100)
    let E = 1.0
    let U0 = -0.1
    let l = 20
    var ftot:seq[float]
    var antot:seq[float]
    var ftot2:seq[float]
    var antot2:seq[float]
    for i in b:
      let i = i
      proc fhelp(x:float):float =
        f1(i,x)
      
      var rmin = sqrt(i^2/(1-U0/E))

      proc fhelp2(x:float):float =
        f2(i,x,U0,E,rmin)
      ftot.add(gaussQuad(fhelp,min,sqrt(rmax-i),l))
      antot.add(exact(i,rmax))
      ftot2.add(gaussQuad(fhelp,min,sqrt(rmax-i),l)+gaussQuad(fhelp2,min,sqrt(rmax-rmin),l))
      antot2.add(exact2(i,rmax,U0,E))
    echo ftot
    echo antot


    var df1 = seqsToDf({"b":b, "Numerical":ftot, "Analytical":antot})
    df1 = df1.gather(["Numerical", "Analytical"], key="type", value="f value")
    echo df1
    ggplot(df1,aes("b","f value", color = "type")) +
      geom_line() +
      ggtitle("Plot over f") +
      ggsave("fplot.png")

    var df2 = seqsToDf({"b":b, "Numerical":ftot2, "Analytical":antot2})
    df2 = df2.gather(["Numerical", "Analytical"], key="type", value="f value")
    ggplot(df2,aes("b","f value", color = "type")) +
      geom_line() +
      ggtitle("Plot over f") +
      ggsave("fplot2.png")

main()
