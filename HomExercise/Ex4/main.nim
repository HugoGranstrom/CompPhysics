import random
import math
import ggplotnim
#randomize(1)
proc step(pos: (float, float)): (float,float) = 
    let num = rand(3)
    if num == 0:
        return (pos[0]+1,pos[1])
    elif num == 1:
        return (pos[0]-1,pos[1])
    elif num == 2:
        return (pos[0],pos[1]+1)
    elif num == 3:
        return (pos[0],pos[1]-1)

proc main() =
    const steps = 10000
    const pers = 100000
    var path: seq[(float,float)] = newSeq[(float, float)](pers)
    var xList: seq[float]
    var yList: seq[float]
    var rList: seq[float]
    var nList: seq[float]

    xList.add 0.0
    yList.add 0.0
    for i in 0 ..< steps:
        var x = 0.0
        var y = 0.0
        var r2 = 0.0
        for j in 0 ..< pers:
            path[j] = step(path[j])
            x += path[j][0]
            y += path[j][1]
            r2 += x^2 + y^2 
        xList.add x / pers
        yList.add y / pers
        if i > 10:
            rList.add sqrt(r2 / pers)
            nList.add float(i)
    let df1 = seqsToDf(xList, yList)
    ggplot(df1, aes("xList","yList")) + geom_line() + ggsave("xy.png")
    var df2 = seqsToDf({"r": rList, "N": nList})
    df2 = df2.mutate(f{"sqrt(N)" ~ sqrt(`N`)*(400/3)})
    df2 = df2.gather(["r", "sqrt(N)"], key = "type", value = "R")
    ggplot(df2, aes("N","R",color = "type")) + 
        geom_line() + 
        ggtitle("loglog plot of r and sqrt(N) plotted against N") +
        scale_x_log10() +
        scale_y_log10() +
        ylab("sqrt(<r²>)") +
        ggsave("rs.png")
main()