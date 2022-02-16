import random

proc step(pos: (int, int)): (int,int) = 
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
    var path: seq[(int,int)]
    var start = (0,0)
    path.add(start)
    for i in 1..100:
        start = step(start)
        path.add(start)
    echo path

main()