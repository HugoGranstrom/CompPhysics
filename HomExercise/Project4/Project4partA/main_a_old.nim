import math, random, sequtils, ggplotnim, numericalnim
import locks, taskpools, cpuinfo

type
  Lattice* = object
    width*, height*: int
    data*: seq[float]
    J*, B*, T*: float 

const kb = 8.617e-5
const J = -0.000189574

proc `[]`*(lattice: Lattice, row, col: int): float =
  # use modulus to wrap around periodic BC. 
  let modRow = (row + lattice.height - 1) mod (lattice.height - 1)
  let modCol = (col + lattice.width - 1) mod (lattice.width - 1)
  #echo modRow, " ", modCol
  result = lattice.data[(lattice.width-1)*modRow + modCol]

proc `[]=`*(lattice: var Lattice, row, col: int, val: float) =
  let modRow = (row + lattice.height - 1) mod (lattice.height - 1)
  let modCol = (col + lattice.width - 1) mod (lattice.width - 1)
  lattice.data[(lattice.width-1)*modRow + modCol] = val

proc flip*(lattice: var Lattice, row, col: int) =
  lattice[row, col] = lattice[row, col] * -1.0

proc newLattice*(height, width: int, J, B, T: float): Lattice =
  result.height = height
  result.width = width
  result.J = J
  result.B = B
  result.T = T
  result.data = newSeqWith[float]((height-1)*(width-1), sample([-1.0, 1.0]))


proc calcHamiltonian(lattice: Lattice, row, col: int): float =
  let current = lattice[row, col]
  result += current * lattice[row-1, col] #
  result += current * lattice[row+1, col]
  result += current * lattice[row, col-1] #
  result += current * lattice[row, col+1]

proc calcHamiltonian*(lattice: Lattice): float =
  result = -lattice.B * sum(lattice.data)
  var pairsSum: float
  for row in 0 ..< lattice.height - 1:
    for col in 0 ..< lattice.width - 1:
      pairsSum += calcHamiltonian(lattice, row, col)
  result += -lattice.J * pairsSum

proc M_calc(lattice: Lattice): float =
  return sum(lattice.data)/lattice.data.len.toFloat

proc Heat_calc(E: float, eSquare: float, T: float): float =
  return (eSquare - E*E)/(kb*T*T)

proc Sus_calc(M: float, mSquare: float, T:float): float =
  return (mSquare-M*M)/(kb*T)

proc Cumulant(mSquare: float, m4: float): float =
  return 1 - m4/(3*mSquare*mSquare)

proc mean(x: seq[float]): float =
  return sum(x)/x.len.float


proc normSpace(xmin: float, xmax:float, xsum:int, xsize:float, xmid:float): seq[float] =
  
  
  var c1 = 0.0
  var c2 = 0.0
  
  let temp = linspace(xmin,xmax,10000)
  for i in temp:
    if i in xmin..xmid-(xsize/2):
      c1 += 1
    elif i in xmid+(xsize/2)..xmax:
      c2 += 1
  c1 /= 10000.0
  c2 /= 10000.0

  c1 = c1/(c1+c2)
  c2 = 1-c1
  var coarse1:seq[float]
  var coarse2:seq[float]
  var fine:seq[float]
  if int(c1*(xsum/2)) != 0:
    coarse1 = linspace(xmin,(xmid-(xsize/2)),int(c1*(xsum/2)))
  if int(c2*(xsum/2)) != 0:
    coarse2 = linspace((xmid-(xsize/2)), xmid+(xsize/2),int(c2*(xsum/2)))
  fine = linspace((xmid-(xsize/2)), xmid+(xsize/2), int(xsum/2))
  return concat(coarse1,fine,coarse2)




proc ising*(lattice: var Lattice): (float, float, float, float, float) =
  var hamiltonian = calcHamiltonian(lattice)
  var dE = 1e10
  var iters: int
  var eTot: seq[float]
  var mTot: seq[float]
  var msquareTot: seq[float]
  var esquareTot: seq[float]
  var m4Tot: seq[float]
  while abs(dE) > 1e-10 or iters < 100:
    iters += 1
    var innerHamiltonian = hamiltonian
    for row in 0 ..< lattice.height - 1:
      for col in 0 ..< lattice.width - 1:
        lattice.flip(row, col)
        let newHamiltonian = calcHamiltonian(lattice)
        let dEinner = newHamiltonian - innerhamiltonian
        if dEinner > 0:
          let r = rand(1.0)
          #if iters < 10: echo -dEinner / (lattice.T * kb)
          if r > exp(-dEinner / (lattice.T * kb)):
            #echo "Rejected! ", dEinner
            lattice.flip(row, col) # flip back
          else:
            innerHamiltonian = newHamiltonian
        else:
          innerhamiltonian = newHamiltonian
    if iters > 50:
      let M = M_calc(lattice)
      let msquare = M*M
      let m4 = msquare*msquare
      let esquare = innerHamiltonian*innerHamiltonian
      eTot.add innerHamiltonian
      mTot.add M
      msquareTot.add msquare
      esquareTot.add esquare
      m4Tot.add m4
    
    dE = lattice.calcHamiltonian() - hamiltonian
    #echo dE
    hamiltonian = innerHamiltonian
  let M = mean(mTot)
  let msquare = mean(msquareTot)
  let m4 = mean(m4Tot)
  let e = mean(eTot)
  let esquare = mean(esquareTot)
  #echo "Iterations: ", iters
  return (M,msquare,m4,e,esquare)
  

var modLock: Lock

proc thread_func(task_id: int, c_len: int, cs: ptr UncheckedArray[float], avgHeat: ptr UncheckedArray[float], avgSus: ptr UncheckedArray[float], avgCumul: ptr UncheckedArray[float], avgM: ptr UncheckedArray[float],latticeZise: int): bool = 
  randomize(task_id)
  echo "THREADFUNC!!!!"
  for i in 0 ..< c_len:
    echo "ITERATION: ",i
    let c = cs[i]

    let T = J*c/kb
    #let c = kb*T/J
      
    #echo "c = ", c
    let B = 0.001
    
    var lattice = newLattice(latticeZise, latticeZise, J, B, T)
    #echo lattice.calcHamiltonian
    echo "Innan Ising"
    let (M,msquare,m4,e,esquare) = ising(lattice)
    echo "Efter Isiag"
    #echo lattice.calcHamiltonian
    let specHeat = Heat_calc(e,esquare,T)
    let sus = Sus_calc(M,msquare,T)
    let cumul = Cumulant(msquare,m4)

    #echo "Specific heat = ", specHeat
    #echo "Suseptibility = ", sus
    #echo "Cumulant = ", cumul
    #echo "Magnetization = ", M
    echo "INNAN LOCK"
    withLock(modLock):
      avgHeat[i] += specHeat
      avgSus[i] += sus
      avgCumul[i] += abs(cumul)
      avgM[i] += abs(M)
    echo "EFTER LOCK"


proc main() =

#SIZE 8
  echo "Running for size 8"
  let c_len = 20
  var avgheat1: seq[float] = newSeq[float](c_len)
  var avgsus1: seq[float] = newSeq[float](c_len)
  var avgcumul1: seq[float] = newSeq[float](c_len)
  var avgM1: seq[float] = newSeq[float](c_len)
  #randomize()
  var cs = linspace(-10,-0.1,10)
  let times_run = 20
  var latticeZise = 8
  var nthreads = countProcessors()
  var tp = Taskpool.new(num_threads = nthreads)
  var pendingTasks = newSeq[FlowVar[bool]](times_run)
  modLock.initLock
  for n in 0 ..< times_run:
    pendingTasks[n] = tp.spawn thread_func(n, c_len, cast[ptr UncheckedArray[float]](cs[0].addr), cast[ptr UncheckedArray[float]](avgheat1[0].addr), cast[ptr UncheckedArray[float]](avgsus1[0].addr), cast[ptr UncheckedArray[float]](avgcumul1[0].addr), cast[ptr UncheckedArray[float]](avgM1[0].addr),latticeZise)
  
  for n in 0 ..< times_run:
    discard sync pendingTasks[n]

  for i in 0 .. cs.high:
    avgheat1[i] /= float(times_run)
    avgsus1[i] /= float(times_run)
    avgcumul1[i] /= float(times_run)
    avgM1[i] /= float(times_run)

# SIZE 16
  echo "Running for size 16"
  var avgheat2: seq[float] = newSeq[float](c_len)
  var avgsus2: seq[float] = newSeq[float](c_len)
  var avgcumul2: seq[float] = newSeq[float](c_len)
  var avgM2: seq[float] = newSeq[float](c_len)
  latticeZise = 16
  modLock.initLock
  for n in 0 ..< times_run:
    pendingTasks[n] = tp.spawn thread_func(n, c_len, cast[ptr UncheckedArray[float]](cs[0].addr), cast[ptr UncheckedArray[float]](avgheat2[0].addr), cast[ptr UncheckedArray[float]](avgsus2[0].addr), cast[ptr UncheckedArray[float]](avgcumul2[0].addr), cast[ptr UncheckedArray[float]](avgM2[0].addr),latticeZise)
  
  for n in 0 ..< times_run:
    discard sync pendingTasks[n]

  for i in 0 .. cs.high:
    avgheat2[i] /= float(times_run)
    avgsus2[i] /= float(times_run)
    avgcumul2[i] /= float(times_run)
    avgM2[i] /= float(times_run)

# SIZE 32

  echo "Running for size 32"
  var avgheat3: seq[float] = newSeq[float](c_len)
  var avgsus3: seq[float] = newSeq[float](c_len)
  var avgcumul3: seq[float] = newSeq[float](c_len)
  var avgM3: seq[float] = newSeq[float](c_len)
  latticeZise = 32
  modLock.initLock
  for n in 0 ..< times_run:
    pendingTasks[n] = tp.spawn thread_func(n, c_len, cast[ptr UncheckedArray[float]](cs[0].addr), cast[ptr UncheckedArray[float]](avgheat3[0].addr), cast[ptr UncheckedArray[float]](avgsus3[0].addr), cast[ptr UncheckedArray[float]](avgcumul3[0].addr), cast[ptr UncheckedArray[float]](avgM3[0].addr),latticeZise)
  
  for n in 0 ..< times_run:
    discard sync pendingTasks[n]
    echo "Progress: ", 100 * (n+1) / times_run, "%" 

  for i in 0 .. cs.high:
    avgheat3[i] /= float(times_run)
    avgsus3[i] /= float(times_run)
    avgcumul3[i] /= float(times_run)
    avgM3[i] /= float(times_run)
    
  echo "Plotting"
  var dfall = seqsToDf({"c":cs,"L=8M1":avgM1, "L=16M1":avgM2, "L=32M1":avgM3,"L=8S1":avgsus1, "L=16S2":avgsus2,"L=32S3":avgsus3,
 "L=8Cu1":avgcumul1,"L=16Cu2":avgcumul2,"L=32Cu3":avgcumul3, "L=8H1":avgheat1,"L=16H2":avgheat2,"L=32H3":avgheat3})
  var df1 = seqsToDf({"c":cs,"L=8":avgM1, "L=16":avgM2, "L=32":avgM3})
  var df2 = seqsToDf({"c":cs,"L=8":avgsus1, "L=16":avgsus2,"L=32":avgsus3})
  var df3 = seqsToDf({"c":cs, "L=8":avgcumul1,"L=16":avgcumul2,"L=32":avgcumul3})
  var df4 = seqsToDf({"c":cs,"L=8":avgheat1,"L=16":avgheat2,"L=32":avgheat3})
  df1 = df1.gather(["L=8","L=16","L=32"], key = "size", value="Magnitization")
  df2 = df2.gather(["L=8","L=16","L=32"], key = "size", value="Suseptibility")
  df3 = df3.gather(["L=8","L=16","L=32"], key = "size", value="Cumulant")
  df4 = df4.gather(["L=8","L=16","L=32"], key = "size", value="Specific Heat")
  ggplot(df1, aes("c","Magnitization",color = "size")) +
    geom_point() +
    geom_line() +
    ggtitle("Magnitization") +
    ggsave("mag.png")
  ggplot(df2, aes("c","Suseptibility",color = "size")) +
    geom_point() +
    geom_line() +
    ggtitle("Susceptability") +
    ggsave("sus.png")
  ggplot(df3, aes("c","Cumulant",color = "size")) +
    geom_point() + 
    geom_line() +
    ggtitle("Cumulant") +
    ggsave("cumul.png")
  ggplot(df4, aes("c","Specific Heat",color = "size")) +
    geom_point() +
    geom_line() +
    ggtitle("Specific Heat") +
    ggsave("heat.png")
  dfall.writeCsv("test.csv")

main()