import math, random, sequtils, ggplotnim, numericalnim
import locks, taskpools, cpuinfo

type
  Lattice* = object
    width*, height*: int
    data*: seq[float]
    J*, B*, T*: float 

const kb = 8.617e-5
const J = 0.000189574

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

proc thread_func(task_id: int, c_len: int, cs: ptr UncheckedArray[float], avgHeat: ptr UncheckedArray[float], avgSus: ptr UncheckedArray[float], avgCumul: ptr UncheckedArray[float], avgM: ptr UncheckedArray[float]): bool = 
  randomize(task_id)
  for i in 0 ..< c_len:
    let c = cs[i]
    let T = J*c/kb
    #let c = kb*T/J
      
    #echo "c = ", c
    let B = 0.0
    
    var lattice = newLattice(10, 10, J, B, T)
    #echo lattice.calcHamiltonian
    let (M,msquare,m4,e,esquare) = ising(lattice)
    #echo lattice.calcHamiltonian
    let specHeat = Heat_calc(e,esquare,T)
    let sus = Sus_calc(M,msquare,T)
    let cumul = Cumulant(msquare,m4)

    #echo "Specific heat = ", specHeat
    #echo "Suseptibility = ", sus
    #echo "Cumulant = ", cumul
    #echo "Magnetization = ", M
    withLock(modLock):
      avgHeat[i] += specHeat
      avgSus[i] += sus
      avgCumul[i] += abs(cumul)
      avgM[i] += abs(M)


proc main() =
  let c_len = 500
  var avgheat: seq[float] = newSeq[float](c_len)
  var avgsus: seq[float] = newSeq[float](c_len)
  var avgcumul: seq[float] = newSeq[float](c_len)
  var avgM: seq[float] = newSeq[float](c_len)
  var cs = linspace(2.0,7.0,c_len)
  let times_run = 100

  var nthreads = countProcessors()
  var tp = Taskpool.new(num_threads = nthreads)
  var pendingTasks = newSeq[FlowVar[bool]](times_run)
  modLock.initLock
  for n in 0 ..< times_run:
    pendingTasks[n] = tp.spawn thread_func(n, c_len, cast[ptr UncheckedArray[float]](cs[0].addr), cast[ptr UncheckedArray[float]](avgheat[0].addr), cast[ptr UncheckedArray[float]](avgsus[0].addr), cast[ptr UncheckedArray[float]](avgcumul[0].addr), cast[ptr UncheckedArray[float]](avgM[0].addr))
  
  for n in 0 ..< times_run:
    discard sync pendingTasks[n]

  for i in 0 .. cs.high:
    avgheat[i] /= float(times_run)
    avgsus[i] /= float(times_run)
    avgcumul[i] /= float(times_run)
    avgM[i] /= float(times_run)

  let df1 = seqsToDf({"c":cs,"magnetization":avgM})
  let df2 = seqsToDf({"c":cs,"suseptibility":avgsus})
  let df3 = seqsToDf({"c":cs, "cumulant":avgcumul})
  let df4 = seqsToDf({"c":cs,"specific heat":avgheat})
  ggplot(df1, aes("c","magnetization")) + geom_point() +  ggsave("mag.png")
  ggplot(df2, aes("c","suseptibility")) + geom_point() +  ggsave("sus.png")
  ggplot(df3, aes("c","cumulant")) + geom_point() +  ggsave("cumul.png")
  ggplot(df4, aes("c","specific heat")) + geom_point() +  ggsave("heat.png")

main()