import math, random, sequtils, ggplotnim, numericalnim
import locks, taskpools, cpuinfo

type
  Lattice* = object
    width*, height*: int
    data*: seq[float]
    J*, B*, T*: float 

const kb* = 8.617e-5

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

proc newLattice*(height, width: int, J, B: float, rnd: var Rand): Lattice =
  result.height = height
  result.width = width
  result.J = J
  result.B = B
  result.data = newSeqWith[float]((height-1)*(width-1), rnd.sample([-1.0, 1.0]))


proc calcHamiltonian*(lattice: Lattice, row, col: int): float =
  let current = lattice[row, col]
  #result += current * lattice[row-1, col] #
  result += current * lattice[row+1, col]
  #result += current * lattice[row, col-1] #
  result += current * lattice[row, col+1]

proc calcHamiltonian*(lattice: Lattice): float =
  result = -lattice.B * sum(lattice.data)
  var pairsSum: float
  for row in 0 ..< lattice.height - 1:
    for col in 0 ..< lattice.width - 1:
      pairsSum += calcHamiltonian(lattice, row, col)
  result += -lattice.J * pairsSum

proc calcHamiltonianFull(lattice: Lattice, row, col: int): float =
  let current = lattice[row, col]
  result += current * lattice[row-1, col] #
  result += current * lattice[row+1, col]
  result += current * lattice[row, col-1] #
  result += current * lattice[row, col+1]
  result = -lattice.B * current + -lattice.J * result

proc M_calc*(lattice: Lattice): float =
  return sum(lattice.data)/lattice.data.len.toFloat

proc Heat_calc*(E: float, eSquare: float, T: float): float =
  return (eSquare - E*E)/(kb*T*T)

proc Sus_calc*(M: float, mSquare: float, T:float): float =
  return (mSquare-M*M)/(kb*T)

proc Cumulant*(mSquare: float, m4: float): float =
  return 1 - m4/(3*mSquare*mSquare)

proc mean*(x: seq[float]): float =
  return sum(x)/x.len.float

const mcIterations = 10000

proc ising*(lattice: var Lattice, rnd: var Rand): (float, float, float, float, float) =
  var hamiltonian = calcHamiltonian(lattice)
  var dE = 1e10
  var eTot: seq[float]
  var mTot: seq[float]
  var msquareTot: seq[float]
  var esquareTot: seq[float]
  var m4Tot: seq[float]
  for iters in 0 .. mcIterations:
    #var innerHamiltonian = hamiltonian
    for row in 0 ..< lattice.height - 1:
      for col in 0 ..< lattice.width - 1:
        let energyBefore = calcHamiltonianFull(lattice, row, col)
        lattice.flip(row, col)
        let energyDiff = calcHamiltonianFull(lattice, row, col) - energyBefore
        #let newHamiltonian = calcHamiltonian(lattice)
        #let dEinner = newHamiltonian - innerhamiltonian
        let dEinner = energyDiff
        #echo "Full: ", deInner, " Local: ", energyDiff, " Rel err: ", (energyDiff - deInner) / deInner
        if dEinner > 0:
          let r = rnd.rand(1.0)
          #if iters < 10: echo -dEinner / (lattice.T * kb)
          if r > exp(-dEinner / (lattice.T * kb)):
            #echo "Rejected! ", dEinner
            lattice.flip(row, col) # flip back
    let innerHamiltonian = lattice.calcHamiltonian
    if iters > 3000 and iters mod 100 == 0:
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
  

proc calcProbability*(lattice: Lattice, row, col: int): float =
  let beta = 1 / (kb * lattice.T)
  let s = lattice[row-1, col] + lattice[row, col-1] + lattice[row, col+1] + lattice[row+1, col]
  let e = exp(2 * lattice.J * beta * s)
  result = e / (1 + e)

proc isingHeatBath*(lattice: var Lattice, rnd: var Rand): (float, float, float, float, float) =
  var hamiltonian = calcHamiltonian(lattice)
  var dE = 1e10
  var eTot: seq[float]
  var mTot: seq[float]
  var msquareTot: seq[float]
  var esquareTot: seq[float]
  var m4Tot: seq[float]
  for iters in 0 .. mcIterations:
    for row in 0 ..< lattice.height - 1:
      for col in 0 ..< lattice.width - 1:
        let p = calcProbability(lattice, row, col)
        let r = rnd.rand(1.0)
        if r < p:
          lattice[row, col] = 1
        else:
          lattice[row, col] = -1
    let innerHamiltonian = calcHamiltonian(lattice)
    if iters > 3000 and iters mod 100 == 0:
      let innerHamiltonian = calcHamiltonian(lattice)
      let M = M_calc(lattice)
      let msquare = M*M
      let m4 = msquare*msquare
      let esquare = innerHamiltonian*innerHamiltonian
      eTot.add innerHamiltonian
      mTot.add M
      msquareTot.add msquare
      esquareTot.add esquare
      m4Tot.add m4
    
    dE = innerHamiltonian - hamiltonian
    
    hamiltonian = innerHamiltonian
  let M = mean(mTot)
  let msquare = mean(msquareTot)
  let m4 = mean(m4Tot)
  let e = mean(eTot)
  let esquare = mean(esquareTot)
  
  return (M,msquare,m4,e,esquare)