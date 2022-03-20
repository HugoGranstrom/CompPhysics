import math, random, sequtils, ggplotnim, numericalnim, ginger, strformat, algorithm
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

const mcIterations = 500
const startSaving = 0
const saveEvery = 1

proc ising*(lattice: var Lattice, rnd: var Rand): (float, float, float, float, float, seq[float]) =
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
    if iters > startSaving and iters mod saveEvery == 0:
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
  return (M,msquare,m4,e,esquare, eTot)
  

proc calcProbability*(lattice: Lattice, row, col: int): float =
  let beta = 1 / (kb * lattice.T)
  let s = lattice[row-1, col] + lattice[row, col-1] + lattice[row, col+1] + lattice[row+1, col]
  let e = exp(2 * lattice.J * beta * s)
  result = e / (1 + e)

proc isingHeatBath*(lattice: var Lattice, rnd: var Rand): (float, float, float, float, float, seq[float]) =
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
    if iters > startSaving and iters mod saveEvery == 0:
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
  
  return (M,msquare,m4,e,esquare, eTot)

proc produceSpinPlots*(isingModel: proc (lattice: var Lattice, rnd: var Rand): (float, float, float, float, float, seq[float]), lattice: Lattice, cs: seq[float], fileName: string, title: string) =
  # Tile plot for different values of c
  var rnd = initRand(1337)
  var df = newDataFrame()
  var xs, ys: seq[int]
  for x in 0 ..< lattice.width-1:
    for y in 0 ..< lattice.height-1:
      xs.add x
      ys.add y
  df["x"] = xs
  df["y"] = ys
  df["Initial State"] = lattice.data
  var colNames = @["Initial State"]

  var tempLattice = lattice # reuse this so the end state of one iteration is used as start state of the next one.
  for c in cs:
    let T = tempLattice.J * c / kb
    tempLattice.T = T

    let _ = isingModel(tempLattice, rnd)
    let colName = fmt"c = {c:2.2f}"
    df[colName] = tempLattice.data
    colNames.add colName
  
  df = df.gather(colNames, value="spin", key="state")
  ggplot(df, aes("x", "y")) +
    scale_fill_discrete() +
    geom_tile(aes(fill="spin")) +
    scale_fill_gradient(magma()) +
    facet_wrap("state") +
    hideLegend() +
    ggtitle(title) +
    ggsave(filename)

proc plotEnergyIntervals*(lattice: Lattice, c: float, filename: string) =
  var df = newDataFrame()
  for (p, name) in [(ising, "metropolis"), (isingHeatBath, "heat bath")]:
    var startLattice = lattice
    let T = startLattice.J * c / kb
    startLattice.T = T
    var rnd = initRand(1337)
    let (M,msquare,m4,e,esquare, eTot) = p(startLattice, rnd)
    df[name] = eTot
  df["Iterations"] = linspace(startSaving.float, mcIterations.float, df.len)
  df = df.gather(["metropolis", "heat bath"], value="Energy", key="Method")
  ggplot(df, aes("Iterations", "Energy", color="Method")) +
    geom_line() +
    ggtitle("Energy against number of iterations for Metropolis and heat bath") +
    ggsave(filename)





when isMainModule:
  var rnd = initRand(1338)
  let startLattice = newLattice(32, 32, 0.000189574, 0.0, rnd)
  let cs = linspace(1.0, 5.0, 8).reversed
  produceSpinPlots(isingHeatBath, startLattice, cs, "heatbath.png", "Spin configurations for different c using heat bath")
  produceSpinPlots(ising, startLattice, cs, "metropolis.png", "Spin configurations for different c using Metropolis")
  plotEnergyIntervals(startLattice, 1.0, "metroVSheat.png")

  # Energy plot as function of iterations?


  
