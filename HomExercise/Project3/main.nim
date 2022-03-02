import math, random, sequtils

type
  Lattice* = object
    width*, height*: int
    data*: seq[float]
    J*, B*, T*: float 

const kb = 8.617e-5

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

proc ising*(lattice: var Lattice, steps: int) =
  var hamiltonian = calcHamiltonian(lattice)
  var dE = 1e10
  var iters: int
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
          innerhamiltonian = newHamiltonian
    dE = lattice.calcHamiltonian() - hamiltonian
    #echo dE
    hamiltonian = innerHamiltonian

  echo "Iterations: ", iters

proc main() =
  let J = 0.000189574
  let T = 2.2
  let B = 0.0
  randomize(10)
  var lattice = newLattice(10, 10, J, B, T)
  echo lattice.calcHamiltonian
  ising(lattice, 100)
  echo lattice.calcHamiltonian

  
main()

# c = kT/J (1 → 4)
# t = 2.2 - 2.4
# J = k*(2.2 → 2.4) / (1 → 4)
# k * 2.2 / 1 → k*2.4 / 4 = k * 2.2 → k * 0.6

