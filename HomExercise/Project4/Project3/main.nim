import math, random, sequtils, ggplotnim, numericalnim, algorithm
import locks, taskpools, cpuinfo
import ../utils

const J = 0.000189574

var modLock: Lock

proc thread_func(task_id: int, c_len: int, cs: ptr UncheckedArray[float], avgHeat: ptr UncheckedArray[float], avgSus: ptr UncheckedArray[float], avgCumul: ptr UncheckedArray[float], avgM: ptr UncheckedArray[float],latticeZise: int): bool = 
  var rnd = initRand(task_id)
  let B = 0.0
  var lattice = newLattice(latticeZise, latticeZise, J, B, rnd)
  for i in 0 ..< c_len:
    let c = cs[i]

    let T = J*c/kb
    lattice.T = T
    
    #echo lattice.calcHamiltonian
    let (M,msquare,m4,e,esquare) = ising(lattice, rnd)
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

#SIZE 8
  echo "Running for size 8"
  let c_len = 40
  var avgheat1: seq[float] = newSeq[float](c_len)
  var avgsus1: seq[float] = newSeq[float](c_len)
  var avgcumul1: seq[float] = newSeq[float](c_len)
  var avgM1: seq[float] = newSeq[float](c_len)
  var cs = linspace(1.0,6.0,c_len).reversed
  let times_run = 6 # 550 = 7h
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