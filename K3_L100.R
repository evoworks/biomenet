library(BiomeNet)
load("./results/network.RData")
newY = NULL
newZ = NULL
K = 3
L = 100
max.iter = 500
iter.lag = 20
burnin = 100
N = length(Compounds)
sample.results = .Call("doGibbs", S, R, 1:N, K, L, max.iter, newY, newZ, burnin, iter.lag, PACKAGE="BiomeNet")
save.image("./results/K3_L100.RData")
