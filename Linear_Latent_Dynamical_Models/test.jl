import Confidence
reload("Confidence.jl")
using PyPlot

covar = [2 2;3 4]
mean = [1;2]

p1, p2 = Confidence.plot95(mean, covar)
plot(p1, p2)
