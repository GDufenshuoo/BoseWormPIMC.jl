
using FileIO

include("PI_MonteCarlo.jl")
include("Worm.jl")
include("Potential.jl")

Particle,System,Worm = Init_MonteCarlo()

global const Dimension = System.Dimension
MonteCarlo_move(1,Particle)
