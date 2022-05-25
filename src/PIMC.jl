
using FileIO
using BenchmarkTools
using StatProfilerHTML
using TypedPolynomials


include("PI_MonteCarlo.jl")
include("Worm.jl")
include("Potential.jl")

Particle,System,Worm = Init_MonteCarlo("MC.in")

global const Dimension = System.Dimension



@benchmark for i in 1:10
    for time_slices in 1:Particle.Timeslices
        # for particle_type in 1:size(Particle.Type)[1]
            MonteCarlo_move!(1,Particle)
            Bisection_Move!(2,time_slices,Particle)
            Worm_move!(1)
            MonteCarlo_Rotation_move!(2,Particle)
        # end
    end
end

