
using FileIO
using BenchmarkTools
using StatProfilerHTML
using TypedPolynomials


include("PI_MonteCarlo.jl")
include("Worm.jl")
include("Potential.jl")

Particle,System,Worm = Init_MonteCarlo("MC.in")

global const Dimension = System.Dimension



# @benchmark for i in 1:10
        # for particle_type in 1:size(Particle.Type)[1]
        @benchmark MonteCarlo_move!(1,Particle)
        # @benchmark for time_slices in 1:Particle.Timeslices
        #     Bisection_Move!(2,time_slices,Particle)
        # end
        # @benchmark Worm_move!(1)
        # @benchmark MonteCarlo_Rotation_move!(2,Particle)
        # end
# end

