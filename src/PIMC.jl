
using FileIO
using BenchmarkTools
using StatProfilerHTML
using TypedPolynomials


include("PI_MonteCarlo.jl")
include("Worm.jl")
include("Potential.jl")

Particle,System,Worm = Init_MonteCarlo("MC.in")

global const Dimension = System.Dimension



@benchmark for i in 1:100
        # for particle_type in 1:size(Particle.Type)[1]
        # @benchmark MonteCarlo_move!(1,Particle)
        # @benchmark MonteCarlo_move!(2,Particle)
            for time_slices in 1:Particle.Timeslices
                Worm_move!(1)
                if time_slices == 1
                    MonteCarlo_move!(1,Particle)
                else
                    Bisection_Move!(2,time_slices,Particle)
                    MonteCarlo_Rotation_move!(2,Particle)
                end
            end

        end
# end

