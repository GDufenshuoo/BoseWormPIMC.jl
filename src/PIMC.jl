
using FileIO
using BenchmarkTools
using StatProfilerHTML
using TypedPolynomials


include("PI_MonteCarlo.jl")
include("Worm.jl")
include("Potential.jl")

Particle,System,Worm = Init_MonteCarlo()

global const Dimension = System.Dimension

for i in 1:10
    for time_slices in 1:Particle.Timeslices
        for particle_type in 1:size(Particle.Type)[1]
            MonteCarlo_move!(1,Particle)
            Bisection_Move!(2,time_slices,Particle)
            Worm_move!(1)
            MonteCarlo_Rotation_move!(2,Particle)
        end
    end
end
"""
10 × 32,16,2
BenchmarkTools.Trial: 9 samples with 1 evaluation.
 Range (min … max):  562.708 ms … 604.295 ms  ┊ GC (min … max): 7.34% … 6.87%
 Time  (median):     583.689 ms               ┊ GC (median):    7.26%
 Time  (mean ± σ):   583.202 ms ±  12.062 ms  ┊ GC (mean ± σ):  7.40% ± 0.38%

  ▁       ▁                   █ ▁     █  ▁                    ▁  
  █▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁█▁▁▁▁▁█▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  563 ms           Histogram: frequency by time          604 ms <

 Memory estimate: 564.73 MiB, allocs estimate: 18951370.
"""
