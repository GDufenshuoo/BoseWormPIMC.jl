
using FileIO
using BenchmarkTools
using StatProfilerHTML
using TypedPolynomials


include("PI_MonteCarlo.jl")
include("Worm.jl")
include("Potential.jl")

# @benchmark Potential_Energy(2, Particle.Coords, Worm,Particle)
# for particle_type in 1:size(Particle.Type)[1:1]
    # for i in 1:100
    #     MonteCarlo_pass[1] += MonteCarlo_move!(1,Particle)
    #     MonteCarlo_pass[2] += MonteCarlo_move!(2,Particle)
    # end

"""
BenchmarkTools.Trial: 3 samples with 1 evaluation.
 Range (min … max):  2.216 s …   2.292 s  ┊ GC (min … max): 2.61% … 3.51%
 Time  (median):     2.280 s              ┊ GC (median):    2.54%
 Time  (mean ± σ):   2.262 s ± 40.688 ms  ┊ GC (mean ± σ):  2.83% ± 0.60%

  █                                              █        █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁█ ▁
  2.22 s         Histogram: frequency by time        2.29 s <

 Memory estimate: 507.78 MiB, allocs estimate: 23930696.
"""
function PIMC(input)
    Particle,System,Worm = Init_MonteCarlo(string(input))

    global Dimension = System.Dimension

    Worm_pass = zeros(Int,10)
    Bisection_step = zeros(Int,2)
    Bisection_pass = zeros(Int,2)
    MonteCarlo_step = 0
    MonteCarlo_pass = 0
    Rotation_step = 0
    Rotation_pass = 0

    for i in 1:10
        for i in 1:100
            for time_slices in 1:Particle.Timeslices
                """
                BenchmarkTools.Trial: 10000 samples with 6 evaluations.
                Range (min … max):  3.556 μs …  1.019 ms  ┊ GC (min … max): 0.00% … 99.02%
                Time  (median):     6.636 μs              ┊ GC (median):    0.00%
                Time  (mean ± σ):   7.431 μs ± 20.162 μs  ┊ GC (mean ± σ):  5.35% ±  1.98%

                    ██▅▂                                                      
                ▂▇████▇▆▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▆▇▆▅▅▄▅▆▅▄▄▄▄▄▃▃▂▂▂▂▂▂▃▃▂▂▁▁▁▁▁▁▁▁ ▃
                3.56 μs        Histogram: frequency by time        14.3 μs <

                Memory estimate: 2.22 KiB, allocs estimate: 102.
                """
                Worm_pass .+= Worm_move!(1)
                if !Worm.state
                    Bisection_step[1] += Particle.multilevel[2]*Particle.Number[2]
                    """
                    BenchmarkTools.Trial: 10000 samples with 1 evaluation.
                    Range (min … max):  190.916 μs …   5.690 ms  ┊ GC (min … max): 0.00% … 95.83%
                    Time  (median):     197.861 μs               ┊ GC (median):    0.00%
                    Time  (mean ± σ):   210.072 μs ± 189.963 μs  ┊ GC (mean ± σ):  3.22% ±  3.44%

                    ▆██▇▆▄▃▂▁▁▁                                                  ▂
                    ██████████████▇▇▇▆▇█▇█▇▇▆▇▇▆▆▆▆▆▅▆▅▆▇▇▆▆▆▅▆▅▅▃▆▅▅▅▅▅▄▂▅▄▄▄▄▄▅ █
                    191 μs        Histogram: log(frequency) by time        312 μs <

                    Memory estimate: 76.06 KiB, allocs estimate: 3563.
                    """
                    Bisection_pass[1] += Bisection_Move_Exchange!(1,rand(1:Particle.Timeslices),Particle)
                end
                if time_slices == 1
                    # MonteCarlo_pass[1] += MonteCarlo_move!(1,Particle)
                    MonteCarlo_step += Particle.Number[2]
                    """
                    BenchmarkTools.Trial: 10000 samples with 1 evaluation.
                    Range (min … max):  26.873 μs …  4.551 ms  ┊ GC (min … max): 0.00% … 99.14%
                    Time  (median):     28.084 μs              ┊ GC (median):    0.00%
                    Time  (mean ± σ):   30.345 μs ± 77.160 μs  ┊ GC (mean ± σ):  4.39% ±  1.72%

                    ▁▇█▇▆▆▆▆▅▅▄▃▂▂▂▁▁▁▁                                         ▂
                    ████████████████████▇██▇▇▇▇▆▇█▇▇▆▇█▇▅▆▇▇▆▅▇▇▇▇▂▅▇▇▇▆▅▃▄▄▇▇▆ █
                    26.9 μs      Histogram: log(frequency) by time        42 μs <

                    Memory estimate: 18.20 KiB, allocs estimate: 1001.
                    """
                    MonteCarlo_pass += MonteCarlo_move!(2,Particle)
                else
                Bisection_step[2] += Particle.multilevel[2]*Particle.Number[2]
                """
                BenchmarkTools.Trial: 10000 samples with 1 evaluation.
                Range (min … max):  81.663 μs …   5.488 ms  ┊ GC (min … max): 0.00% … 97.88%
                Time  (median):     83.730 μs               ┊ GC (median):    0.00%
                Time  (mean ± σ):   90.996 μs ± 147.991 μs  ┊ GC (mean ± σ):  4.57% ±  2.77%

                ▇█▇▆▄▃▂▁▁▁   ▁    ▁                                          ▂
                ███████████████▇▆███▆█▇▆▆██▅▇█▇▆▇▆▅▄▆▆▅▅▆▆▆▅▅▇▅▄▅▃▅▃▁▄▅▅▆▅▄▄ █
                81.7 μs       Histogram: log(frequency) by time       152 μs <

                Memory estimate: 46.27 KiB, allocs estimate: 2449.
                """
                Bisection_pass[2] += Bisection_Move!(2,time_slices,Particle)
                Rotation_step += Particle.Rotation_Timeslices
                """
                BenchmarkTools.Trial: 10000 samples with 7 evaluations.
                Range (min … max):  4.320 μs … 781.533 μs  ┊ GC (min … max): 0.00% … 98.83%
                Time  (median):     4.729 μs               ┊ GC (median):    0.00%
                Time  (mean ± σ):   5.672 μs ±  15.026 μs  ┊ GC (mean ± σ):  3.94% ±  1.71%

                ▁▇█▇▅▃▃▂▁▁                   ▁▁▁                            ▁
                ████████████▇▇▅▇██▇█▇▇▇██▇▇▇███████▇▇▇▆▆▅▆▅▆▄▅▆▅▅▄▆▅▅▅▅▅▅▅▄ █
                4.32 μs      Histogram: log(frequency) by time      12.9 μs <

                Memory estimate: 2.32 KiB, allocs estimate: 132.
                """
                Rotation_pass += MonteCarlo_Rotation_move!(2,Particle)
                end
            end
        end
        println(Particle.Name[1])
        println("Worm   total:  open ", Worm_pass[1], "   close ", Worm_pass[3], "    advance ", Worm_pass[5], "   recede ", Worm_pass[7], "     swap", Worm_pass[9])
        println("Worm   pass:   open ", Worm_pass[2], "   close ", Worm_pass[4], "    advance ", Worm_pass[6], "   recede ", Worm_pass[8], "     swap", Worm_pass[10])
        println("Worm   radio:  open ", Float16(Worm_pass[1]\Worm_pass[2]), "   close ", Float16(Worm_pass[3]\Worm_pass[4]), "    advance ", Float16(Worm_pass[5]\Worm_pass[6]), "   recede ", Float16(Worm_pass[7]\Worm_pass[8]), "     swap ", Float16(Worm_pass[9]\Worm_pass[10]))
        println("Bisection  Total:  ", Bisection_step[1])
        println("Bisection  pass:   ", Bisection_pass[1])
        println("Bisection  radio:  ", Bisection_step[1]\Bisection_pass[1])
        println(Particle.Name[2])
        println("MonteCarlo Total:  ", MonteCarlo_step)
        println("MonteCarlo pass:   ", MonteCarlo_pass)
        println("MonteCarlo radio:  ", MonteCarlo_step\MonteCarlo_pass)
        println("Bisection  Total:  ", Bisection_step[2])
        println("Bisection  pass:   ", Bisection_pass[2])
        println("Bisection  radio:  ", Bisection_step[2]\Bisection_pass[2])
        println("Rotation   Total:  ", Rotation_step)
        println("Rotation   pass:   ", Rotation_pass)
        println("Rotation   radio:  ", Rotation_step\Rotation_pass)

    end
end

# end

