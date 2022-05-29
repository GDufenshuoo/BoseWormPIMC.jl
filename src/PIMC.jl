
using FileIO
using BenchmarkTools


include("PI_MonteCarlo.jl")
include("Worm.jl")
include("Potential.jl")
include("Observe.jl")

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

Particle, System, Worm = Init_MonteCarlo("MC.in")

global Dimension = System.Dimension

Worm_pass = zeros(Int,10)
Bisection_step = zeros(Int,2)
Bisection_pass = zeros(Int,2)
MonteCarlo_step = 0
MonteCarlo_pass = 0
Rotation_step = 0
Rotation_pass = 0

g2d = zeros(Float32,501,401)
# d2d = zeros(Float32,501,401)
Angle_Cosine = Angle_Cosine_convert(2,Particle::Particle_)

for i in 1:100
    MonteCarlo_move!(1,Particle)
    MonteCarlo_move!(2,Particle)
end

println("Begin BLOCK")
for BLOCK in 1:System.Block_total
    for pass in 1:System.Pass_perBlock
        for time_slices in 1:Particle.Timeslices
            global Worm_pass .+= Worm_move!(1,Worm,Particle)
            if !Worm.state
                global Bisection_step[1] += Particle.multilevel[2]*Particle.Number[2]
                global Bisection_pass[1] += Bisection_Move_Exchange!(1,rand(1:Particle.Timeslices),Particle)
            end
            if time_slices == 1
                # MonteCarlo_pass[1] += MonteCarlo_move!(1,Particle)
                global MonteCarlo_step += Particle.Number[2]
                global MonteCarlo_pass += MonteCarlo_move!(2,Particle)
            else
            global Bisection_step[2] += Particle.multilevel[2]*Particle.Number[2]
            global Bisection_pass[2] += Bisection_Move!(2,time_slices,Particle)
            global Rotation_step += Particle.Rotation_Timeslices
            global Rotation_pass += MonteCarlo_Rotation_move!(2,Angle_Cosine,Particle)
            end
        end
        Observe_2d!(g2d,Particle)
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

    Observe_Plots(g2d,BLOCK)
    savefig("H2_worm_sum_gr")
end
return g2d
# end

