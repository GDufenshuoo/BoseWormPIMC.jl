module PIMC

export pimc

using JLD2

    include("initMC.jl")
    include("PI_MonteCarlo.jl")
    include("Worm.jl")
    include("Potential.jl")
    include("Observe.jl")

    function pimc(MC)
        global Particle, System, Worm = Init_MonteCarlo(MC)

        global Dimension = System.Dimension

        begin
            r = (Particle.Total_particle-1)*Particle.Timeslices
            a = collect(0:Particle.Number[1]-1)
            b = a.*2pi./Particle.Number[1]
            c = 0.2*rand(Particle.Number[1])
            local j = 1
            for i in 1:Particle.Number[1]
                Particle.Coords[(i-1)*Particle.Timeslices+1:i*Particle.Timeslices,1] .= 5 *cos.(b[j])
                Particle.Coords[(i-1)*Particle.Timeslices+1:i*Particle.Timeslices,2] .= 5 *sin.(b[j])
                Particle.Coords[(i-1)*Particle.Timeslices+1:i*Particle.Timeslices,3] .= 5 *c[j]
                j += 1
            end
        end

        for i in 1:100

            MonteCarlo_move!(1)
            for i in 0:Particle.Timeslices-1
                Bisection_Move!(1,i)
                # MonteCarlo_move!(2,Particle)
                Bisection_Move!(2,i)
            end

        end

        println("Begin BLOCK")
        for BLOCK in 1:System.Block_total
            Worm_pass = zeros(Int,10)
            Bisection_step = zeros(Int,2)
            Bisection_pass = zeros(Int,2)
            MonteCarlo_step = 0
            MonteCarlo_pass = 0
            Rotation_step = 0
            Rotation_pass = 0
            g1d = zeros(Int64,151)
            g2d = zeros(Int64,51,41)
            count_Bisection_1 = ((2^Particle.multilevel[1]+1)-1)*Particle.Number[1]
            count_Bisection_2 = ((2^Particle.multilevel[1]+1)-1)*Particle.Number[1]
            for pass in 1:System.Pass_perBlock
                for time_slices in 0:Particle.Timeslices-1
                    # global Worm_pass .+= Worm_move!(1,Worm,Particle)
                    if !Worm.state
                        Bisection_step[1] += count_Bisection_1
                        # global Bisection_pass[1] += Bisection_Move_Exchange!(1,rand(0:Particle.Timeslices-1))
                    end
                    if time_slices == 0
                        # MonteCarlo_pass += MonteCarlo_move!(1,Particle)
                        MonteCarlo_step += Particle.Number[2]
                        MonteCarlo_pass += MonteCarlo_move!(2)
                    else
                        Bisection_step[2] += count_Bisection_2
                        # global Bisection_pass[2] += Bisection_Move!(2,time_slices)
                        Rotation_step += Particle.Rotation_Timeslices
                        # global Rotation_pass += MonteCarlo_Rotation_move!(2,Angle_Cosine)
                    end
                    if !Worm.state
                        Observe_2d!(g1d,g2d,Particle)
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

            Observe_Plots(g2d,BLOCK)
            savefig("H2_worm_sum_gr_$BLOCK")
        end
        jldsave("Backup";Worm,Particle)

    end

end
