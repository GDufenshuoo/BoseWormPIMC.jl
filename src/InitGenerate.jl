"""
This program is based on Moribs-PIMC
code by Zilong Wang <clockwang@icloud.com>
"""

using DelimitedFiles
using Interpolations

export
    Init_MonteCarlo

function Init_MonteCarlo(MC_IN_PATH::String)
    setting = readdlm(MC_IN_PATH; comments=true, comment_char='#')

    System = System_setting(0,0.0,0.0,0,0,0)

    k = 0
    j = 0
    Timeslices = 0
    Total_particle = 0
    for i in setting[:,1]
        j+=1
        if i == "Dimension"
            System.Dimension = setting[j,2]
        elseif i == "Density"
            System.Density = setting[j,2]
        elseif i == "Temperature"
            System.Τ = setting[j,2]
        elseif i == "Atom" || i == "Molecule"
            k += 1
            Total_particle += setting[j,3]
        elseif i == "Log"
            System.Block_total = setting[j,2]
            System.Pass_perBlock = setting[j,3]
        elseif i == "Timeslices"
            Timeslices = setting[j,2]
        end
    end

    PES = Array{Any}(undef, k)

    Particle = Particle_(
        0.0,0.0,0.0,0.0,0.0,0.0,
        System.Dimension,0,0,Timeslices,0,Total_particle,
        Array{String}(undef, k),
        Array{Float64}(undef, k),
        Array{Float64}(undef, k),
        Array{Float64}(undef, k),
        Array{Float64}(undef, k),
        Array{Float64}(undef, k),
        zeros(Float64,Total_particle*Timeslices,3),
        zeros(Float64,Total_particle*Timeslices,3),
        zeros(Float64,Total_particle*Timeslices,3),
        zeros(Float64,Total_particle*Timeslices,3),
        zeros(System.Pass_perBlock,16),
        Array{Int}(undef, k),
        Array{Int}(undef, k),
        zeros(Int,k),
        Array{Int}(undef, k),
        collect(0:Total_particle-1),
        collect(0:Total_particle-1),
        Array{Int}(undef, k),
        Array{Int}(undef, 16),
        Array{Bool}(undef, Total_particle),
        false,
        Array{Any}(undef, k)
        )


    Bose_Femi = Dict("Bose"=>0,"Femi"=>1,"Boltz"=>-1)
    
    Worm_temp = [0,0,0.0]
    j = 1
    k = 0
    for i in setting[:,1]
        if i == "Atom" || i == "Molecule"
            k += 1
            Particle.Name[k] = setting[j,2]
            Particle.Number[k] = setting[j,3]
            Particle.Type[k] = Bose_Femi[setting[j,4]]
            Particle.Step[k] = setting[j,5]
            Particle.multilevel[k] = setting[j,6]

            Particle.Potential[k] = set_potention(load("../PES/"*setting[j,7], setting[j,8]))

            Particle.mass[k] = setting[j,9]
            Particle.λ[k] = 36.568596391365595/(2*setting[j,9])

            Particle.Rotation_constant[k] = setting[j,10]
            Particle.wordline[k] = Int(sum(Particle.Number[1:k])-Particle.Number[k])

            if i == "Molecule"
                Particle.Molcule_check[1+sum(Particle.Number)-setting[j,3]:sum(Particle.Number)] .= true
            else
                Particle.Molcule_check[1+sum(Particle.Number)-setting[j,3]:sum(Particle.Number)] .= false
            end
            

        elseif i == "Rotation"
            Particle.Rotation_Switch = true
            count = 0

            Particle.Rotation_Timeslices = setting[j,3]
            Particle.Rotation_Ratio = Int(Particle.Timeslices/Particle.Rotation_Timeslices)
            Particle.Rotation_step = setting[j,4]
        elseif i == "Worm"
            type = 0
            count = 0
            for l in Particle.Name
                count += 1
                if l == setting[j,2]
                    type = count
                end
                i
            end
            Worm_temp =[type,setting[j,3],setting[j,4]]
        elseif i == "xyzInit"
            if setting[j,2] == "Read"
                Particle.Coords = readdlm(setting[j,3], Float64; comments=true, comment_char='#')
            end
        elseif i == "log"
            System.Block_total = setting[j,2]
            System.Pass_perBlock = setting[j,3]
        elseif i == "#"
        end
        j += 1
    end


    Particle.β = 1/System.Τ
    Particle.τ =  Particle.β/Particle.Timeslices
    Particle.Rotation_τ = Particle.β/Particle.Rotation_Timeslices 

    for i in 1:k
        Particle.T_wave²[i] = 4*Particle.λ[i] *Particle.τ
        Particle.multilevel_σ[i] = Int(2^Particle.multilevel[i])
    end
    Particle.Angle_Cosine = Angle_Cosine_convert(2,Particle::Particle_)
    Particle.BoxSize = (Total_particle/System.Density)^System.Dimension
    Worm = Worm_Init(Int(Worm_temp[1]),Int(Worm_temp[2]),Worm_temp[3],Particle::Particle_,System::System_setting)
    return Particle,System,Worm
end

function set_potention(Pot)

    if Pot.Dimension == 1
        scale = (Pot.scale_begin[1]:Pot.bin[1]:Pot.scale_end[1])
    elseif Pot.Dimension == 2
        scale = (Pot.scale_begin[1]:Pot.bin[1]:Pot.scale_end[1], Pot.scale_begin[2]:Pot.bin[2]:Pot.scale_end[2])
    end

    return LinearInterpolation(
        scale,
        Pot._PES,
        extrapolation_bc=Line())
end
