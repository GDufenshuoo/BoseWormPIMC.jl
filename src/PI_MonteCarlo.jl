"""
This program is based on Moribs-PIMC
code by Zilong Wang <clockwang@icloud.com>
"""

using QuasiMonteCarlo, Distributions
using DelimitedFiles
using Interpolations

mutable struct System_setting{T,K}
    Dimension::K
    Density::T
    Τ::T
    Block_step::K
    Block_total::K
    Pass_perBlock::K
end

struct Potential_{T}
    _PES
    Dimension::Int
    scale_begin::Vector{T}
    scale_end::Vector{T}
    bin::Vector{T}
end

mutable struct Particle_{T}
    λ::T
    β::T
    τ::T
    BoxSize::T
    Rotation_σ::T
    Rotation_τ::T
    Rotation_step::T
    Rotation_Ratio::Int
    Rotation_Timeslices::Int
    Timeslices::Int
    MonteCarlo_Step::Int
    Name::Vector{String}
    Step::Vector{T}
    mass::Vector{T}
    T_wave²::Vector{T}
    Rotation_constant::Vector{T}
    Coords::Matrix{T}
    Coords_forward::Matrix{T}
    Angles::Matrix{T}
    Quasi_rand::Matrix{T}
    multilevel::Vector{Int}
    Type::Vector{Int}
    Number::Vector{Int}
    wordline::Vector{Int}
    Particle_Index::Vector{Int}
    Ring_Index::Vector{Int}
    multilevel_σ::Vector{Int}
    Accept::Vector{Int}
    Molcule_check::Vector{Bool}
    Rotation_Switch::Bool
    Potential
end

function Init_MonteCarlo()
    setting = readdlm("MC.in"; comments=true, comment_char='#')

    k = 0
    j = 0
    Total_particle = 0
    for i in setting[:,1]
        j+=1
        if i == "Atom" || i == "Molecule"
            k += 1
            Total_particle += setting[j,3]
        elseif i == "Log"
            System.Block_total = setting[j,2]
            System.Pass_perBlock = setting[j,3]
        end
    end

    PES = Array{Any}(undef, k)

    Particle = Particle_(
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,
        0,0,0,0,
        Array{String}(undef, k),
        Array{Float64}(undef, k),
        Array{Float64}(undef, k),
        Array{Float64}(undef, k),
        Array{Float64}(undef, k),
        Array{Float64}(undef, Total_particle,3),
        Array{Float64}(undef, Total_particle,3),
        Array{Float64}(undef, Total_particle,3),
        zeros(System.Pass_perBlock,16),
        Array{Int}(undef, k),
        Array{Int}(undef, k),
        Array{Int}(undef, k),
        Array{Int}(undef, k),
        Array{Int}(undef, k),
        Array{Int}(undef, k),
        Array{Int}(undef, k),
        Array{Int}(undef, 16),
        Array{Bool}(undef, k),
        false,
        Array{Any}(undef, k)
        )

    System = System_setting(0,0.0,0.0,0,0,0)

    Bose_Femi = Dict("Bose"=>0,"Femi"=>1,"Boltz"=>-1)
    
    Worm_temp = [0,0,0.0]
    j = 1
    k = 0
    for i in setting[:,1]
        if i == "Dimension:"
        elseif i == "Density"
            System.Density = setting[j,2]
        elseif i == "Temperature"
            System.Τ = setting[j,2]
        elseif i == "Timeslices"
            Particle.Timeslices = setting[j,2]
        elseif i == "Atom" || i == "Molecule"
            k += 1
            Particle.Name[k] = setting[j,2]
            Particle.Number[k] = setting[j,3]
            Particle.Type[k] = Bose_Femi[setting[j,4]]
            Particle.Step[k] = setting[j,5]
            Particle.multilevel[k] = setting[j,6]
            # jldopen("example.jld2", "r")
            # Potential[k] = read(f,[setting[j,8]])
            Particle.Potential[k] = load(setting[j,7], setting[j,8])

            Particle.mass[k] = setting[j,9]
            Particle.Rotation_constant[k] = setting[j,10]

            Particle.wordline[k] = sum(Particle.Number[1:k]) - Particle.Number[k]
        elseif i == "Rotation"
            Particle.Rotation_Switch = true
            count = 0
            for l in Particle.Name
                count += 1
                if l == setting[j,2]
                    Particle.Molcule_check[count] = true
                else
                    Particle.Molcule_check[count] = false
                end
            end
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
        elseif i == "#"
        end
        j += 1
    end

    Particle.λ = 100*4.850873418093444e-19

    Particle.β = 1/System.Τ
    Particle.τ =  Particle.β/Particle.Timeslices
    Particle.Rotation_τ = Particle.β/Particle.Rotation_Timeslices

    for i in 1:k
        Particle.T_wave²[i] = 4*Particle.λ *Particle.τ
        Particle.multilevel_σ[i] = Int(2^Particle.multilevel[i])
    end
    # Particle.Particle.step = zeros(Int,Particle.Pass_perBlock,7)
    Particle.BoxSize = (Total_particle/System.Density)^System.Dimension
    Worm = Worm_Init(Int(Worm_temp[1]),Int(Worm_temp[2]),Worm_temp[3],Particle::Particle_,System::System_setting)
    return Particle,System,Worm
end


function MonteCarlo_move(particle_type,Particle::Particle_)
    for i in 1:Particle.Number[particle_type]
        particle_wordline = Particle.wordline[particle_type] + i
        for i_t in 1:Particle.Timeslices
            for dim in 1:Dimension
                disp = MonteCarlo_step[particle_type] * QcasiRand(1,Particle.Count)
                offset = particle_wordline * Particle.Timeslices
                Particle.Coords_forward[offset+i_t,dim] = Particle.Coords[offset+i_t,dim] + disp
            end            
        end
        Delta_E_p = Potential_Energy(particle_wordline,Particle.Coords_forward) - Potential_Energy(particle_wordline,Particle.Coords)
        MonteCarlo_accept(Delta_E_p,Particle)
    end
end

"""
multilevel Metropolis
"""
function Bisection_Move(particle_type,time,Particle)
    for i in 1:Particle.Number[particle_type]
        offset = Particle.wordline
        gatom = Int(offset/Partilce.Timeslices)
        pit = Int((time + Particle.multilevel_σ[particle_type])%Particle.Timeslices)
        for dim in Particle.Dimension
            Particle.Coords_forward[offset+time,dim] = Particle.Coords[offset+time,dim]
            Particle.Coords_forward[offset+pit,dim] = Particle.Coords[offset+pit,dim]
        end
        bnorm = 1/(Particle.λ*Particle.τ)
        for level in 0:Particle.multilevel[particle_type]-1
            level_σ = Int(2^(Particle.multilevel[particle_type]-level))
            bkin_norm = bnorm/level_σ
            bpot_norm = 0.5*Particle.τ*level_σ
            E_p0 = 0.0
            E_p1 = 0.0
            t_right = 0
            k = false
            while t_right < Particle.multilevel_σ[particle_type] && k == true
                t_left = t_right
                t_right = t_left + level_σ
                t_middle = Int((t_left+t_right)/2)
                p_t_left    = (time + t_left)   %Particle.Timeslices+1
                p_t_middle  = (time + t_middle) %Particle.Timeslices+1
                p_t_right   = (time + t_right)  %Particle.Timeslices+1
                for dim in 1:Particle.Dimension
                    Particle.Coords_forward[offset+p_t_middle,dim] = 0.5*(Particle.Coords_forward[offset+p_t_left,dim]+Particle.Coords_forward[offset+p_t_right,dim])
                    Particle.Coords_forward[offset+p_t_middle,dim] += gauss(bkin_norm)
                end
                E_p0 += Potential_Energy_single(gatom,Particle.Coords_forward,p_t_middle,Particle) - Potential_Energy_single(gatom,Particle.Coords,p_t_middle,Particle)
                if t_left != 1
                    E_p0 += Potential_Energy_single(gatom,Particle.Coords_forward,p_t_left,Particle) - Potential_Energy_single(gatom,Particle.Coords,p_t_left,Particle)
                end
                k = true
            end
            Δ_v = (E_p0 - 2*E_p1)*bpot_norm
            if Δ_v < 0 || exp(-Δ_v)>rand_3()
                for i_t in time:time+Particle.multilevel[particle_type]
                    for dim in 1:Dimension
                        p_i_t = i_t % Particle.Timeslices + 1
                        Particle.Coords[offset+p_i_t,dim] = Particle.Coords_forward[offset+p_i_t,dim]
                    end 
                end
            end
        end
    end
end

function MonteCarlo_accept(Delta_E,Particle::Particle_)
    if Delta_E < 0 || exp(-Delta_E*Particle.τ) > QcasiRand(2,Particle.Particle.Count)
        for i_t in 1:Particle.Timeslices
            for dim in 1:Dimension
                offset = particle_wordline * Particle.Timeslices
                Particle.Coords[offset+i_t,dim] = Particle.Coords_forward[offset+i_t,dim]
            end            
        end
        return true
    end
    return false
end

function MonteCarlo_Rotation_move(particle_type,Particle::Particle_)
    particle_wordline = Particle.wordline[particle_type]
    offset = particle_wordline * Particle.Timeslices
    for i_t in 1:Particle.Rotation_Timeslices
        i_t_0 = i_t - 1
        i_t_2 = i_t + 1
        Fix_range!(i_t_0,1,Particle.Rotation_Timeslices, Particle.Rotation_Timeslices)
        Fix_range!(i_t_2,1,Particle.Rotation_Timeslices, Particle.Rotation_Timeslices)
        t_0 = offset + i_t_0
        t_1 = offset + i_t
        t_2 = offset + i_t_2
        cosθ = Angles[t_1,2] + QcasiRand(3,Particle.Count)
        ϕ = Angles[t_1,1] + QcasiRand(4,Particle.Count)
        Fix_cos_range!(cosθ)
        sinθ = sqrt(1-cosθ^2)
        Particle.Coords_forward[t_0,1] = sinθ*cos(ϕ)
        Particle.Coords_forward[t_1,1] = sinθ*sin(ϕ)
        Particle.Coords_forward[t_2,1] = cosθ
        
        p_0 = 0.0
        p_1 = 0.0
        for dim in Dimension
            p_0 += Particle.Coords_forward[dim,t_0]*Particle.Coords_forward[dim,t_1]
            p_1 += Particle.Coords_forward[dim,t_1]*Particle.Coords_forward[dim,t_2]
        end
        Density_present = Rotation_Density(p_0,particle_type,Particle.Potential)*Rotation_Density(p_1,particle_type,Particle.Potential)
        if abs(Density_present)<1e-10
            Density_present = 0.0
        end
        i_t_r0 = i_t * RotRatio
        i_t_r1 = i_t_r0 + RotRatio - 1

        E_p_present = 0.0
        for i_t in i_t_r0:i_t_r1
            E_p_present += Potential_rotation_Energy(particle_wordline,Cosine,i_t)
        end

        p_0 = 0.0
        p_1 = 0.0
        for dim in Dimension
            p_0 += Particle.Coords_forward[dim,t_0]*Particle.Coords_forward[dim,t_1]
            p_1 += Particle.Coords_forward[dim,t_1]*Particle.Coords_forward[dim,t_2]
        end
        Density_forword = Rotation_Density(p_0,particle_type,Particle.Potential)*Rotation_Density(p_1,particle_type,Particle.Potential)
        if abs(Density_forword)<1e-10
            Density_forword = 0.0
        end
        i_t_r0 = i_t * RotRatio
        i_t_r1 = i_t_r0 + RotRatio

        E_p_forword = 0.0
        for i_t in i_t_r0:i_t_r1
            E_p_forword += Potential_rotation_Energy(particle_wordline,Particle.Coords_forward,i_t)
        end

        if abs(Density_present)>1e-10
            Rotation_Density =  Density_forword/Density_present * exp(-Particle.τ*(Density_forword-Density_present))
        else
            Rotation_Density = exp(-Particle.τ*(Density_forword-Density_present))
        end

        if Rotation_Density > 1 || Rotation_Density > QcasiRand(5,Particle.Count)
            Angles[2,t_1] = cosθ
            Angles[1,t_1] = ϕ

            for dim in Dimension
                MonteCarlo_Cosine[dim,t_1] = Particle.Coords_forward[dim,t_1]
            end
        end
    end
end

function Fix_range!(i,range_A,range,range_B)
    if i <= range_A
        i += range
    elseif i_t_2 > range_B
        i -= range
    end
end

function Fix_cos_range!(i)
    if i > 1 
        i = range - i
    elseif i < -1
        i = -range - i
    end
end




