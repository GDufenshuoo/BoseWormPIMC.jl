"""
This program is based on Moribs-PIMC
code by Zilong Wang <clockwang@icloud.com>
"""

using QuasiMonteCarlo, Distributions
using DelimitedFiles
using Interpolations

struct System_setting{T}
    Dimension::Int8
    Density::T
    Τ::T
    Block_step::Int
    Pass_perBlock::Int
end


struct Potential_{T}
    _PES
    Dimension::Int8
    scale_begin::Vector{T}
    scale_end::Vector{T}
    bin::Vector{T}
end

struct Particle_{T}
    λ::T
    β::T
    τ::T
    BoxSize::T
    Rotation_σ::T
    Rotation_Ratio::Int
    Rotation_Type::Int
    Rotation_Timeslices::Int
    Timeslices::Int
    Name::Vector{String}
    Step::Vector{T}
    multilevel::Vector{T}
    mass::Vector{T}
    Coords::Vector{T}
    Coords_forward::Vector{T}
    Angles::Vector{T}
    Type::Vector{Int}
    Number::Vector{Int}
    Type_offset::Vector{Int}
    T_wave²::Vector{Int}
    Particle_Index::Vector{Int}
    Ring_Index::Vector{Int}
    Potential::Vector{Potential_}
    Molcule_check::Vector{Bool}
end

function Init_MonteCarlo()
    setting = readdlm("MC.in", Float64; comments=true, comment_char='#')

    k = 0
    for i in setting[:,1]
        if i == "Atom" || i == "Molecule"
            k += 1
        elseif i == "Rotation"
            Particlle.Rotation_state = true
        end
    end

    Particle = Particle_(
        0.0,0.0,0.0,0.0,0.0,
        0,0,0,0,
        Array{String}(undef, k),
        Array{Vector{T}}(undef, k),
        Array{Vector{T}}(undef, k),
        Array{Vector{T}}(undef, k),
        Array{Vector{T}}(undef, k),
        Array{Vector{T}}(undef, k),
        Array{Vector{T}}(undef, k),
        Array{Vector{Int}}(undef, k),
        Array{Vector{Int}}(undef, k),
        Array{Vector{Int}}(undef, k),
        Array{Vector{Int}}(undef, k),
        Array{Vector{Int}}(undef, k),
        Array{Vector{Int}}(undef, k),
        Array{Vector{Potential_}}(undef, k),
        Array{Vector{Booll}}(undef, k)
        )

    j = 1
    k = 1
    for i in setting[:,1]
        if i == "Dimension:"
            System.Dimension = setting[j,2]
        elseif i == "Density"
            System.Density = setting[j,2]
        elseif i == "Temperature"
            System.Τ = setting[j,2]
        elseif i == "Timeslices"
            Particle.Timeslices = setting[j,2]
        elseif i == "Atom" || i == "Molecule"
            Particle.Name[k] = setting[j,2]
            Particle.Number[k] = setting[j,3]
            Particle.Type[k] = setting[j,4]
            Particle.Step[k] = setting[j,5]
            Particle.multilevel[k] = setting[j,6]
            Potential[k]=load(setting[j,7]+".pot")
            Particle.mass[i] = setting[j,8]
            Particle.Rotation_constant[i] = setting[j,9]
            k += 1
        elseif i == "Rotation"
            Particlle.Rotation_state = true
            Particle.Rotation_Type = setting[j,2]
            Particle.Rotation_Timeslices = setting[j,4]
            Particle.Rotation_Ratio = Int(Particle.Timeslices/Particle.Rotation_Timeslices)
            Particle.Rotation_step = setting[j,3]
        elseif i == "Worm"
            System.Worm_setting = [true,setting[j,2],setting[j,3]]
        elseif i == "Log"
            System.Block_total = setting[j,2]
            System.Pass_perBlock = setting[j,3]
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
    Particle.τ =  β/Particle.Timeslices
    Particle.Rotation_τ = Particle.β/Particle.Rotation_Timeslices
    for i in 1:k-1
        Particle.T_wave²[i] = 4*Particle.λ *Particle.τ
        Particle.multilevel_σ = Int(pow(2,Particle.multilevel[i]))
    end
    Particle.step = zeros(Int,Block_step,7)
    Particle.BoxSize = pow(sum(Particle.Number)/Particle.Density,Dimension)
    
    return Particle
end


function MonteCarlo_move(particle_type,particle::Particle_)
    for i in 1:Particle.Number[particle_type]
        particle_wordline = Particle.Type_offset[particle_type] + i
        for i_t in 1:Particle.Timeslices
            for dim in 1:Dimension
                disp = MonteCarlo_step[particle_type] * QcasiRand!(1,step)
                offset = particle_wordline * Particle.Timeslices
                Particle.Coords_forward[offset+i_t,dim] = Particle.Coords[offset+i_t,dim] + disp
            end            
        end
        Delta_E_p = Potential_Energy(particle_wordline,Particle.Coords_forward) - Potential_Energy(particle_wordline,Coords)
        MonteCarlo_accept(Delta_E_p,Particle.Coords_forward)
    end
end

"""
multilevel Metropolis
"""
function Bisection_Move(particle_type,time,Particle)
    for i in 1:Particle.Number[particle_type]
        offset = Particle.Type_offset
        gatom = Int(offset/Partilce.Timeslices)
        pit = Int((time + Particle.multilevel_σ[particle_type])%Particle.Timeslices)
        for dim in Particle.Dimension
            Particle.Coords_forward[offset+time,dim] = Particle.Coords[offset+time,dim]
            Particle.Coords_forward[offset+pit,dim] = Particle.Coords[offset+pit,dim]
        end
        bnorm = 1/(Particle.λ*particle.τ)
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
    if Delta_E < 0 || exp(-Delta_E*System.τ) > QcasiRand!(2,step)
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
    particle_wordline = Particle.Type_offset[particle_type]
    offset = particle_wordline * Particle.Timeslices
    for i_t in 1:Particle.Rotation_Timeslices
        i_t_0 = i_t - 1
        i_t_2 = i_t + 1
        Fix_range!(i_t_0,1,Particle.Rotation_Timeslices, Particle.Rotation_Timeslices)
        Fix_range!(i_t_2,1,Particle.Rotation_Timeslices, Particle.Rotation_Timeslices)
        t_0 = offset + i_t_0
        t_1 = offset + i_t
        t_2 = offset + i_t_2
        cosθ = Angles[t_1,2] + QcasiRand!(3,step)
        ϕ = Angles[t_1,1] + QcasiRand!(4,step)
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
            Rotation_Density =  Density_forword/Density_present * exp(-System.τ*(Density_forword-Density_present))
        else
            Rotation_Density = exp(-System.τ*(Density_forword-Density_present))
        end

        if Rotation_Density > 1 || Rotation_Density > QcasiRand!(5,step)
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




