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
    Dimension::Int
    Rotation_Ratio::Int
    Rotation_Timeslices::Int
    Timeslices::Int
    MonteCarlo_Step::Int
    Total_particle::Int
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

function Init_MonteCarlo(MC_IN_PATH::String)
    setting = readdlm(MC_IN_PATH; comments=true, comment_char='#')

    System = System_setting(0,0.0,0.0,0,0,0)

    k = 0
    j = 0
    Timeslices = 0
    Total_particle = 0
    for i in setting[:,1]
        j+=1
        if i == "Atom" || i == "Molecule"
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
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,
        System.Dimension,0,0,Timeslices,0,Total_particle,
        Array{String}(undef, k),
        Array{Float64}(undef, k),
        Array{Float64}(undef, k),
        Array{Float64}(undef, k),
        Array{Float64}(undef, k),
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
        if i == "Dimension"
            System.Dimension = setting[j,2]
        elseif i == "Density"
            System.Density = setting[j,2]
        elseif i == "Temperature"
            System.Τ = setting[j,2]
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

    Particle.λ = 100*4.850873418093444e-19

    Particle.β = 1/System.Τ
    Particle.τ =  Particle.β/Particle.Timeslices
    Particle.Rotation_τ = Particle.β/Particle.Rotation_Timeslices 

    for i in 1:k
        Particle.T_wave²[i] = 4*Particle.λ *Particle.τ
        Particle.multilevel_σ[i] = Int(2^Particle.multilevel[i])
    end
    # Particle.Particle.Step = zeros(Int,Particle.Pass_perBlock,7)
    Particle.BoxSize = (Total_particle/System.Density)^System.Dimension
    Worm = Worm_Init(Int(Worm_temp[1]),Int(Worm_temp[2]),Worm_temp[3],Particle::Particle_,System::System_setting)
    return Particle,System,Worm
end


function MonteCarlo_move!(particle_type,Particle::Particle_)
    #Threads.@threads 
    pass = 0
    for atom in 0:Particle.Number[particle_type]-1
        particle_0 = Particle.wordline[particle_type] + atom
        offset = particle_0 * Particle.Timeslices
        for dim in 1:Dimension
            for i_t in 1:Particle.Timeslices
                # disp = MonteCarlo_step[particle_type] * QcasiRand(1,Particle.Count)
                Particle.Coords_forward[offset+i_t,dim] = Particle.Coords[offset+i_t,dim] + Particle.Step[particle_type]*(rand() - 0.5)
            end
        end
        pass +=  MonteCarlo_accept!(
            (Potential_Energy(particle_0,Particle.Coords_forward) - Potential_Energy(particle_0,Particle.Coords)),
            offset,
            Particle)
    end
    return pass
end

"""
multilevel Metropolis
"""
function Bisection_Move!(particle_type,time_slices,Particle::Particle_)
    pass = 0
    # Threads.@threads 
    for i in 0:Particle.Number[particle_type]-1
        gatom = Particle.wordline[particle_type] + i
        offset = gatom*Particle.Timeslices
        pit = Int((time_slices + Particle.multilevel_σ[particle_type])%Particle.Timeslices)

        for dim in 1:Dimension
            Particle.Coords_forward[offset+time_slices,dim] = Particle.Coords[offset+time_slices,dim]
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
            k = true
            while t_right < Particle.multilevel_σ[particle_type] && k
                t_left = t_right
                t_right = t_left + level_σ
                t_middle = Int((t_left+t_right)/2)
                p_t_left    = (time_slices + t_left)   %Particle.Timeslices
                p_t_middle  = (time_slices + t_middle) %Particle.Timeslices
                p_t_right   = (time_slices + t_right)  %Particle.Timeslices
                for dim in 1:Dimension
                    Particle.Coords_forward[offset+p_t_middle,dim] = 0.5*(Particle.Coords_forward[offset+p_t_left,dim]+Particle.Coords_forward[offset+p_t_right,dim])
                    Particle.Coords_forward[offset+p_t_middle,dim] += Gaussian(bkin_norm)
                end
                E_p0 += Potential_Energy_single(gatom,Particle.Coords_forward,p_t_middle) - Potential_Energy_single(gatom,Particle.Coords,p_t_middle)
                if t_left != 1
                    E_p0 += Potential_Energy_single(gatom,Particle.Coords_forward,p_t_left) - Potential_Energy_single(gatom,Particle.Coords,p_t_left)
                end
                k = false
            end
            Δ_v = (E_p0 - 2*E_p1)*bpot_norm
            # if Δ_v < 0 || exp(-Δ_v)>rand_3()
            if Δ_v < 0 || exp(-Δ_v)>rand()
                for i_t in time_slices:time_slices+Particle.multilevel[particle_type]
                    for dim in 1:Dimension
                        p_i_t = i_t % Particle.Timeslices + 1
                        Particle.Coords[offset+p_i_t,dim] = Particle.Coords_forward[offset+p_i_t,dim]
                    end 
                end
                pass += 1 
            end
        end
    end
    return pass
end

"""
multilevel Metropolis
"""
function Bisection_Move_Exchange!(particle_type,time_slices,Particle::Particle_)
    pass = 0
    for particle in 0:Particle.Number[particle_type]-1
        offset_0 = (Particle.wordline[particle_type]+particle)*Particle.Timeslices
        offset_1 = offset_0

        time_1 = time_slices + Particle.multilevel_σ[particle_type]
        time_p = time_1 % Particle.Timeslices
        if time_p != time_1
            offset_1 = (Particle.wordline[particle_type] + Particle.Particle_Index[particle_type]) *Particle.Timeslices
        end

        for dim in 1:Dimension
            Particle.Coords_forward[offset_0+time_slices,dim] = Particle.Coords[offset_0+time_slices,dim]
            Particle.Coords_forward[offset_1+time_p+1,dim] = Particle.Coords[offset_1+time_p+1,dim]
        end

        bnorm = 1/(Particle.λ*Particle.τ)
        for level in 0:Particle.multilevel[particle_type]-1
            level_σ = Int(2^(Particle.multilevel[particle_type]-level))
            bkin_norm = bnorm/level_σ
            bpot_norm = 0.5*Particle.τ*level_σ
            E_p0 = 0.0
            E_p1 = 0.0
            t_right = 0
            k = true
            while t_right < Particle.multilevel_σ[particle_type] && k
                t_left = t_right
                t_right = t_left + level_σ
                t_middle = Int((t_left+t_right)/2)
                p_t_left    = (time_slices + t_left)   %Particle.Timeslices
                p_t_middle  = (time_slices + t_middle) %Particle.Timeslices
                p_t_right   = (time_slices + t_right)  %Particle.Timeslices
                offset_left = offset_0
                offset_middle = offset_0
                offset_right = offset_0
                if offset_left != time_slices + t_left
                    offset_left = offset_1
                elseif offset_middle != time_slices + t_left
                    offset_middle = offset_1
                elseif offset_right != time_slices + t_left
                    offset_right = offset_1
                end
                for dim in 1:Dimension
                    Particle.Coords_forward[offset_middle+p_t_middle+1,dim] = 0.5*(Particle.Coords_forward[offset_left+p_t_left+1,dim]+Particle.Coords_forward[offset_right+p_t_right+1,dim])
                    Particle.Coords_forward[offset_middle+p_t_middle+1,dim] += Gaussian(bkin_norm)
                end
                gatom_left = Int(offset_left/Particle.Timeslices)
                # gatom_middle = Int(offset_middle/Particle.Timeslices)
                E_p0 += Potential_Energy_single(gatom_left,Particle.Coords_forward,p_t_middle) - Potential_Energy_single(gatom_left,Particle.Coords,p_t_middle)
                if t_left != 1
                    E_p0 += Potential_Energy_single(gatom_left,Particle.Coords_forward,p_t_left) - Potential_Energy_single(gatom_left,Particle.Coords,p_t_left)
                end
                k = false
            end
            Δ_v = (E_p0 - 2*E_p1)*bpot_norm
            # if Δ_v < 0 || exp(-Δ_v)>rand_3()
            if Δ_v <= 0 || exp(-Δ_v)>rand()
                for i_t in time_slices:time_slices+Particle.multilevel[particle_type]
                    p_i_t = (i_t % Particle.Timeslices) + 1
                    offset = offset_0
                    if p_i_t != i_t
                        offset = offset_1
                    end
                    for dim in 1:Dimension
                         Particle.Coords[offset+p_i_t,dim] = Particle.Coords_forward[offset+p_i_t,dim]
                    end 
                end
                pass += 1
            end
        end
    end
    return pass
end

function MonteCarlo_accept!(Delta_E,offset,Particle::Particle_)
    # if Delta_E < 0 || exp(-Delta_E*Particle.τ) > QcasiRand(2,Particle.Particle.Count)
    if Delta_E <= 0 || exp(-Delta_E*Particle.τ) > rand()
        for i_t in 1:Particle.Timeslices
            for dim in 1:Dimension
                Particle.Coords[offset+i_t,dim] = Particle.Coords_forward[offset+i_t,dim]
            end            
        end
        return 1
    else
        return 0
    end
end

function MonteCarlo_Rotation_move!(particle_type::Int,Angle_Cosine,Particle::Particle_)
    pass = 0

    particle_wordline = Particle.wordline[particle_type]
    offset = particle_wordline * Particle.Timeslices

    for i_t in 1:Particle.Rotation_Timeslices
        i_t_0 = i_t - 1
        i_t_2 = i_t + 1

        # This part is more or less too much to use a function
        i_t_0 = Fix_range(i_t_0,1,Particle.Rotation_Timeslices, Particle.Rotation_Timeslices)
        i_t_2 = Fix_range(i_t_2,1,Particle.Rotation_Timeslices, Particle.Rotation_Timeslices)

        # t_0 = offset + i_t_0
        t_1 = offset + i_t
        # t_2 = offset + i_t_2
        # cosθ = Angles[t_1,2] + QcasiRand(3,Particle.Count)
        # ϕ = Angles[t_1,1] + QcasiRand(4,Particle.Count)
        cosθ = Particle.Angles[t_1,2] + Particle.Step[particle_type]*(rand() -0.5)
        ϕ = Particle.Angles[t_1,1] + Particle.Step[particle_type]*(rand() -0.5)
        cosθ = Fix_cos_range(cosθ)
        sinθ = sqrt(1-cosθ^2)
        Particle.Coords_forward[t_1,1] = sinθ*cos(ϕ)
        Particle.Coords_forward[t_1,2] = sinθ*sin(ϕ)
        Particle.Coords_forward[t_1,3] = cosθ

        p_0 = 0.0
        p_1 = 0.0
        @simd for dim in 1:Dimension
            p_0 += Angle_Cosine[i_t_0,dim]*Angle_Cosine[i_t,dim]
            p_1 += Angle_Cosine[i_t,dim]*Angle_Cosine[i_t_2,dim]
        end
        
        Density_present = Rotation_Density(p_0,particle_type)*Rotation_Density(p_1,particle_type)

        if abs(Density_present)<1e-10
            Density_present = 0.0
        end
        i_t_r1 = i_t * Particle.Rotation_Ratio
        i_t_r0 = i_t_r1 + Particle.Rotation_Ratio

        E_p_present = 0.0
        for i_t_r in i_t_r0:i_t_r1-1
            E_p_present += Potential_rotation_Energy(particle_wordline,Angle_Cosine,i_t_r)
        end

        p_0 = 0.0
        p_1 = 0.0
        @simd for dim in 1:Dimension
            p_0 += Angle_Cosine[i_t_0,dim]*Particle.Coords_forward[t_1,dim]
            p_1 += Particle.Coords_forward[t_1,dim]*Angle_Cosine[i_t_2,dim]
        end
        Density_forword = Rotation_Density(p_0,particle_type,Particle.Potential)*Rotation_Density(p_1,particle_type,Particle.Potential)
        if abs(Density_forword)<1e-10
            Density_forword = 0.0
        end
        # i_t_r0 = i_t * Particle.Rotation_Ratio
        # i_t_r1 = i_t_r0 + Particle.Rotation_Ratio

        E_p_forword = 0.0
        for i_t_r in (i_t_r0+1):i_t_r1
            E_p_forword += Potential_rotation_Energy(particle_wordline,Particle.Coords_forward,i_t_r)
        end

        Density = 0.0
        if abs(Density_present)>1e-10
            Density =  Density_forword/Density_present * exp(-Particle.τ*(Density_forword-Density_present))
        else
            Density = 1.0
        end

        Density *= exp(-Particle.τ*(E_p_forword-E_p_present))
        # if Rotation_Density > 1 || Rotation_Density > QcasiRand(5,Particle.Count)
        if Density > 1 || Density > rand()
            Particle.Angles[t_1,2] = cosθ
            Particle.Angles[t_1,1] = ϕ

            for dim in 1:Dimension
                Angle_Cosine[i_t,dim] = Particle.Coords_forward[t_1,dim]
            end
            pass += 1
        end
    end
    return pass
end

function Fix_range(i,range_A,range,range_B)
    if i < range_A
        i += range
    elseif i > range_B
        i -= range
    end
    return i
end

function Fix_cos_range(i)
    if i > 1 
        i = 2 - i
    elseif i < -1
        i = -2 - i
    end
    return i
end

function Angle_Cosine_convert(particle_type,Particle::Particle_)
    offset = Particle.wordline[particle_type] * Particle.Timeslices
    range = Particle.Number[particle_type]*Particle.Rotation_Timeslices
    Angle_Cosine = zeros(range,Dimension)
    @simd for i_t in 1:range
        cosθ = Particle.Angles[i_t+offset,2]
        ϕ = Particle.Angles[i_t+offset,1]
        sinθ = sqrt(1-cosθ^2)
        i_t_r = fld(i_t,Particle.Rotation_Ratio)+1
        Angle_Cosine[i_t_r,1] = sinθ*cos(ϕ)
        Angle_Cosine[i_t_r,2] = sinθ*sin(ϕ)
        Angle_Cosine[i_t_r,3] = cosθ
    end
    return Angle_Cosine
end

"""
@benchmark for i in 1:100000
    Gaussian(rand())
end

BenchmarkTools.Trial: 1264 samples with 1 evaluation.
Range (min … max):  3.601 ms …  11.334 ms  ┊ GC (min … max): 0.00% … 0.00%
Time  (median):     3.892 ms               ┊ GC (median):    0.00%
Time  (mean ± σ):   3.948 ms ± 389.519 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

    ▃▅ ▂▃▄▆▇█▄▄▆▄▃ ▁                                        
▄▇▄▅▇████████████████▇▆▆▇▅▅▄▃▃▃▃▃▃▃▃▂▃▃▃▃▃▁▃▂▃▂▃▂▂▂▁▁▂▂▁▂▁▂ ▄
3.6 ms          Histogram: frequency by time        4.84 ms <

Memory estimate: 0 bytes, allocs estimate: 0.
"""
function Gaussian(α)
    sqrt(-log(rand()))*cos(2*pi*rand())/α
end


