
function MonteCarlo_move!(particle_type;Worm::Worm_=Worm,Particle::Particle_=Particle)
    pass = 0
    for atom in 0:Particle.Number[particle_type]-1
        particle_0 = Particle.wordline[particle_type] + atom
        offset = particle_0 * Particle.Timeslices

        @. Particle.Coords_forward[offset+1:offset+Particle.Timeslices,:] = Particle.Coords[offset+1:offset+Particle.Timeslices,:] + Particle.Step[particle_type]*(rand() - 0.5)

        Delta_E = Potential_Energy(particle_0,Particle.Coords_forward) - Potential_Energy(particle_0,Particle.Coords)

        if Delta_E <= 0 || exp(-Delta_E*Particle.τ) > rand()
            @. Particle.Coords[offset+1:offset+Particle.Timeslices,:] = Particle.Coords_forward[offset+1:offset+Particle.Timeslices,:]     
            pass += 1
        end
    end
    return pass
end

"""
multilevel Metropolis
"""
function Bisection_Move!(particle_type,time_slices;Worm::Worm_=Worm,Particle::Particle_=Particle)
    pass = 0
    for i in 0:Particle.Number[particle_type]-1
        gatom = Particle.wordline[particle_type] + i
        offset = gatom*Particle.Timeslices
        pit = Int((time_slices + Particle.multilevel_σ[particle_type])%Particle.Timeslices)

        for dim in 1:Dimension
            Particle.Coords_forward[offset+time_slices+1,dim] = Particle.Coords[offset+time_slices+1,dim]
            Particle.Coords_forward[offset+pit+1,dim] = Particle.Coords[offset+pit+1,dim]
        end

        bnorm = 1/(Particle.λ[particle_type]*Particle.τ)
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
                t_middle = Int(div((t_left+t_right),2))
                p_t_left    = (time_slices + t_left)   %Particle.Timeslices
                p_t_middle  = (time_slices + t_middle) %Particle.Timeslices
                p_t_right   = (time_slices + t_right)  %Particle.Timeslices
                for dim in 1:Dimension
                    Particle.Coords_forward[offset+p_t_middle+1,dim] = 0.5*(Particle.Coords_forward[offset+p_t_left+1,dim]+Particle.Coords_forward[offset+p_t_right+1,dim])
                    Particle.Coords_forward[offset+p_t_middle+1,dim] += Gaussian(bkin_norm)

                end
                E_p0 += Potential_Energy_single(gatom,Particle.Coords_forward,p_t_middle) - Potential_Energy_single(gatom,Particle.Coords,p_t_middle)
                if t_left != 1
                    E_p0 += Potential_Energy_single(gatom,Particle.Coords_forward,p_t_left) - Potential_Energy_single(gatom,Particle.Coords,p_t_left)
                end
                k = false
            end
            Δ_v = (E_p0 - 2*E_p1)*bpot_norm
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
function Bisection_Move_Exchange!(particle_type,time_slices;Worm::Worm_=Worm,Particle::Particle_=Particle)
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
            Particle.Coords_forward[offset_0+time_slices+1,dim] = Particle.Coords[offset_0+time_slices+1,dim]
            Particle.Coords_forward[offset_1+time_p+1,dim] = Particle.Coords[offset_1+time_p+1,dim]
        end

        bnorm = 1/(Particle.λ[particle_type]*Particle.τ)
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
                E_p0 += Potential_Energy_single(gatom_left,Particle.Coords_forward,p_t_middle) - Potential_Energy_single(gatom_left,Particle.Coords,p_t_middle)
                if t_left != 0
                    E_p0 += Potential_Energy_single(gatom_left,Particle.Coords_forward,p_t_left) - Potential_Energy_single(gatom_left,Particle.Coords,p_t_left)
                end
                k = false
            end
            Δ_v = (E_p0 - 2*E_p1)*bpot_norm
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

function MonteCarlo_accept!(Delta_E,offset;Worm::Worm_=Worm,Particle::Particle_=Particle)
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

function MonteCarlo_Rotation_move!(particle_type::Int;Worm::Worm_=Worm,Particle::Particle_=Particle)
    pass = 0

    particle_0 = Particle.wordline[particle_type]
    offset = particle_0 * Particle.Timeslices

    for i_t in 1:Particle.Rotation_Timeslices
        i_t_0 = i_t - 1
        i_t_2 = i_t + 1

        # This part is more or less too much to use a function
        i_t_0 = Fix_range(i_t_0,1,Particle.Rotation_Timeslices, Particle.Rotation_Timeslices)
        i_t_2 = Fix_range(i_t_2,1,Particle.Rotation_Timeslices, Particle.Rotation_Timeslices)
        t_1 = offset + i_t

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
            p_0 += Particle.Angle_Cosine[i_t_0,dim]*Particle.Angle_Cosine[i_t,dim]
            p_1 += Particle.Angle_Cosine[i_t,dim]*Particle.Angle_Cosine[i_t_2,dim]
        end
        
        Density_present = Rotation_Density(p_0,particle_type)*Rotation_Density(p_1,particle_type)

        if abs(Density_present)<1e-10
            Density_present = 0.0
        end
        i_t_r1 = i_t * Particle.Rotation_Ratio
        i_t_r0 = i_t_r1 - Particle.Rotation_Ratio

        E_p_present = 0.0
        for i_t_r in i_t_r0:i_t_r1-1
            E_p_present += Potential_rotation_Energy(particle_0,i_t_r)
        end
        p_0 = 0.0
        p_1 = 0.0
        @simd for dim in 1:Dimension
            p_0 += Particle.Angle_Cosine[i_t_0,dim]*Particle.Coords_forward[t_1,dim]
            p_1 += Particle.Coords_forward[t_1,dim]*Particle.Angle_Cosine[i_t_2,dim]
        end
        Density_forword = Rotation_Density(p_0,particle_type,Particle.Potential)*Rotation_Density(p_1,particle_type,Particle.Potential)
        if abs(Density_forword)<1e-10
            Density_forword = 0.0
        end

        E_p_forword = 0.0
        for i_t_r in i_t_r0:i_t_r1-1
            E_p_forword += Potential_rotation_Energy(particle_0,Particle.Coords_forward,i_t_r)
        end

        Density = 0.0
        if abs(Density_present)>1e-10
            Density =  Density_forword/Density_present * exp(-Particle.τ*(Density_forword-Density_present))
        else
            Density = 1.0
        end

        Density *= exp(-Particle.τ*(E_p_forword-E_p_present))
        if Density > 1 || Density > rand()
            Particle.Angles[t_1,2] = cosθ
            Particle.Angles[t_1,1] = ϕ

            for dim in 1:Dimension
                Particle.Angle_Cosine[i_t,dim] = Particle.Coords_forward[t_1,dim]
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
    range = Particle.Number[particle_type]*Particle.Timeslices#Particle.Rotation_Timeslices
    Angle_Cosine = zeros(range,Particle.Dimension)
    for i_t in 1:range
        cosθ = Particle.Angles[i_t+offset,2]
        ϕ = Particle.Angles[i_t+offset,1]
        sinθ = sqrt(1-cosθ^2)
        Angle_Cosine[i_t,1] = sinθ*cos(ϕ)
        Angle_Cosine[i_t,2] = sinθ*sin(ϕ)
        Angle_Cosine[i_t,3] = cosθ
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
    gaussian = sqrt(-log(rand()))*cos(2*pi*rand())/α
end
