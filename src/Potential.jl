

function Potential_Energy(atom_0, Worm::Worm_=Worm,Particle::Particle_=Particle)
    Name_0 = Particle.Name[atom_0]
    offset_0 = Particle.Timeslices * atom_0

    E_p = 0.0

    for atom_i in 1:Number_Atoms
        if atom_i != atom_0
            Name_i = Particle.Name[atom_i]
            offset_i = Particle.Timeslices * atom_i
            for i_t in 1:Particle.Timeslices
                Potential_Energy_loop_i_t(Name_0, position, offset_0, Name_i, offset_i, i_t, Worm::Worm_,Particle::Particle_)
            end
        end
    end

    return E_p
end

function Potential_Energy_single(atom_0, position, i_t, Worm::Worm_=Worm,Particle::Particle_=Particle)
    Name_0 = Particle.Name[atom_0]
    offset_0 = Particle.Timeslices * atom_0

    E_p = 0.0

    for atom_i in 1:Number_Atoms
        if atom_i != atom_0
            Name_i = Particle.Name[atom_i]
            offset_i = Particle.Timeslices * atom_i
            Potential_Energy_loop_i_t(Name_0, position, offset_0, Name_i, offset_i, i_t, Worm::Worm_, Particle::Particle_)
        end
    end

    return E_p
end

function Potential_Energy_loop_i_t(Name_0, position, offset_0, Name_i, offset_i, i_t, Worm::Worm_=Worm,Particle::Particle_=Particle)
    Wordline = true
            
    if WORM_switch && Worm.exists && (Worm.Particle_Type == Name_i) 
        Wordline = Wordline_check((atom_i-Type_offset[Name_i]),i_)
    end

    dr² = 0.0
    if Wordline
        t_0 = offset_0 + i_t
        t_i = offset_i + i_t
        for dim in 1:Dimension 
            dr[dim] = position[dim,t_0] - Particle.Coords[dim][t_i]
            dr² += dr[dim]^2
        end
        r = sqrt(dr²)

        if Particle.Molcule_check[Name_0] || Particle.Molcule_check[Name_i]
            return
        else
            E_p += PES_1D(r,Name_i,Particle.Potential)
        end
        
    end
end

function Potential_rotation_Energy(atom_0,Cosine,i_t, Worm::Worm_=Worm,Particle::Particle_=Particle)
    Name_0 = Particle.Name[atom_0]
    offset_0 = Particle.Timeslices * atom_0
    E_p = 0.0

    for atom_i in 1:Number_Atoms
        if atom_i != atom_0
            Name_i = Particle.Name[atom_i]
            offset_i = Particle.Timeslices * atom_i
            Potential_rotation_Energy_loop_i_t(Cosine, Name_0, offset_0, Name_i, offset_i, i_t, Worm, Particle)
        end
    end

    return E_p
end

function Potential_rotation_Energy_loop_i_t(Cosine, Name_0, offset_0, Name_i, offset_i, i_t, Worm::Worm_=Worm,Particle::Particle_=Particle)
    Wordline = true
            
    if Worm.exists && (Worm.Name == Name_i) 
        Wordline = Wordline_check((atom_i-Particle.Type_offset[Name_i]),i_)
    end

    dr = Array{Vector{Int}}(undef, Dimension)

    dr² = 0.0
    if Wordline
        t_0 = offset_0 + i_t
        t_i = offset_i + i_t
        for dim in 1:Dimension 
            dr[dim] = Particle.Coords[dim,t_0] - Particle.Coords[dim][t_i]
            dr² += dr[dim]^2
        end
        r = sqrt(dr²)

        t_m = offset_0 + i_t/Rotation_Ratio
        Cos_t = 0.0
        for dim in 1:Dimension
            Cos_t += Cosine[dim,t_m] * dr[dim]
        end

        Cos_t /= -r
        
        E_rp +=  PES_2D(r,Cos_t,Name_0,Particle.Potential)
    end

    return E_rp
end

function Rotation_Density(γ,type,Potential=Particle.Potential)
    
end

function PES_1D(r,Name_i,Potential=Particle.Potential)
    PES_1D
end

function PES_2D(r,Cos_t,Name_0,Potential=Particle.Potential)
    
end
