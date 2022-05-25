

export 
Potential_Energy,
Potential_Energy_single,
Potential_rotation_Energy,
Rotation_Density


"""
after benchmark this part is not suitable for multi-threads
"""
function Potential_Energy(atom_0, position, Worm::Worm_=Worm,Particle::Particle_=Particle)
    # Name_0 = Particle.Name[atom_0]
    offset_0 = Particle.Timeslices * atom_0

    E_p = 0.0#Threads.Atomic{Float64}(0.0)

    # Threads.@threads 
    for atom_i in 1:Particle.Total_particle
        if atom_i != atom_0
            # Name_i = Particle.Name[atom_i]
            offset_i = Particle.Timeslices * (atom_i-1)
            # if WORM_switch && Worm.exists && (Worm.Particle_Type == Name_i) 
            #     Wordline = Wordline_check((atom_i-Type_offset[Name_i]),i_)
            # end
            # if Particle.Molcule_check[Name_0] || Particle.Molcule_check[Name_i]

            for i_t in 1:Particle.Timeslices
                # Threads.atomic_add!(E_p,Potential_Energy_loop_i_t(position, offset_0, offset_i, i_t, Worm::Worm_,Particle::Particle_))
                E_p+=Potential_Energy_loop_i_t(position, offset_0, offset_i, i_t, Worm::Worm_,Particle::Particle_)
            end
        end
    end

    return E_p[]
end

function Potential_Energy_single(atom_0, position, i_t, Worm::Worm_=Worm,Particle::Particle_=Particle)
    # Name_0 = Particle.Name[atom_0]
    offset_0 = Particle.Timeslices * (atom_0-1)

    E_p = 0.0

    for atom_i in 1:Particle.Total_particle
        if atom_i != atom_0
            # Name_i = Particle.Name[atom_i]
            offset_i = Particle.Timeslices * (atom_i-1)
            # if WORM_switch && Worm.exists && (Worm.Particle_Type == Name_i) 
            #     Wordline = Wordline_check((atom_i-Type_offset[Name_i]),i_)
            # end
            # if Particle.Molcule_check[Name_0] || Particle.Molcule_check[Name_i]
            E_p += Potential_Energy_loop_i_t(position, offset_0, offset_i, i_t, Worm::Worm_, Particle::Particle_)
        end
    end

    return E_p[]
end

function Potential_Energy_loop_i_t(position, offset_0, offset_i, i_t, Worm::Worm_=Worm,Particle::Particle_=Particle)
    Wordline = true
    E_p = 0.0
    if Wordline

        dr = zeros(Dimension)
        dr² = 0.0
        t_0 = offset_0 + i_t
        t_i = offset_i + i_t
        for dim in 1:Dimension 
            dr[dim] = position[t_0,dim] - Particle.Coords[t_i,dim]
            dr² += dr[dim]^2
        end
        r = sqrt(dr²)
        Name_i = 1
        # if Particle.Molcule_check[Name_0] || Particle.Molcule_check[Name_i]
            E_p = PES_1D(r,Name_i,Particle.Potential)
        # end
    end
    return E_p
end

function Potential_rotation_Energy(atom_0,Angle_Cosine,i_t, Worm::Worm_=Worm,Particle::Particle_=Particle)
    # Name_0 = Particle.Name[atom_0]
    offset_0 = Particle.Timeslices * atom_0
    E_p = 0.0
    for atom_i in 1:Particle.Total_particle
        if atom_i != atom_0
            # Name_i = Particle.Name[atom_i]
            offset_i = Particle.Timeslices * (atom_i-1)
            E_p += Potential_rotation_Energy_loop_i_t(Angle_Cosine, offset_0, offset_i, i_t, Worm, Particle)
        end
    end

    return E_p[]
end

function Potential_rotation_Energy_loop_i_t(Angle_Cosine, offset_0, offset_i, i_t, Worm::Worm_=Worm,Particle::Particle_=Particle)
    Wordline = true
            
    # if Worm.exists && (Worm.Name == Name_i) 
    #     Wordline = Wordline_check((atom_i-Particle.Type_offset[Name_i]),i_)
    # end

    E_rp = 0.0

    if Wordline
    dr² = 0.0
    dr = zeros(Dimension)
    # dr = Array{Int}(undef, Dimension)
        t_0 = offset_0 + i_t
        t_i = offset_i + i_t
        for dim in 1:Dimension 
            dr[dim] = Particle.Coords[t_0,dim] - Particle.Coords[t_i,dim]
            dr² += dr[dim]^2
        end
        Cos_t = 0.0
        r = sqrt(dr²)
        Cos_t /= -r
        if r == 0
            Cos_t = 0
        end

        i_t = Int(i_t ÷ Particle.Rotation_Ratio)
        for dim in 1:Dimension
            Cos_t += Angle_Cosine[i_t,dim] * dr[dim]
        end

        E_rp =  PES_2D(r,Cos_t,2,Particle.Potential)
    end

    return E_rp
end

function Rotation_Density(γ,type,Potential=Particle.Potential)
    rand()
end

function PES_1D(r,Name_i,Potential=Particle.Potential)
    Potential[1](r)
end

function PES_2D(r,Cos_t,Name_0,Potential=Particle.Potential)
    Potential[2](r,Cos_t)
end
