

export 
Potential_Energy,
Potential_Energy_single,
Potential_rotation_Energy,
Rotation_Density


"""
after benchmark this part is not suitable for multi-threads
"""
function Potential_Energy(Particle_0, position, Worm::Worm_=Worm,Particle::Particle_=Particle)
    # Name_0 = Particle.Name[atom_0]
    offset_0 = Particle.Timeslices * Particle_0

    E_p = Threads.Atomic{Float64}(0.0)

    relation_currect = Particle.Molcule_check[Particle_0+1]
    Threads.@threads for particle_type in 1:size(Particle.wordline)[1]
        Particle_number = Particle.Number[particle_type]
        Molcule_check = Particle.Molcule_check[Particle_number+1] || relation_currect

        Threads.@threads for particle in 0:Particle_number
            # if particle != Particle_0 && in_Wordline 

                offset_i = Particle.Timeslices * particle
                # if Worm.states && (Worm.Particle_Type == Name_i) 
                #     in_Wordline = Wordline_check((atom_i-Type_offset[Name_i]),i_)
                # end

                if  Molcule_check
                    for i_t in 1:Particle.Timeslices
                        Threads.atomic_add!(E_p,Potential_Energy_multi_loop_i_t(position, offset_0, offset_i, i_t, relation_currect))
                        # E_p+=Potential_Energy_loop_i_t(position, offset_0, offset_i, i_t, Worm::Worm_,Particle::Particle_)
                    end
                else
                    for i_t in 1:Particle.Timeslices
                        Threads.atomic_add!(E_p,Potential_Energy_loop_i_t(position, offset_0, offset_i, i_t))
                        # E_p+=Potential_Energy_loop_i_t(position, offset_0, offset_i, i_t, Worm::Worm_,Particle::Particle_)
                    end
                end

            # end
        end
    end

    return E_p[]
end

function Potential_Energy_single(Particle_0, position, i_t, Worm::Worm_=Worm,Particle::Particle_=Particle)
    offset_0 = Particle.Timeslices * Particle_0

    E_p = Threads.Atomic{Float64}(0.0)

    relation_currect = Particle.Molcule_check[Particle_0+1]
    Threads.@threads for particle_type in 1:size(Particle.wordline)[1]
        Particle_number = Particle.Number[particle_type]
        Molcule_check = Particle.Molcule_check[Particle_number+1] || relation_currect

        for particle in 0:Particle_number
            # if particle != Particle_0 && in_Wordline 
                offset_i = Particle.Timeslices * particle
                # if Worm.states && (Worm.Particle_Type == Name_i) 
                #     in_Wordline = Wordline_check((atom_i-Type_offset[Name_i]),i_)
                # end

                if  Molcule_check
                        Threads.atomic_add!(E_p,Potential_Energy_multi_loop_i_t(position, offset_0, offset_i, i_t, relation_currect))
                        # E_p+=Potential_Energy_loop_i_t(position, offset_0, offset_i, i_t, Worm::Worm_,Particle::Particle_)
                else
                        Threads.atomic_add!(E_p,Potential_Energy_loop_i_t(position, offset_0, offset_i, i_t))
                        # E_p+=Potential_Energy_loop_i_t(position, offset_0, offset_i, i_t, Worm::Worm_,Particle::Particle_)
                end
            # end
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
        E_p = PES_1D(r)
    end

    return E_p
end

function Potential_Energy_multi_loop_i_t(position, offset_0, offset_i, i_t, relation_currect, Worm::Worm_=Worm,Particle::Particle_=Particle)
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
        
        sign = 1
        t_m = Int(offset_i + i_t ÷ Particle.Rotation_Ratio)
        if relation_currect
            sign = -1
            t_m = Int(offset_0 + i_t ÷ Particle.Rotation_Ratio)
        end
        Cos_t = 0.0

        cosθ = Particle.Angles[t_m+1,2]
        ϕ = Particle.Angles[t_m+1,1]
        sinθ = sqrt(1-cosθ^2)
        Cos_t += sinθ*cos(ϕ)
        Cos_t += sinθ*sin(ϕ)
        Cos_t += cosθ

        Cos_t /= sign*r
        if r == 0
            Cos_t = 0
        end

        E_p = PES_2D(r,Cos_t)
        
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

    return E_p
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

        E_rp = PES_2D(r,Cos_t,Particle.Potential)
    end

    return E_rp
end

function Rotation_Density(γ,type,Potential=Particle.Potential)
    rand()
end

function PES_1D(r,Potential=Particle.Potential)
    Potential[1](r)
end

function PES_2D(r,Cos_t,Potential=Particle.Potential)
    Potential[2](r,Cos_t)
end
