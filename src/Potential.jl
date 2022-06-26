"""
This program is based on Moribs-PIMC
code by Zilong Wang <clockwang@icloud.com>
"""

export 
Potential_Energy,
Potential_Energy_single,
Potential_rotation_Energy,
Rotation_Density


"""
after benchmark this part is not suitable for multi-threads
"""
function Potential_Energy(particle_0::Int, Pos; Worm::Worm_=Worm,Particle::Particle_=Particle)
    offset_0 = Particle.Timeslices * particle_0

    E_p = 0.0

    relation_currect = Particle.Molcule_check[particle_0+1]
    for particle_type in 1:size(Particle.wordline)[1]

        Particle_number = Particle.Number[particle_type]
        wordline_type = Particle.wordline[particle_type]
        Molcule_check = Particle.Molcule_check[wordline_type+1] || relation_currect

        for particle in 0:Particle_number-1
            if particle + wordline_type != particle_0 
                offset_i = Particle.Timeslices * (particle + wordline_type)

                Worm_check = false
                if Worm.state && (Worm.type == particle_type)
                    Worm_check = true
                end
                if  Molcule_check
                    E_p += Potential_Energy_2d_loop_i_t(Pos, offset_0, offset_i, relation_currect)
                else
                    if Worm_check
                        E_p += Potential_Energy_loop_i_t(Pos, offset_0, offset_i)
                    else
                        E_p += Potential_Energy_loop_i_t(Pos, offset_0, offset_i, Worm::Worm_,Particle::Particle_)
                    end
                end
            end
        end
    end
    return E_p
end

function Potential_Energy_single(particle_0::Int, Pos, i_t::Int; Worm::Worm_=Worm,Particle::Particle_=Particle)
    offset_0 = Particle.Timeslices * particle_0
    i_t += 1
    E_p = 0.0
    relation_currect = Particle.Molcule_check[particle_0+1]
    for particle_type in 1:size(Particle.wordline)[1]
        Particle_number = Particle.Number[particle_type]
        wordline_type = Particle.wordline[particle_type]
        Molcule_check = Particle.Molcule_check[wordline_type+1] || relation_currect

        for particle in 0:Particle_number-1
            if particle + wordline_type != particle_0 
                offset_i = Particle.Timeslices * (particle + wordline_type)
                Worm_check = false
                if Worm.state && (Worm.type == particle_type)
                    Worm_check = true
                end
                if  Molcule_check
                    dc = Pos[offset_0 + i_t,:].-Particle.Coords[offset_i + i_t,:]
                    r = R(dc[1],dc[2],dc[3])
                
                    sign = 1
                    if relation_currect
                        sign = -1
                    end
    
                    Cos_t = sum(Particle.Angle_Cosine[i_t,:].*dc)/r
                    
                    E_p += sum(PES_2D.(r,Cos_t))
                else
                    if Worm_check && Wordline_check(particle-Particle.wordline[particle_type],i_t)
                        dc = Pos[offset_0 + i_t,:].-Particle.Coords[offset_i + i_t,:]
                        E_p += PES_1D(R(dc[1],dc[2],dc[3]))
                    else
                        dc = Pos[offset_0 + i_t,:].-Particle.Coords[offset_i + i_t,:]
                        E_p += PES_1D(R(dc[1],dc[2],dc[3]))
                    end
                end 
            end

        end
    end
    return E_p
end

function Potential_Energy_loop_i_t(Pos::Matrix{Float64}, offset_0::Int, offset_i::Int; Worm::Worm_=Worm,Particle::Particle_=Particle)
    dc = .-(Pos[(offset_0 +1):(offset_0 +Particle.Timeslices),:],Particle.Coords[(offset_i +1):(offset_i +Particle.Timeslices),:])
    return sum(PES_1D(R.(dc[:,1],dc[:,2],dc[:,3])))
end

function Potential_Energy_2d_loop_i_t(Pos::Matrix{Float64}, offset_0::Int, offset_i::Int,relation_currect::Bool; Worm::Worm_=Worm,Particle::Particle_=Particle)
    dc = .-(Pos[(offset_0 +1):(offset_0 +Particle.Timeslices),:],Particle.Coords[(offset_i +1):(offset_i +Particle.Timeslices),:])
    r = R.(dc[:,1],dc[:,2],dc[:,3])

    sign = 1

    if relation_currect
        sign = -1
    end

    Cos_t = sum(Particle.Angle_Cosine.*dc./r,dims=2)

    return sum(PES_2D.(r,Cos_t))
end

function Potential_rotation_Energy(atom_0,i_t; Worm::Worm_=Worm,Particle::Particle_=Particle)
    offset_0 = Particle.Timeslices * atom_0
    in_Wordline = true
    E_p = 0.0
    for particle_type in 1:size(Particle.wordline)[1]
        for atom_i in 1:Particle.Number[particle_type]
            if atom_i != atom_0
                if Worm.state && (Worm.type == particle_type) 
                    in_Wordline = Wordline_check((atom_i-Type_offset[Name_i]),particle)
                end
                if in_Wordline

                    offset_i = Particle.Timeslices * (atom_i-1)
                    E_p += Potential_rotation_Energy_loop_i_t( offset_0, offset_i, i_t, Worm, Particle)
                end
            end
        end
    end
    return E_p
end

function Potential_rotation_Energy_loop_i_t( offset_0, offset_i, i_t; Worm::Worm_=Worm,Particle::Particle_=Particle)
    E_rp = 0.0
    dr² = 0.0
    dr = Array{Float64}(undef, Dimension)

    t_0 = offset_0 + i_t +1
    t_i = offset_i + i_t +1
    for dim in 1:Dimension 
        dr[dim] = Particle.Coords[t_0,dim] - Particle.Coords[t_i,dim]
        dr² += dr[dim]^2
    end
    
    Cos_t = 0
    i_t = fld(i_t,Particle.Rotation_Ratio)+1
    for dim in 1:Dimension
        Cos_t += Particle.Angle_Cosine[i_t,dim] * dr[dim]
    end
    r = sqrt(dr²)
    Cos_t /= -r
    if dr² == 0
        Cos_t = 0
    end
    E_rp = PES_2D(r,Cos_t,Particle.Potential)
    println(E_rp)
    return E_rp
end


function Rotation_Density(γ,type,Potential=Particle.Potential)
    0
end

function PES_1D(r,Potential=Particle.Potential)
    Potential[1](r)
end

function PES_2D(r,Cos_t,Potential=Particle.Potential)
    Potential[2](r,Cos_t)
end

function i_Times(A_C,R)
    A_C*R
end

function Cost(sinθ,cosθ,ϕ,sign,r,dc)
    @. (sinθ*cos(ϕ)*dc[:,1] + sinθ.*sin(ϕ)*dc[:,2] + cosθ*dc[:,3])/sign*r
end

function CosSinTrans(cosθ)
    sqrt(1-cosθ^2)
end

function dx(A,B)
    A-B
end

function R(x,y,z)
    sqrt(x^2+y^2+z^2)
end

