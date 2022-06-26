"""
This program is based on Moribs-PIMC
code by Zilong Wang <clockwang@icloud.com>
"""

using Plots

function Observe_2d!(g1d,g2d,Particle::Particle_=Particle)
    for Particle_0 in 1:Particle.Total_particle-1
        offset_0 = Particle.Timeslices * Particle_0
        relation_currect = Particle.Molcule_check[Particle_0+1]
        for particle in 0:Particle_0-1
            offset_i = Particle.Timeslices * particle
            Molcule_check = Particle.Molcule_check[particle+1] || relation_currect
            if  Molcule_check
                for i_t in 0:Particle.Timeslices-1
                    Observe_2d_loop_i_t!(g2d, offset_0, offset_i, i_t, relation_currect)
                end
            else
                for i_t in 0:Particle.Timeslices-1
                    Observe_1d_loop_i_t!(g1d, offset_0, offset_i, i_t, relation_currect)
                end
            end 
        end
    end
end

function Observe_1d_loop_i_t!(g1d, offset_0::Int, offset_i::Int, i_t::Int, relation_currect::Bool, Particle::Particle_=Particle)

    dr = zeros(Float32,Dimension)
    dr² = 0.0
    t_0 = offset_0 + i_t
    t_i = offset_i + i_t
    for dim in 1:Dimension 
        dr[dim] = Particle.Coords[t_0+1,dim] - Particle.Coords[t_i+1,dim]
        dr² += dr[dim]^2
    end
    if dr² < 1
        println(t_0,t_i)
    end
    
    r = sqrt(dr²)

    if r < 15
        g1d[Int(div(r,0.1))+1] += 1
    end
end

function Observe_2d_loop_i_t!(g2d, offset_0::Int, offset_i::Int, i_t::Int, relation_currect::Bool, Particle::Particle_=Particle)

    dr = zeros(Float32,Dimension)
    dr² = 0.0
    t_0 = offset_0 + i_t
    t_i = offset_i + i_t
    for dim in 1:Dimension 
        dr[dim] = Particle.Coords[t_0+1,dim] - Particle.Coords[t_i+1,dim]
        dr² += dr[dim]^2
    end
    r = sqrt(dr²)
    
    sign = 1
    t_m = fld(offset_i + i_t,Particle.Rotation_Ratio)
    if relation_currect
        sign = -1
        t_m = fld(offset_0 + i_t,Particle.Rotation_Ratio)
    end

    cosθ = Particle.Angles[t_m+1,2]
    ϕ = Particle.Angles[t_m+1,1]
    sinθ = sqrt(1-cosθ^2)
    Cos_t = (sinθ*cos(ϕ)*dr[1] + sinθ*sin(ϕ)*dr[2] + cosθ*dr[3])/(sign*r)
    if r == 0
        Cos_t = 0
    end
    if 14 > r > 2 && 14 > Cos_t > -14
        g2d[Int(div((r-2),0.24))+1,Int(div(Cos_t+1,0.05))+1] += 1
    end
end

function Observe_Plots(g2d,BLOCK::Int;threads::Int=-1)
    temp = g2d/sum(g2d)
    itp = LinearInterpolation((2:0.24:14,-1:0.05:1),temp,extrapolation_bc=0)
    f(x,y) = itp(hypot(x,y),x/hypot(x,y))
    if threads == -1
        the_title = "BLOCK $BLOCK"
    else
        the_title = "threads $threads  BLOCK $BLOCK"
    end
    heatmap(-7:0.2:7,0:0.2:7, f)
end