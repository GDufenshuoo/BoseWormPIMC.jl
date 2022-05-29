using Plots

function Observe_2d!(g2d,Particle::Particle_=Particle)
    E = 0
    for Particle_0 in 0:Particle.Total_particle-2
        offset_0 = Particle.Timeslices * Particle_0
        relation_currect = Particle.Molcule_check[Particle_0+1]
        for particle in 0:Particle.Total_particle-1
            Molcule_check = Particle.Molcule_check[particle+1] || relation_currect
            offset_i = Particle.Timeslices * particle
            if  Molcule_check
                for i_t in 1:Particle.Timeslices
                    Observe_2d_loop_i_t!(g2d, offset_0, offset_i, i_t, relation_currect)
                end
            else
                # Observe_g1d_loop_i_t!(g1d, offset_0, offset_i, i_t, relation_currect)
            end
        end
    end
end


function Observe_2d_loop_i_t!(g2d, offset_0::Int, offset_i::Int, i_t::Int, relation_currect::Bool, Particle::Particle_=Particle)

    dr = zeros(Float32,Dimension)
    dr² = 0.0
    t_0 = offset_0 + i_t
    t_i = offset_i + i_t
    for dim in 1:Dimension 
        dr[dim] = Particle.Coords[t_0,dim] - Particle.Coords[t_i,dim]
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
    # println(Cos_t)
    if r == 0
        Cos_t = 0
    end
    if 14 > r > 2 && 14 > Cos_t > -14
        g2d[Int(div((r-2),0.024))+1,Int(div(Cos_t+1,0.005))+1] += 1E-5
    end
end

function Observe_Plots(g2d,BLOCK::Int;threads::Int=-1)
    g2d ./= sum(g2d)
    itp = LinearInterpolation((2:0.024:14,-1:0.005:1),g2d,extrapolation_bc=0)
    f(x,y) = itp(hypot(x,y),x/hypot(x,y))
    if threads == -1
        the_title = "BLOCK $BLOCK"
    else
        the_title = "threads $threads  BLOCK $BLOCK"
    end
    plot(-7:0.05:7,0:0.05:7, f, st = [:contourf], title = the_title, levels=25)
end
    # x = -7:0.05:7
    # y = 0:0.05:7
    # X = repeat(reshape(x, 1, :), length(y), 1)
    # Y = repeat(y, 1, length(x))
    # Z = map(f, X, Y)
    # plot(-7:0.05:7,0:0.05:7, f, st = [:contourf], title = the_title, levels=25)
    # contour(x, y, Z)
    # contour(-7:0.05:7,0:0.05:7, f, fill = true)