using PyCall
@pyimport scipy.interpolate as si
@pyimport numpy
using Plots
using DelimitedFiles

function analysis_g2d(θrρ)
    θrρ
    s = size(θrρ)
    xyρ = zeros(Float64, s[1], s[2])

    x_trans!(xyρ, θrρ, s[1])

    one!(xyρ)
#    sort!(xyρ, dims = 1)
#    sort!(xyρ, dims = 1) error & big mistake
    return xyρ
end

function x_trans!(xyρ, θrρ, s)

    xyρ[:,1] = cos.((θrρ[:,1].-0.9)./180 .* pi).*θrρ[:,2]
    xyρ[:,2] = sin.((θrρ[:,1].-0.9)./180 .* pi).*θrρ[:,2]
    xyρ[:,3] = θrρ[:,3]./(θrρ[:,2].^2 .* sin.(θrρ[:,1]./180 .* pi))

end

function one!(xyρ)

    xyρ[:,3] =  (xyρ[:,3] ./  sum(xyρ[:,3]))
    return xyρ
end

# Location of random points to sample the function at
xmin = -15
xmax = 15
ymin = -1
ymax = 15

# Create a uniform grid to interpolate onto
nx = 300
ny = 160
xgrid = collect(numpy.linspace(xmin,xmax,nx))
ygrid = collect(numpy.linspace(ymin,ymax,ny))
grid_x = kron(ones(ny),xgrid')
grid_y = kron(ygrid,ones(1,nx))

# Perform the interpolation

gr()
x = range(0,100,step=1)
y = range(-70,70,step=1)

function con(x)
    if x <10
        return "00$x" 
    elseif x <100
        return "0$x"
    else
        return "$x"
    end    
end


function draw_plot(ρ,title_what)
    f(x,y) = ρ[y+14,x+150]
    # contour(y,x,f, fill = true, levels=5, clabels= 0.0025)
    # p = plot(y, x, f, st = [:contourf], title = "10H2 init 1-  BLOCK $i g2d file")
    plot(y, x, f, st = [:contourf], title = title_what, levels=100)
    # plot(y, x, f, st = [:contourf], title = title_what)
    # plot!(p[1])
end

function set_gird(xyρ)
    points = [xyρ[:,1] xyρ[:,2]]
    val = xyρ[:,3]
    grid_val = si.griddata(points,val,(grid_x,grid_y),method="cubic")
    return grid_val
end

function set_gird_diff(xyρ,xyρ2,beads,beads2)
    points = [xyρ[:,1] xyρ[:,2]]
    val = xyρ[:,3]
    grid_val = si.griddata(points,val,(grid_x,grid_y),method="cubic")
    points = [xyρ2[:,1] xyρ2[:,2]]
    val = xyρ2[:,3]
    grid_val2 = si.griddata(points,val,(grid_x,grid_y),method="cubic")
    # return grid_val.*beads .- (grid_val2.*beads2)
    return grid_val .- grid_val2
end

function set_gird_double_diff(xyρ1,xyρ12,xyρ,xyρ2,beads,beads2)
    set_gird_diff(xyρ1,xyρ12,beads,beads2).- set_gird_diff(xyρ,xyρ2,beads,beads2)
end

function draw(l,beads,ys="no") 
    i = 1
    frame = @animate while  i <= l
        println("draw    "*con(i))
        # xyρ = analysis_g2d("$beads-nw_bl/$(beads)_nw_bl_"*con(i)*".g2d")
        # xyρ = analysis_g2d("$(beads)_15/$(beads)__"*con(i)*".g2d")
        if 24 > i > 22 
            i=24
        end
        xyρ = analysis_g2d("0.15/$(beads)/$(beads)__"*con(i)*".g2d")
        # xyρ = analysis_g2d("9-local-s/0.5/local-s__"*con(i)*".g2d")
        ρ = set_gird(xyρ)
        # draw_plot(ρ,"$(beads) H2 $(ys)worm BLOCK$i")
        draw_plot(ρ,"$(beads) H2 $(ys)worm BLOCK$i")
        i += 1
    end
    gif(frame, "0.15/$(beads)_015K_g2d.gif", fps = 3)
end

function draw_f(l,beads,ys="no") 
    i = 1
@gif while  i <= l
    xyρ = analysis_g2d("$beads-f__"*con(i)*".g2d")
    ρ = set_gird(xyρ)
    println("draw    "*con(i))
    draw_plot(ρ,"$(beads)H2 f_$(ys)worm BLOCK$i")
    i += 1

end
end

function draw_f_2(beads,ys="no") 
    i=1
    @gif for opfile in readdir("$beads-f/4/")
        
        xyρ = analysis_g2d("$beads-f/4/" * opfile)
        ρ = set_gird(xyρ)

        println("draw    "*con(i))
        draw_plot(ρ,"$(beads)H2 f_$(ys)worm BLOCK$i")
        i += 1

    end
end

function draw_diff(l,beads) 
    i = 1
    b = 1
    xyρ = analysis_g2d("$beads-worm/$beads-worm__"*con(1)*".g2d") 
    xyρ2 = analysis_g2d("$beads-noworm/$beads-noworm__"*con(1)*".g2d")
    ρ = set_gird_diff(xyρ,xyρ2)
@gif while  i < l
    if (i<=20 && (b<(21-i)^1.2 +2)) || (i>=l-20 && b<=3) || (i>=l-7 && b<=7) || (i == l && b <= 30)
        b+=1
    else
        i += 1
        b = 1
        xyρ = analysis_g2d("$beads-worm/$beads-worm__"*con(i)*".g2d") 
        xyρ2 = analysis_g2d("$beads-noworm/$beads-noworm__"*con(i)*".g2d")
        ρ = set_gird_diff(xyρ,xyρ2)
    end 
    println("draw    "*con(i))
    draw_plot(ρ,"$beads H2  BLOCK$i")
end
end

function draw_single_diff(title,beads,beads2,ys1="",ys2="")
    xyρ = analysis_g2d("$beads-$(ys1)worm/$beads-$(ys1)worm___sum.g2d")
    xyρ2 = analysis_g2d("$beads2-$(ys2)worm/$beads2-$(ys2)worm___sum.g2d")
    ρ = set_gird_diff(xyρ,xyρ2,beads,beads2)
    # draw_plot(ρ,"$beads/$beads2 H2 $(ys1)/$(ys2)worm diff sum")
    draw_plot(ρ,title)
end


function draw_single_diff_f(title,beads,beads2,ys="")
    xyρ = analysis_g2d("$beads-$(ys)worm/$beads-$(ys)worm___sum.g2d") 
    xyρ2 = analysis_g2d("$beads2-$(ys)worm/$beads2-$(ys)worm___sum.g2d")
    ρ = set_gird_diff(xyρ,xyρ2,beads,beads2)
    # draw_plot(ρ,"$beads/$beads2 H2 $(ys)worm diff sum")
    draw_plot(ρ,title)
end



pyplot()

function Observe(g2d,BLOCK)
    xyρ = analysis_g2d(g2d)
    ρ = set_gird(xyρ)
    draw_plot(ρ,"$BLOCK")
    savefig("$BLOCK")
end


