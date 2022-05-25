

export 
    Worm_Init,
    Worm_move

mutable struct Worm_{T}
    ira::Int
    masha::Int
    type::Int
    atom_i::Int
    atom_m::Int
    Timeslices::Int
    m::Int
    c::T 
    Cutoff²::T
    _norm::T
    _dr²::Vector{T}
    _atom::Vector{T}
    neighbor::Matrix{T}
    _Particle_table::Vector{Int}
    state::Bool
end

function Worm_Init(type,worm_m,worm_c,Particle::Particle_,System::System_setting)
    Worm = Worm_(
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0.0,
        0.0,
        0.0,
        zeros(Particle.Number[type]),
        zeros(Particle.Number[type]),
        zeros(Particle.Number[type],2),
        Array{Int}(undef, Particle.Number[type]),
        false
        )

    worm_Cutoff = 100

    Worm.m = worm_m
    Worm.Timeslices = Particle.Timeslices

    Worm.c = worm_c*System.Density/Particle.Number[type]*Worm.Timeslices*Worm.m
    Worm._norm = Worm.c*Worm.Timeslices*Worm.m
    Worm.Cutoff² = worm_Cutoff^2 * Worm.m*4*Particle.λ[type]*Particle.τ

    return Worm
end

function Worm_move!(particle_type, Worm::Worm_=Worm, Particle::Particle_=Particle)
    for i in 1:Particle.Number[particle_type]
        if Worm.state
            Worm_open!(Worm::Worm_,Particle::Particle_)
        else
            Worm_close!(Worm::Worm_,Particle.Coords)
        end
        if Worm.state
            # if QcasiRand(6,Particle.Count) > 0.5
            if rand() > 0.5
                Worm_advance!(Worm::Worm_,Particle::Particle_)
            else
                Worm_recede!(Worm::Worm_,Particle::Particle_)
            end
        else
            Worm_close!(Worm::Worm_,Particle::Particle_)
        end
        if Particle.Type[particle_type]==0 && Worm.state
            Worm_swap!(Worm::Worm_,Particle::Particle_)
        end
    end
end 

function Worm_open!(Worm::Worm_,Particle::Particle_)
    # Worm.atom_i = QcasiRand(7,Particle.Count)
    # Worm.ira = QcasiRand(8,Particle.Count)
    # σ = QcasiRand(9,Particle.Count)
    Worm.atom_i = rand(1:Particle.Number[type])
    Worm.ira = rand(1:Particle.Timeslices)
    σ = rand(1:Worm.m)

    Worm.masha = (Worm.ira + σ) % Worm.Timeslices
    Worm.atom_m = Worm.atom_i

    if Worm.masha != (Worm.ira + σ)
        Worm.atom_m = Particle_Index[Worm.atom_i]
    end
    Probability = Worm_open_p(σ, Worm::Worm_,Particle::Particle_)
    if Probability >= 1 || Probability > QcasiRand(10,Particle.Count)
        Worm.state = true
    end
end

function Worm_close!(Worm::Worm_,Coords)
    σ = Worm.masha - Worm.ira
    if σ <= 0
        σ += Worm.Timeslices
    end
    if σ <= Worm.m
        i_t_0 = Worm.ira
        i_t_2 = Worm.ira + σ 
        Sample_middle!(i_t_0,i_t_2,Worm.atom_i,Worm.atom_m,Coords)
        Probability = 1/Worm_open_p(σ, Worm::Worm_,Particle::Particle_)
        # if Probability >= 1 || Probability > QcasiRand(11,Particle.Count)
        if Probability >= 1 || Probability > rand()
            Worm.state = false
        end
    end
end

function Worm_advance!(Worm::Worm_,Particle::Particle_)
    σ += Worm.masha - Worm.ira
    if σ < 0
        σ += Worm.Timeslices
    end
    # advance = QcasiRand(12,Particle.Count) + 1 
    advance = rand(1:Worm.m) + 1 
    if σ - advance > 0
        type = Worm.type
        offset =  Particle.Type_offset[type]
        i_t_0 = Worm.ira
        i_t_2 = Worm.ira + advance
        ira_new = i_t_2 % Worm.Timeslices
        atom_i_new = Worm.atom_i
        gvar = 1/advance*Particle.T_wave²[type]
        p_t_0 = offset + (Worm.atom_i+i_t_0)*Worm.Timeslices
        p_t_2 = offset + (atom_i_new+i_t_2)*Worm.Timeslices
        for dim in 1:Dimension
            # Particle.Coords[p_t_2+1,dim] = Particle.Coords[p_t_0+1,dim] + QcasiRand(12,Particle.Count)
            Particle.Coords[p_t_2+1,dim] = Particle.Coords[p_t_0+1,dim] +  Gaussian(gvar)###
        end
        Sample_middle!(i_t_0,i_t_2,Worm.atom_i,atom_i_new,Particle)
        E_p = Worm_Potential_Energy(i_t_0,i_t_2+1,Worm.atom_i,atom_i_new, Worm::Worm_,Particle::Particle_)
        if E_p < 0 || exp(-E_p*τ) > QcasiRand(13,Particle.Count)
            Worm.ira = ira_new
            Worm.atom_i = atom_i_new
        end
    end
end

function Worm_recede!(Worm::Worm_,Particle::Particle_)
    σ = Worm.ira - Worm.masha
    if segm < 0
        σ += Worm.Timeslices
    end
    recede = QcasiRand(14,Particle.Count)
    if σ >= recede + 1
        i_t_0 = Worm.ira - recede
        i_t_1 = Worm.ira
        atom_0 = Worm.atom_i
        atom_1 = Worm.atom_i
        if i_t_0 < 0
            i_t_0 += Worm.Timeslices
            i_t_1 += Worm.Timeslices
            atom_0 = Ring_Index[atom_0]
        end
        E_p = Worm_Potential_Energy(i_t_0,i_t_1,atom_0,atom_1, Worm::Worm_,Particle::Particle_)
        # if E_p > 0 || exp(E_p*Particle.τ) > QcasiRand(15,Particle.Count)
        if E_p > 0 || exp(E_p*Particle.τ) > rand()
            Worm.ira = i_t_0 % Worm.Timeslices
            Worm.atom_i = atom_0
        end
    end
end

function Worm_Potential_Energy(i_t_0,i_t_1,atom_0,atom_1, Worm::Worm_,Particle::Particle_)
    atom_offset = Particle.Type_offset[Worm.type]/Worm.Timeslices
    p_i_t_0 = i_t_0 % Worm.Timeslices
    # p_i_t_1 = i_t_1 % Worm.Timeslices
    E_p = 0.0
    atom = atom_0
    for i_t in (i_t_0+1):i_t_1
        p_i_t = i_t % Worm.Timeslices
        if p_i_t != i_t && p_i_t_0 == i_t_0
            atom = atom_1
        end
        E_p += Potential_Energy(atom_offset+atom,p_i_t,Particle::Particle_)
    end
end

function Worm_swap!(Worm::Worm_,Particle::Particle_)
    σ = Worm.m
    i_t_0 = Worm.ira
    i_t_1 = i_t_0 + σ
    p_i_t_0 = i_t_0
    p_i_t_1 = i_t_1 % Worm.Timeslices
    atom_w = Worm.atom_i
    count = Particle_table(atom_w,p_i_t_0,p_i_t_1,σ,i_t_1, Worm::Worm_,Particle::Particle_)
    if count > 0
        atom_1,pnorm_old = atom_Swap(count, Worm::Worm_)
        if atom_1 > 0
            atom_0 = atom_i
            if p_i_t_1 != i_t_1
                atom_0 = Ring_Index[atom_1]
            end
            type = Worm.type
            offset_0 = offset_1 = Particle.Type_offset[type] + atom_0*Worm.Timeslices
            offset_1 = offset_1 = Particle.Type_offset[type] + atom_1*Worm.Timeslices
            offset_w = offset_1 = Particle.Type_offset[type] + atom_w*Worm.Timeslices
            for dim in 1:Dimension
                Particle.Coords_forward[dim,offset_0+p_i_t_0+1] = Particle.Coords[dim,offset_w+p_i_t_0+1]
                Particle.Coords_forward[dim,offset_1+p_i_t_1+1] = Particle.Coords[dim,offset_1+p_i_t_1+1]
            end
            Sample_middle(i_t_0,i_t_1,atom_0,atom_1, Worm::Worm_,Particle::Particle_)
            gatom_0 = offset_0/Worm.Timeslices
            gatom_1 = offset_1/Worm.Timeslices
            E_p = 0.0
            gatom = gatom_0
            for i_t in (i_t_0+1):(i_t_1-1)
                p_i_t = i_t % Worm.Timeslices
                if p_i_t != i_t
                    gatom = gatom_1
                end
                E_p = Potential_Energy(gatom,p_i_t,Particle::Particle_) - Potential_Energy(gatom,p_i_t,Particle::Particle_)
            end
            Probability = exp(-E_p*Particle.τ)
            count = Particle_table(atom_0,p_i_t_0,p_i_t_1,σ,i_t_1, Worm::Worm_,Particle::Particle_)
            pnorm_new = 0.0
            for i_c in 1:count
                pnorm_new += Worm._Particle_table[i_c]
            end
            Probability *= pnorm_old/pnorm_new
            # if Probability > QcasiRand(16,Particle.Count)
            # if Probability > 1 || Probability > QcasiRand(16,Particle.Count)
            if Probability > 1 || Probability > rand()
            for i_t in (i_t_0+1):(i_t_1-1)
                    offset = offset_0
                    for dim in 1:Dimension
                        p_i_t = i_t%Worm.Timeslices
                        if p_i_t != i_t
                            offset = offset_1
                        end
                        Coords[dim,offset+p_i_t] = Coords_forward[dim,offset+p_i_t]
                    end
                end
                for i_t in 1:i_t_0+1
                    for dim in 1:Dimension
                        Coords_forward[dim,offset_0+i_t] = Coords[dim,offset_0+i_t]
                        Particle.Coords[dim,offset_0+i_t] = Particle.Coords[dim,offset_w+i_t]
                        Particle.Coords_forward[dim,offset_w+i_t] = Particle.Coords_forward[dim,offset_0+i_t]
                    end
                end
                Ring_atom_w = Ring_Index[atom_w]
                Ring_atom_0 = Ring_Index[atom_0]
                Particle_Index[Ring_atom_w] = atom_0
                Particle_Index[atom_0] = Ring_atom_w
                Particle_Index[Ring_atom_w] = atom_w
                Particle_Index[atom_w] = Ring_atom_0
                
                if Worm.ira > Worm.masha
                    if atom_0 == Worm.atom_m
                        Worm.atom_m = Worm.atom_i
                    elseif Worm.atom_i == Worm.atom_m
                        Worm.atom_m = atom_0
                    end
                end
            end
        end
    end
end

function Sample_middle(i_t_0,i_t_2,atom_0,atom_2, Worm::Worm_,Particle::Particle_)
    if i_t_2 - i_t_0 >= 2
        i_t_1 = Int((i_t_0+i_t_2) ÷ 2)
        p_t_0 = i_t_0 % Worm.Timeslices
        p_t_1 = i_t_1 % Worm.Timeslices
        p_t_2 = i_t_2 % Worm.Timeslices
        if p_t_1 != i_t_1 && p_t_0 == i_t_0
            atom_1 = atom_2
        else
            atom_1 = atom_0
        end
        offset = Particle.Type_offset[Worm.type]
        p_t_0 += offset + atom_0*Worm.Timeslices
        p_t_1 += offset + atom_1*Worm.Timeslices
        p_t_2 += offset + atom_2*Worm.Timeslices
        s0 = i_t_1 - i_t_0
        s2 = i_t_2 - i_t_1
        gkin = (s0 + s2)/(Particle.T_wave²[Worm.type])
        for dim in 1:Dimension
            Particle.Coords[dim,p_t_1] = (s2*Particle.Coords[dim,p_t_0]+s0*Particle.Coords[dim,p_t_2])/(s2+s2)
            Particle.Coords[dim,p_t_1] += Gaussian(gkin)
        end
        Sample_middle(i_t_0,i_t_1,atom_0,atom_1, Worm::Worm_,Particle::Particle_)
        Sample_middle(i_t_1,i_t_2,atom_1,atom_2, Worm::Worm_,Particle::Particle_)
    end
end

function Particle_table(atom_w,p_t_0,p_t_1,σ,t_1, Worm::Worm_,Particle::Particle_)
    type = Worm.type
    offset = Particle.Type_offset[type]
    i_t_w = offset + atom_w * Worm.Timeslices + p_t_0
    count = 0
    for atom_1 in 0:Particle.Number[type]-1
        if Wordline_check(atom_1,p_t_1)
            if t_1 != p_t_1
                atom_0 = Ring_Index[atom_1]
            else
                atom_0 = atom_1
            end
            if atom_0 != Worm.atom_i
                i_t_1 = offset + atom_1*Worm.Timeslices + p_t_1
                dr² = 0.0
                for dim in 1:Dimension
                    dx = Particle.Coords[dim,i_t_w+1] - Particle.Coords[dim,i_t_1+1]
                    dr² += dx^2
                end
                if dr² < Cutoff²
                    count += 1
                    Worm.neighbor[1,count] = dr²
                    Worm.neighbor[2,count] = atom_1
                end
            end
        end
    end
    sort!(neighbor)
    if count > Max_neighbors
        count = Max_neighbors
    end
    norm = 1/σ*Particle.T_wave²[Worm.type]

    for i_c in 1:count
        Worm._Particle_table[i_c] = exp(-norm*Worm.neighbor[1,i_c])
    end

    return count
end

function atom_Swap(count, Worm::Worm_)
    pnorm = 0.0
    for i_c in 1:count
        pnorm += Worm._Particle_table[i_c]
    end
    prand = pnorm*rand()
    sum = 0.0
    i_c = 0
    while i_c < count && sum < prand
        i_c += 1
        sum += Worm._Particle_table[i_c]
    end
    return Int(Worm.neighbor[2,i_c]),pnorm
end

function Worm_open_p(σ, Worm::Worm_,Particle::Particle_)
    offset = atom_offset[Worm.type]
    kin = 0.0
    p_t_0 = offset + Worm.atom_i*Worm.Timeslices + Worm.ira
    p_t_1 = offset + Worm.atom_m*Worm.Timeslices + Worm.masha
    for dim in 1:Dimension
        dr = Particle.Coords[dim,p_t_0+1] - Particle.Coords[dim,p_t_1+1]
        kin += dr²
    end
    kin /= Particle.T_wave²[Worm.type]*σ
    E_p = Worm_Potential_Energy(Worm.ira,Worm.ira+σ,Worm.atom_i,Worm.atom_m, Worm::Worm_,Particle::Particle_)
    E_p *= Particle.τ

    return Worm.c*Particle.Number[Worm.type]*Worm.Timeslices*Worm.m^(σ,Dimension*0.5)*exp(kin+pot)
end


