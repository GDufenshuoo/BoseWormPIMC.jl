"""
This program is based on Moribs-PIMC
code by Zilong Wang <clockwang@icloud.com>
"""

export 
    Worm_Init,
    Worm_move

function Worm_Init(type::Int,worm_m::Int,worm_c::Float64,Particle::Particle_,System::System_setting)
    Worm = Worm_(
        0,
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
        zeros(Float64, Particle.Number[type]),
        zeros(Float64, Particle.Number[type]),
        zeros(Float64, Particle.Number[type]),
        false
        )

    Worm.Cutoff = 100

    Worm.m = worm_m
    Worm.Timeslices = Particle.Timeslices
    Worm.type = Int(type)
    Worm.c = worm_c*System.Density/(Particle.Number[type]*Worm.Timeslices*Worm.m)
    Worm._norm = Worm.c*Particle.Number[type]*Worm.Timeslices*Worm.m
    Worm.Cutoff² = Worm.Cutoff^2 * Worm.m*4*Particle.λ[type]*Particle.τ

    return Worm
end

function Worm_move!(particle_type::Int, Worm::Worm_=Worm, Particle::Particle_=Particle)
    Worm_pass = zeros(Int,10)
    for i in 1:Particle.Number[particle_type]
        if Worm.state
            Worm_pass[3] += 1
            Worm_pass[4] += Worm_close!(Worm::Worm_,Particle::Particle_)
        else
            Worm_pass[1] += 1
            Worm_pass[2] += Worm_open!(particle_type,Worm::Worm_,Particle::Particle_)
        end
        if Worm.state
            # if QcasiRand(6,Particle.Count) > 0.5
            if rand() > 0.5
                Worm_pass[5] += 1
                Worm_pass[6] += Worm_advance!(Worm::Worm_,Particle::Particle_)
            else
                Worm_pass[7] += 1
                Worm_pass[8] += Worm_recede!(Worm::Worm_,Particle::Particle_)
            end
        else
            Worm_pass[3] += 1
            Worm_pass[4] += Worm_close!(Worm::Worm_,Particle::Particle_)
        end
        if Particle.Type[particle_type]==0 && Worm.state
            Worm_pass[9] += 1
            Worm_pass[10] += Worm_swap!(Worm::Worm_,Particle::Particle_)
        end
    end
    return Worm_pass
end 

function Worm_open!(particle_type::Int,Worm::Worm_,Particle::Particle_)
    # Worm.atom_i = QcasiRand(7,Particle.Count)
    # Worm.ira = QcasiRand(8,Particle.Count)
    # σ = QcasiRand(9,Particle.Count)
    Worm.atom_i = rand(0:Particle.Number[particle_type]-1)
    Worm.ira = rand(0:Worm.Timeslices-1)
    σ = rand(1:Worm.m)

    Worm.masha = (Worm.ira + σ) % Worm.Timeslices
    Worm.atom_m = Worm.atom_i

    if Worm.masha != (Worm.ira + σ)
        Worm.atom_m = Particle.Particle_Index[Worm.atom_i+1]
    end
    Probability = Worm_open_p(σ, Worm::Worm_,Particle::Particle_)
    # if Probability >= 1 || Probability > QcasiRand(10,Particle.Count)
    if Probability >= 1 || Probability > rand()
        Worm.state = true
        return 1
    end
    return 0
end

function Worm_close!(Worm::Worm_,Particle::Particle_)
    σ = Worm.masha - Worm.ira
    if σ < 0
        σ += Worm.Timeslices
    end
    if σ <= Worm.m
        i_t_0 = Worm.ira
        i_t_2 = Worm.ira + σ 
        Sample_middle!(i_t_0,i_t_2,Worm.atom_i,Worm.atom_m,Worm::Worm_,Particle::Particle_)
        Probability = 1/Worm_open_p(σ, Worm::Worm_,Particle::Particle_)
        # if Probability >= 1 || Probability > QcasiRand(11,Particle.Count)
        if Probability >= 1 || Probability > rand()
            Worm.state = false
            return 1
        end
    end
    return 0
end

function Worm_advance!(Worm::Worm_,Particle::Particle_)
    σ = Worm.masha - Worm.ira
    if σ < 0
        σ += Worm.Timeslices
    end
    # advance = QcasiRand(12,Particle.Count) + 1 
    advance = rand(1:Worm.m)
    if σ > advance
        type = Worm.type
        wordline =  Particle.wordline[type]#* Worm.Timeslices
        i_t_0 = Worm.ira
        i_t_2 = Worm.ira + advance
        ira_new = i_t_2 % Worm.Timeslices
        atom_i_new = Worm.atom_i
        gvar = 1/(advance*Particle.T_wave²[type])
        # println(gvar)
        if gvar < 1e-1
            gvar = 1e-1
        end
        p_t_0 = (wordline + Worm.atom_i)*Worm.Timeslices + i_t_0 % Worm.Timeslices+1
        p_t_2 = (wordline + atom_i_new)*Worm.Timeslices + i_t_2 % Worm.Timeslices+1
        for dim in 1:Dimension
            # Particle.Coords[p_t_2+1,dim] = Particle.Coords[p_t_0+1,dim] + QcasiRand(12,Particle.Count)
            Particle.Coords[p_t_2,dim] = Particle.Coords[p_t_0,dim] +  Gaussian(gvar)
        end
        Sample_middle!(i_t_0,i_t_2,Worm.atom_i,atom_i_new, Worm::Worm_,Particle::Particle_)
        E_p = Worm_Potential_Energy(i_t_0,i_t_2,Worm.atom_i,atom_i_new, Worm::Worm_,Particle::Particle_)
        # if E_p < 0 || exp(-E_p*τ) > QcasiRand(13,Particle.Count)
        if E_p < 0 || exp(-E_p*Particle.τ) > rand()
            Worm.ira = ira_new
            Worm.atom_i = atom_i_new
            return 1
        end
    end
    return 0
end

function Worm_recede!(Worm::Worm_,Particle::Particle_)
    σ = Worm.ira - Worm.masha
    if σ < 0
        σ += Worm.Timeslices
    end
    # recede = QcasiRand(14,Particle.Count)
    recede = rand(1:Worm.m)
    if σ >= recede + 1
        i_t_0 = Worm.ira - recede
        i_t_1 = Worm.ira
        atom_0 = Worm.atom_i
        atom_1 = Worm.atom_i
        if i_t_0 < 0
            i_t_0 += Worm.Timeslices
            i_t_1 += Worm.Timeslices
            atom_0 = Particle.Ring_Index[atom_0+1]
        end
        E_p = Worm_Potential_Energy(i_t_0,i_t_1,atom_0,atom_1, Worm::Worm_,Particle::Particle_)
        # if E_p > 0 || exp(E_p*Particle.τ) > QcasiRand(15,Particle.Count)
        if E_p > 0 || exp(E_p*Particle.τ) > rand()
            Worm.ira = i_t_0 % Worm.Timeslices
            Worm.atom_i = atom_0
            return 1
        end
    end
    return 0
end

function Worm_Potential_Energy(i_t_0::Int,i_t_1::Int,atom_0::Int,atom_1::Int, Worm::Worm_,Particle::Particle_)
    atom_offset = Particle.wordline[Worm.type]
    p_i_t_0 = i_t_0 % Worm.Timeslices
    # p_i_t_1 = i_t_1 % Worm.Timeslices
    E_p = 0.0
    atom = atom_0
    for i_t in i_t_0+1:(i_t_1-1)
        p_i_t = i_t % Worm.Timeslices
        if p_i_t != i_t && p_i_t_0 == i_t_0
            atom = atom_1
        end
        E_p += Potential_Energy_single(atom_offset+atom,Particle.Coords,p_i_t)
    end
    E_p
end

function Worm_swap!(Worm::Worm_,Particle::Particle_)
    σ = Worm.m
    i_t_0 = Worm.ira
    i_t_1 = i_t_0 + σ
    p_i_t_0 = i_t_0
    p_i_t_1 = i_t_1 % Worm.Timeslices
    atom_w = Worm.atom_i
    Count_index = Particle_table(atom_w,p_i_t_0,p_i_t_1,σ,i_t_1, Worm::Worm_,Particle::Particle_)
    if Count_index > 0
        atom_1,pnorm_old = atom_Swap(Count_index, Worm::Worm_)
        if atom_1 > 0
            atom_0 = atom_1
            if p_i_t_1 != i_t_1
                atom_0 = Particle.Ring_Index[atom_1+1]
            end
            type = Worm.type
            offset_0 = (Particle.wordline[type] + atom_0)*Worm.Timeslices
            offset_1 = (Particle.wordline[type] + atom_1)*Worm.Timeslices
            offset_w = (Particle.wordline[type] + atom_w)*Worm.Timeslices
            for dim in 1:Dimension
                Particle.Coords_forward[offset_0+p_i_t_0+1,dim] = Particle.Coords[offset_w+p_i_t_0+1,dim]
                Particle.Coords_forward[offset_1+p_i_t_1+1,dim] = Particle.Coords[offset_1+p_i_t_1+1,dim]
            end
            Sample_middle!(i_t_0,i_t_1,atom_0,atom_1, Worm::Worm_,Particle::Particle_)
            gatom_0 = div(offset_0,Worm.Timeslices)
            gatom_1 = div(offset_1,Worm.Timeslices)
            E_p = 0.0
            gatom = gatom_0
            for i_t in (i_t_0+1):(i_t_1-1)
                p_i_t = Int(i_t % Worm.Timeslices)
                if p_i_t != i_t
                    gatom = gatom_1
                end
                E_p = Potential_Energy_single(gatom,Particle.Coords,p_i_t) - Potential_Energy_single(gatom,Particle.Coords,p_i_t)
            end
            Probability = exp(-E_p*Particle.τ)
            Count_index = Particle_table(atom_0,p_i_t_0,p_i_t_1,σ,i_t_1, Worm::Worm_,Particle::Particle_)
            pnorm_new = 0.0
            for i_c in 1:Count_index
                pnorm_new += Worm.Particle_table[i_c]
            end
            Probability *= pnorm_old/pnorm_new
            # if Probability > QcasiRand(16,Particle.Count)
            # if Probability > 1 || Probability > QcasiRand(16,Particle.Count)
            if Probability >= 1 || Probability > rand()
                for i_t in (i_t_0+1):(i_t_1-1)
                    offset = offset_0
                    for dim in 1:Dimension
                        p_i_t = i_t%Worm.Timeslices
                        if p_i_t != i_t
                            offset = offset_1
                        end
                        Particle.Coords[offset+p_i_t+1,dim] = Particle.Coords_forward[offset+p_i_t+1,dim]
                    end
                end
                for i_t in 1:i_t_0+1
                    for dim in 1:Dimension
                        Particle.Coords_forward[offset_0+i_t,dim] = Particle.Coords[offset_0+i_t,dim]
                        Particle.Coords[offset_0+i_t,dim] = Particle.Coords[offset_w+i_t,dim]
                        Particle.Coords_forward[offset_w+i_t,dim] = Particle.Coords_forward[offset_0+i_t,dim]
                    end
                end
                Ring_atom_w = Particle.Ring_Index[atom_w+1]
                Ring_atom_0 = Particle.Ring_Index[atom_0+1]
                Particle.Particle_Index[Ring_atom_w+1] = atom_0
                Particle.Ring_Index[atom_0+1] = Ring_atom_w
                Particle.Particle_Index[Ring_atom_w+1] = atom_w
                Particle.Ring_Index[atom_w+1] = Ring_atom_0
                
                if Worm.ira > Worm.masha
                    if atom_0 == Worm.atom_m
                        Worm.atom_m = Worm.atom_i
                    elseif Worm.atom_i == Worm.atom_m
                        Worm.atom_m = atom_0
                    end
                end
                return 1
            end
        end
    end
    return 0
end

function Sample_middle!(i_t_0::Int,i_t_2::Int,atom_0::Int,atom_2::Int, Worm::Worm_=Worm,Particle::Particle_=Particle)
    if i_t_2 >= i_t_0 + 2
        i_t_1 = div(i_t_0+i_t_2,2)
        p_t_0 = i_t_0 % Worm.Timeslices
        p_t_1 = i_t_1 % Worm.Timeslices
        p_t_2 = i_t_2 % Worm.Timeslices
        if p_t_1 != i_t_1 && p_t_0 == i_t_0
            atom_1 = atom_2
        else
            atom_1 = atom_0
        end
        offset = Particle.wordline[Worm.type]* Worm.Timeslices
        p_t_0 += offset + atom_0*Worm.Timeslices +1
        p_t_1 += offset + atom_1*Worm.Timeslices +1
        p_t_2 += offset + atom_2*Worm.Timeslices +1
        s0 = i_t_1 - i_t_0
        s2 = i_t_2 - i_t_1
        gkin = (s0 + s2)/(Particle.T_wave²[Worm.type])
        for dim in 1:Dimension
            Particle.Coords[p_t_1,dim] = (s2*Particle.Coords[p_t_0,dim]+s0*Particle.Coords[p_t_2,dim])/(s0+s2)
            Particle.Coords[p_t_1,dim] += Gaussian(gkin)
        end
        Sample_middle!(i_t_0,i_t_1,atom_0,atom_1, Worm::Worm_,Particle::Particle_)
        Sample_middle!(i_t_1,i_t_2,atom_1,atom_2, Worm::Worm_,Particle::Particle_)
    end
end

function Particle_table(atom_w::Int,p_t_0::Int,p_t_1::Int,σ::Int,t_1::Int, Worm::Worm_,Particle::Particle_)
    type = Worm.type
    offset = Particle.wordline[type]*Worm.Timeslices
    i_t_w = offset + atom_w * Worm.Timeslices + p_t_0
    Count_index = 0
    for atom_1 in 0:Particle.Number[type]-1
        if Wordline_check(atom_1,p_t_1)
            if t_1 != p_t_1
                atom_0 = Particle.Ring_Index[atom_1+1]
            else
                atom_0 = atom_1
            end
            if atom_0 != Worm.atom_i
                i_t_1 = offset + atom_1*Worm.Timeslices + p_t_1
                dr² = 0.0
                for dim in 1:Dimension
                    dr² += (Particle.Coords[i_t_w+1,dim] - Particle.Coords[i_t_1+1,dim])^2
                end
                if dr² < Worm.Cutoff
                    Count_index += 1
                    Worm.dr²[Count_index] = dr²
                    Worm.atom_list[Count_index] = atom_1
                end
            end
        end
    end

    if Count_index > Worm.Cutoff
        Count_index = Worm.Cutoff
    end
    Count_index = ceil(Int,Count_index)
    MM_Sort!(Worm.dr²,Worm.atom_list,Count_index)
    Worm.Particle_table .= exp.(Worm.dr²./(-σ* Particle.T_wave²[type]))

    return Count_index
end

function atom_Swap(Count_index::Int, Worm::Worm_)
    pnorm = 0.0
    for i_c in 1:Count_index
        pnorm += Worm.Particle_table[i_c]
    end
    prand = pnorm*rand()
    sum = 0.0
    i_c = 1
    while i_c <= Count_index && sum <= prand
        sum += Worm.Particle_table[i_c]
        i_c += 1
    end
    return Int(Worm.atom_list[i_c-1]),pnorm
end

function Worm_open_p(σ::Int, Worm::Worm_,Particle::Particle_)
    offset = Particle.wordline[Worm.type]*Worm.Timeslices
    kin = 0.0
    p_t_0 = offset + Worm.atom_i*Worm.Timeslices + Worm.ira
    p_t_1 = offset + Worm.atom_m*Worm.Timeslices + Worm.masha
    for dim in 1:Dimension
        dr = Particle.Coords[p_t_0+1,dim] - Particle.Coords[p_t_1+1,dim]
        kin += dr^2
    end
    kin /= Particle.T_wave²[Worm.type]*σ
    E_p = Worm_Potential_Energy(Worm.ira,Worm.ira+σ,Worm.atom_i,Worm.atom_m, Worm::Worm_,Particle::Particle_)

    return Worm._norm*σ^(Dimension*0.5)*exp(kin+E_p*Particle.τ)
end

function Wordline_check(atom::Int,p_t::Int)
    wordline = true
    if atom == Worm.atom_m || atom == Worm.atom_i
        if Worm.atom_i != Worm.atom_m || Worm.ira > Worm.masha
            if (atom == Worm.atom_m && p_t < Worm.masha) || (atom == Worm.atom_i && p_t > Worm.ira)
                wordline = false                
            end
        elseif p_t > Worm.ira && p_t < Worm.masha
            wordline = false
        end
    end
    return wordline
end


"""
Some Tools
"""
function MM_Sort!(dist,labels,Count_index::Int)
    for j in 2:Count_index
        dr_temp = dist[j]
        i_temp = labels[j]
        i = j - 1
        while i>0 && dist[i]>dr_temp
            dist[i+1] = dist[i]
            labels[i+1] = labels[i]
            i -= 1
        end
        dist[i+1] = dr_temp
        labels[i+1] = i_temp
    end
end


