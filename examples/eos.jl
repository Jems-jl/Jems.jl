using Jems 
using ForwardDiff
using DelimitedFiles


struct Tρ_table_EOS <: Jems.EOS.AbstractEOS
    X :: Float64 #redundant 
    Z :: Float64
    logTs :: Vector{Float64}
    logQs :: Vector{Float64}
    eos_data :: Array{Float64, 3}

    col_names :: Vector{String}
end 

struct EOS_Z_slice
    Z :: Float64
    Xs :: Vector{Float64}
    tables :: Vector{Tρ_table_EOS}
end
struct EOS_table_collector <: Jems.EOS.AbstractEOS
    Zs :: Vector{Float64}
    slices :: Vector{EOS_Z_slice}
    include_radiation :: Bool
end 


function Tρ_table_EOS(filepath :: String)
    lines = readlines(filepath)
    header_idx = findfirst(l -> occursin("version", l), lines)

# The actual metadata values are on the line directly after the header
    meta_idx = header_idx + 1
    parts = split(lines[meta_idx])

    X_val = parse(Float64, parts[2])
    Z_val = parse(Float64, parts[3])
    # Extract Grid Dimensions
    num_Ts   = parse(Int, parts[4])
    logT_min = parse(Float64, parts[5])
    # logT_max = parse(Float64, parts[6])
    logT_del = parse(Float64, parts[7])
    
    num_Qs    = parse(Int, parts[8])
    logQ_min  = parse(Float64, parts[9])
    # logQ_max  = parse(Float64, parts[10])
    logQ_del  = parse(Float64, parts[11])

    logTs = collect(range(logT_min, step = logT_del, length = num_Ts))
    logQs= collect(range(logQ_min, step = logQ_del, length = num_Qs))
    
    # Finding the header row // can be changed to a number later but this is a safety check working for MESA EOS data 
    header_row_idx = findfirst(l -> occursin("logE", l), lines)
    col_names = String.(split(lines[header_row_idx]))
    num_vars = length(col_names)

    # Allocating 3D array 
    eos_data_3D = fill(NaN, (num_vars, num_Ts, num_Qs))

    gap_size = 6

    for iQ in 1:num_Qs
        for iT in 1:num_Ts
            current_line_idx = header_row_idx + (iQ-1)*(num_Ts + gap_size) + iT 
            # Safety check
            if current_line_idx > length(lines)
                error("Unexpected end of file at Q index $iQ, T index $iT")
            end

            line = lines[current_line_idx]
            row_values = parse.(Float64, split(line))
            eos_data_3D[:, iT, iQ] .= row_values

        end
    end 
    return Tρ_table_EOS(X_val, Z_val, logTs, logQs, eos_data_3D, col_names)
    
end 


function EOS_table_collector(directory :: String)
    files = readdir(directory; join = true)
    candidate_files = filter(f -> endswith(f, ".data"), files)
    #Load all tables
    all_tables = [Tρ_table_EOS(f) for f in candidate_files]

    unique_Zs = sort(unique([t.Z for t in all_tables]))

    slices = EOS_Z_slice[]

    for z_val in unique_Zs
        tables_at_z = filter(t -> t.Z == z_val, all_tables)

        #Sort by X in each Z slice
        sort!(tables_at_z, by = t -> t.X)
        #X values in each Z slice 
        xs_at_z = [t.X for t in tables_at_z]

        push!(slices, EOS_Z_slice(z_val, xs_at_z, tables_at_z))
    end
    println("EOS tables loaded")
    return EOS_table_collector(unique_Zs, slices, true)
end 

@inline function cubic_weights(t :: T) where T
    t2 = t * t
    t3 = t2 * t 

    # Standard Catmull-Rom weights: [1 t t^2 t^3] M_SH (M_SH is the  Catmull-Rom Cubic Basis matrix for 1D)
    # https://en.wikipedia.org/wiki/Catmull–Rom_spline
    w_m1 = -0.5 * t3 + 1.0 * t2 - 0.5 * t
    w_0  =  1.5 * t3 - 2.5 * t2 + 1.0
    w_1  = -1.5 * t3 + 2.0 * t2 + 0.5 * t
    w_2  =  0.5 * t3 - 0.5 * t2
    
    return w_m1, w_0, w_1, w_2
end

"""
The function takes the parameter for which we want to get the 1D weights for the bicubic interpolation 
and returns the weights, and the closest index in the data from the binary search 
"""
@inline function get_cubic_weights(grid:: Vector{Float64}, val :: T) where T 
    N = length(grid)
    if N == 1
        return 1, (zero(T), one(T), zero(T), zero(T))
    end
    val_safe = clamp(val, grid[1], grid[end])
    i = searchsortedlast(grid, ForwardDiff.value(val))
    i = clamp(i, 2, max(2, N-2)) #for small neighbours less than 4 
    denom = grid[i+1] - grid[i]
    u = (denom == 0) ? zero(T) : (val_safe - grid[i]) / denom

    return i, cubic_weights(u) # Returns tuple (wm1, w0, w1, w2)
end 

# function get_eos_data_per_table(table::Tρ_table_EOS, col_idx :: Int, val_logT :: Float64, val_logQ :: Float64, logT :: T, logQ :: T ) where T 

#     #Table boundary 
#     min_T, max_T = table.logTs[1], table.logTs[end]
#     min_Q, max_Q = table.logQs[1], table.logQs[end]

#     eff_logT = clamp(logT, min_T, max_T)
#     eff_logQ = clamp(logQ, min_Q, max_Q)

#     # Also clamp the Float values for index searching
#     eff_val_logT = clamp(val_logT, table.logTs[1], table.logTs[end])
#     eff_val_logQ = clamp(val_logQ, table.logQs[1], table.logQs[end])

#     i_T = searchsortedlast(table.logTs, eff_val_logT)
#     i_Q = searchsortedlast(table.logQs, eff_val_logQ)
    
#     #length of the tables
#     N_T = length(table.logTs)
#     N_Q = length(table.logQs)

#     #clamping the indexes ()
#     i_T = clamp(i_T, 1, N_T - 1)
#     i_Q = clamp(i_Q, 1, N_Q - 1)

#     T_min, T_max = table.logTs[i_T], table.logTs[i_T+1]
#     Q_min, Q_max = table.logQs[i_Q], table.logQs[i_Q+1]

#     u = (eff_logT - T_min) / (T_max - T_min)
#     v = (eff_logQ - Q_min) / (Q_max - Q_min)

#     wt_m1, wt_0, wt_1, wt_2 = cubic_weights(u)
#     wr_m1, wr_0, wr_1, wr_2 = cubic_weights(v)

#     weights_T = (wt_m1, wt_0, wt_1, wt_2)
#     weights_Q = (wr_m1, wr_0, wr_1, wr_2)

#     result = zero(T)

#     for (idx_T, offset_T) in enumerate(-1:2)
#         for (idx_Q, offset_Q) in enumerate(-1:2)
#             curr_T_idx = clamp(i_T + offset_T, 1, N_T)
#             curr_Q_idx = clamp(i_Q + offset_Q, 1, N_Q)

#             val = table.eos_data[col_idx, curr_T_idx, curr_Q_idx] # Accessing the value for that particular iT and iQ index, and col_idx provides the variable required

#             result += val * weights_T[idx_T] * weights_Q[idx_Q]
#         end 
#     end 
#     return  result
# end 

##

# function get_eos_data(collection :: EOS_table_collector, lnT::TT, lnRho::TT, xa::AbstractVector{<:TT}, species::Vector{Symbol}) where {TT<:Real}
#     inv_ln10 = 0.4342944819032518
#     logT = lnT * inv_ln10
#     logRho = lnRho * inv_ln10
#     logQ = logRho - 2*logT + 12

#     val_logT = ForwardDiff.value(logT)
#     val_logQ = ForwardDiff.value(logQ)

#     # iH1 = findfirst(==(:H1), species)
#     # X_dual = xa[iH1]
#     # iHe4 = findfirst(==(:He4), species)
#     # Y_dual = xa[iHe4]
#     # Z_dual = 1 - X_dual -Y_dual

    
#     iH1  = findfirst(==(:H1), species)
#     iH2  = findfirst(==(:H2), species)  # Might be nothing
#     iHe3 = findfirst(==(:He3), species) # Common in PPI chain
#     iHe4 = findfirst(==(:He4), species)

#     val_H1  = xa[iH1]
#     val_H2  = isnothing(iH2)  ? 0.0 : xa[iH2]
#     val_He3 = isnothing(iHe3) ? 0.0 : xa[iHe3]
#     val_He4 = xa[iHe4]

#     X_dual = val_H1 + val_H2
#     Y_dual = val_He3 + val_He4
#     Z_dual = 1.0 - X_dual - Y_dual
        
#     val_X = ForwardDiff.value(X_dual)
#     val_Z = ForwardDiff.value(Z_dual)

#     Nz = length(collection.Zs)

#     iz = clamp(searchsortedlast(collection.Zs, val_Z), 1, Nz - 1)

#     u_Z = (Z_dual - collection.Zs[iz]) / (collection.Zs[iz+1] - collection.Zs[iz])

    
#     wZ_m1, wZ_0, wZ_1, wZ_2 = cubic_weights(u_Z)
#     weights_Z = (wZ_m1, wZ_0, wZ_1, wZ_2)

#     interp_vals = ntuple(Val(18)) do k
#         col_idx = k + 1  
        
#         col_val_accum = zero(TT)

#         for (idx_Z, offset_Z) in enumerate(-1:2)
#             real_iz = clamp(iz + offset_Z, 1, Nz) #getting the correct Z and clamping if required 
#             slice = collection.slices[real_iz]

#             Nx = length(slice.Xs)

#             if Nx > 1 
#                 ix = searchsortedlast(slice.Xs, val_X)
#                 ix = clamp(ix, 1, Nx - 1)

#                 X_min, X_max = slice.Xs[ix], slice.Xs[ix+1]
#                 u_X = (X_dual - X_min) / (X_max - X_min)

#             else
#                 # Edge case: Slice has only 1 table (e.g., Z=1.0, X=0.0)
#                 ix = 1
#                 u_X = zero(TT)
#             end 

#             wX_m1, wX_0, wX_1, wX_2 = cubic_weights(u_X)
#             weights_X_local = (wX_m1, wX_0, wX_1, wX_2)

#             for (idx_X, offset_X) in enumerate(-1:2)
                
#                 real_ix = clamp(ix + offset_X, 1, Nx)
#                 table = slice.tables[real_ix]
                
#                 val_TQ = get_eos_data_per_table(table, col_idx, val_logT, val_logQ, logT, logQ)
                
#                 col_val_accum += val_TQ * weights_X_local[idx_X] * weights_Z[idx_Z]
#             end
#         end
#         return col_val_accum

# end
#     return (lnT, lnRho, interp_vals...)
# end

function get_eos_data(c::EOS_table_collector, lnT::TT, lnRho::TT, xa::AbstractVector{TT}, species::Vector{Symbol}) where {TT<:Real}
 
    inv_ln10 = 0.4342944819032518
    logT = lnT * inv_ln10
    logRho = lnRho * inv_ln10
    logQ = logRho - 2*logT + 12

   
    iH1  = findfirst(==(:H1), species)
    iH2  = findfirst(==(:H2), species)
    iHe3 = findfirst(==(:He3), species)
    iHe4 = findfirst(==(:He4), species)
    
    val_H1  = isnothing(iH1) ? 0.0 : xa[iH1]
    val_H2  = isnothing(iH2)  ? 0.0 : xa[iH2]
    val_He3 = isnothing(iHe3) ? 0.0 : xa[iHe3]
    val_He4 = isnothing(iHe4) ? 0.0 : xa[iHe4]
    
    X_val = val_H1 + val_H2
    Y_val = val_He3 + val_He4
    Z_val = 1.0 - X_val - Y_val


    iz, wZ = get_cubic_weights(c.Zs, Z_val)
    Nz = length(c.Zs)

    results = zeros(TT, 18)

    for (dz, weight_Z) in enumerate(wZ)
        real_iz = clamp(iz + dz - 2, 1, Nz)
        slice = c.slices[real_iz]

        ix, wX = get_cubic_weights(slice.Xs, X_val)
        Nx = length(slice.Xs)
        for (dx, weight_X) in enumerate(wX)
            real_ix = clamp(ix + dx - 2, 1, Nx)
            table = slice.tables[real_ix]

            #Calculating the weights for each (T,Q) table once 
            iT, wT = get_cubic_weights(table.logTs, logT)
            iQ, wQ = get_cubic_weights(table.logQs, logQ)
            Nt = length(table.logTs)
            Nq = length(table.logQs)

            w_Table = weight_Z * weight_X

            # Loop T Neighbors
            for (dt, weight_T) in enumerate(wT)
                real_iT = clamp(iT + dt - 2, 1, Nt)
                
                # Loop Q Neighbors
                for (dq, weight_Q) in enumerate(wQ)
                    real_iQ = clamp(iQ + dq - 2, 1, Nq)
                    
                    # Final combined weight
                    w_Point = w_Table * weight_T * weight_Q

                    for k in 1:18
                        val = table.eos_data[k, real_iT, real_iQ]
                        results[k] += val * w_Point
                    end
                end
            end
        end
    end

    return results
end

@inline function smooth_step_func(x::T, floor::Float64, ceil::Float64) where T
    if x <= floor 
        return zero(T)
    elseif x >= ceil
        return one(T)
    else
        t = (x - floor)/(ceil - floor)

        return t * t * (3.0 -2.0 * t)
    end 

end


function Jems.EOS.set_EOS_resultsTρ!(eos:: EOS_table_collector, r::EOSResults{TT}, lnT::TT, lnρ::TT,
                           xa::AbstractVector{TT}, species::Vector{Symbol}) where {TT<:Real}


    logT = log10(exp(lnT)) 
    # logρ = log10(exp(lnρ))    
    # logQ = logρ - 2*logT + 12 
    
    r.lnT = lnT
    r.lnρ = lnρ
    r.T = exp(lnT)
    r.ρ = exp(lnρ)
    
    # calculating weight of the step based on the logT value 
    # Transition 1 at low T (3.1 -> 3.2)
    # Transition 2 at high T (7.0 -> 7.2)

    w_low = smooth_step_func(logT, 3.1, 3.2)
    w_high = 1.0 - smooth_step_func(logT, 7.0, 7.4)
    # Combined Weight: 
    w = w_low * w_high 
    iw = 1.0 - w

    if w >= 0.9999
        eos_data = get_eos_data(eos, lnT, lnρ, xa, species)
        r.μ = eos_data[12]
        logP_gas = eos_data[2]
        P_gas = 10^(logP_gas)

        if eos.include_radiation
            r.Prad = Jems.Constants.CRAD * r.T^4 / 3
        else
            r.Prad = 0
        end
        r.P  = P_gas + r.Prad
        r.lnP = log(r.P)
        r.χ_ρ = eos_data[5]
        r.χ_T = eos_data[6]
        r.β = 1 - r.Prad/r.P
        r.α = 1/eos_data[5]
        r.δ = eos_data[6]/eos_data[5]
        r.u = 10^(eos_data[3])
        r.cₚ = eos_data[7]
        r.∇ₐ = eos_data[16]
        r.Γ₁ = eos_data[14]

    elseif w<=0.0001

        r.μ = Jems.EOS.get_μ_IdealEOS(xa, species)
        if eos.include_radiation
            r.Prad = CRAD * r.T^4 / 3
        else
            r.Prad = 0
        end
        r.P = CGAS * r.T * r.ρ / r.μ + r.Prad
        r.lnP = log(r.P)
        r.β = 1 - r.Prad / r.P  # gas pressure fraction
        r.α = 1 / r.β
        r.δ = (4 - 3 * r.β) / r.β
        r.χ_ρ = r.β
        r.χ_T = 4 - 3 * r.β
        r.u = CGAS * r.T / r.μ * (3 / 2 + 3 * (1 - r.β) / r.β)  # specific internal energy
        r.cₚ = CGAS / r.μ * (3 / 2 + 3 * (4 + r.β) * (1 - r.β) * r.α^2 + r.δ * r.α)  # specific heat capacity at constant P
        r.∇ₐ = CGAS * r.δ / (r.β * r.μ * r.cₚ)  # adiabatic temperature gradient
        r.Γ₁ = 1 / (r.α - r.δ * r.∇ₐ)  # first adiabatic exponent dlnP/dlnρ

    else
        # table data 
        eos_data = get_eos_data(eos, lnT, lnρ, xa, species)
        t_μ     = eos_data[12]
        t_P_gas = 10^(eos_data[2])
        t_u     = 10^(eos_data[3])
        t_χ_ρ   = eos_data[5]
        t_χ_T   = eos_data[6]
        t_cₚ    = eos_data[7]
        t_Γ₁    = eos_data[14]
        t_∇ₐ    = eos_data[16]
        # We do NOT need t_α or t_δ here, we derive them later.

        # Ideal Data
        i_μ = Jems.EOS.get_μ_IdealEOS(xa, species)
        i_P_gas = CGAS * r.T * r.ρ / i_μ
        
        if eos.include_radiation
            r.Prad = CRAD * r.T^4 / 3.0
        else
            r.Prad = 0.0
        end

        i_P_tot = i_P_gas + r.Prad
        i_β = i_P_gas / i_P_tot
        i_α = 1.0 / i_β
        i_δ = (4.0 - 3.0 * i_β) / i_β  
        
        i_χ_ρ = i_β
        i_χ_T = 4.0 - 3.0 * i_β
        

        i_u  = CGAS * r.T / i_μ * (1.5 + 3.0 * (1.0 - i_β) / i_β)
        i_cₚ = CGAS / i_μ * (1.5 + 3.0 * (4.0 + i_β) * (1.0 - i_β) * i_α^2 + i_δ * i_α)
        i_∇ₐ = CGAS * i_δ / (i_β * i_μ * i_cₚ)
        i_Γ₁ = 1.0 / (i_α - i_δ * i_∇ₐ)
        
        # Blend
        r.μ   = w * t_μ     + iw * i_μ
        P_gas = w * t_P_gas + iw * i_P_gas # Blended Gas Pressure
        r.u   = w * t_u     + iw * i_u
        r.χ_ρ = w * t_χ_ρ   + iw * i_χ_ρ
        r.χ_T = w * t_χ_T   + iw * i_χ_T
        r.cₚ  = w * t_cₚ    + iw * i_cₚ
        r.∇ₐ  = w * t_∇ₐ    + iw * i_∇ₐ
        r.Γ₁  = w * t_Γ₁    + iw * i_Γ₁
        
        r.P   = P_gas + r.Prad      
        r.lnP = log(r.P)
        r.β   = 1.0 - r.Prad / r.P
        r.α = 1.0 / r.χ_ρ
        r.δ = r.χ_T / r.χ_ρ
    end
   
end  