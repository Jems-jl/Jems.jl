using Jems 
using ForwardDiff
using DelimitedFiles

struct RT_table_opacity <: Jems.Opacity.AbstractOpacity 
    X :: Float64
    Z :: Float64
    logTs :: Vector{Float64}
    logRs :: Vector{Float64}
    kap_data :: Vector{Float64} #Instead of Array{Float64} 
end 

    
function RT_table_opacity(filepath :: String)
    
    #Reading the whole file into a string array 
    lines = readlines(filepath) 

    meta_idx = findfirst(l -> occursin(r"^\s*1\s+\d+", l), lines) # Looks for 1. for the start of metadata
    if meta_idx === nothing 
        error("Could not find Metadata line (starting with ' 1 '). File format unknown.")
    end

    parts = split(lines[meta_idx]) # Split the string into structured data and removing the spaces

    X_val = parse(Float64, parts[3])
    Z_val = parse(Float64, parts[4])
    num_Rs = parse(Int, parts[5])
    num_Ts = parse(Int, parts[8])

    grid_line_idx = 0 
    logRs = Float64[]

    for i in (meta_idx + 1):min(meta_idx + 20, length(lines))
        row_vals = split(lines[i])
        # Check if line is all numbers and has correct length
        if length(row_vals) == num_Rs && all(x -> tryparse(Float64, x) !== nothing, row_vals)
            grid_line_idx = i
            logRs = parse.(Float64, row_vals)
            break
        end
    end

    if grid_line_idx == 0
        error("Could not find Grid Line! Searched 20 lines after metadata for a row with $num_Rs numbers.")
    end
    println("Grid (logRs) found at Line $grid_line_idx.")

    data_start_idx = 0
    
    for i in (grid_line_idx + 1):min(grid_line_idx + 10, length(lines))
        row_vals = split(lines[i])
        if length(row_vals) == (num_Rs + 1) && all(x -> tryparse(Float64, x) !== nothing, row_vals)
            data_start_idx = i
            break
        end
    end

    if data_start_idx == 0
        error("Could not find Data Start! Searched lines after grid for row with $(num_Rs + 1) numbers.")
    end
    println("Data starts at Line $data_start_idx.")

    ## pre-allocating memory 

    logTs = Vector{Float64}(undef, num_Ts)
    data = Vector{Float64}(undef, num_Rs * num_Ts)


    for T_idx in 1:num_Ts
        current_line_idx = data_start_idx + T_idx - 1 
        if current_line_idx > length(lines)
            error("Unexpected End of File at line $line_idx")
        end
        # Format: "3.750  -0.6258  -0.7084 ..."
        row_parts = split(lines[current_line_idx])
        #Safety check
        if length(row_parts) != (num_Rs + 1)
            error("Format Error at line $current_line_idx: Expected $(num_Rs+1) columns, found $(length(row_parts))")
        end
        logTs[T_idx] = parse(Float64, row_parts[1]) #Storing the Temperature

        for R_idx in 1:num_Rs
            val = parse(Float64, row_parts[R_idx + 1])
            flat_idx = R_idx + (T_idx - 1) * num_Rs
            data[flat_idx] = val 

        end
    end 
    return RT_table_opacity(X_val, Z_val, logTs, logRs, data)
end 


function Jems.Opacity.get_opacity_resultsTρ(table :: RT_table_opacity, lnT::TT, lnρ::TT, xa::AbstractVector{<:TT}, species::Vector{Symbol})::TT where {TT<:Real}
    
    inv_ln10 = 0.4342944819032518
    logT = lnT * inv_ln10
    logρ = lnρ * inv_ln10

    #Calculating logR 
    logR = logρ - 3*logT + 18

    #Converting the dual lnT and lnρ to values so that we can search up 
    val_logT = ForwardDiff.value(logT)
    val_logR = ForwardDiff.value(logR)
    """
    The following part was for clamping the value of T,R at the edges so that value of κ outside the boundary is constant and stable, but this constant value 
    causes issues with duals, so now the value is extrapolated smoothly outside the boundary using a very small slope
    """
    # # Clamping data 
    # min_T, max_T = table.logTs[1], table.logTs[end]
    # min_R, max_R = table.logRs[1], table.logRs[end]

    # eff_logT = clamp(logT, min_T, max_T)
    # eff_logR = clamp(logR, min_R, max_R)

    T_idx = searchsortedlast(table.logTs, val_logT)
    R_idx = searchsortedlast(table.logRs, val_logR)

    ## Boundary cases (For now clamping to the edge)
    T_idx = clamp(T_idx, 1, length(table.logTs)-1)
    R_idx = clamp(R_idx, 1, length(table.logRs)-1)

    #logR and logT values
    R0 = table.logRs[R_idx]
    R1 = table.logRs[R_idx+1]
    T0 = table.logTs[T_idx]
    T1 = table.logTs[T_idx+1]

    # Calculating the slopes 
    u = (logR - R0)/ (R1 - R0)
    v = (logT - T0)/ (T1 - T0)

    function damp_slope(x::T) where T
        if x > 1.0
            # Start at edge (1.0) + tiny fraction of the excess
            return 1.0 + 0.05 * (x - 1.0)
        elseif x < 0.0
            # Start at edge (0.0) + tiny fraction of the deficit
            return 0.0 + 0.05 * x
        else
            return x
        end
    end

    u_damped = damp_slope(u)
    v_damped = damp_slope(v)

    # Now we will look into the 1D array since we have the indexes

    N_R = length(table.logRs)

    # calculating the kappa values 
    #Bottom left value for interpolation
    idx_00 = R_idx + (T_idx - 1) * N_R
    val_00 = table.kap_data[idx_00]
    
    #Bottom right value for interpolation
    idx_10 = (R_idx + 1)+ (T_idx - 1) * N_R
    val_10 = table.kap_data[idx_10]

    #Top left value for interpolation
    idx_01 = R_idx + T_idx * N_R
    val_01 = table.kap_data[idx_01]

    #Top right value for interpolation
    idx_11 = (R_idx + 1) + T_idx * N_R
    val_11 = table.kap_data[idx_11]


    #Bilinear Interpolation
    # f(x,y) = (1-u)(1-v)f00 + u(1-v)f10 + (1-u)v f01 + uv f11
    w00 = (1 - u_damped) * (1 - v_damped)
    w10 = u_damped * (1 - v_damped)
    w01 = (1 - u_damped) * v_damped
    w11 = u_damped * v_damped

    interpolated_log_kappa = (w00 * val_00) + (w10 * val_10) + 
                             (w01 * val_01) + (w11 * val_11)

    return 10^interpolated_log_kappa

end