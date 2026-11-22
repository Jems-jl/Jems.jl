"""
    remesher!

Acts on sm.start_step_props using info from sm.prv_step_props to decide whether to merge/split cells
"""
function remesher!(sm::StellarModel)
    psp = sm.prv_step_props  # for brevity
    ssp = sm.start_step_props
    # we first do cell splitting
    do_split = Vector{Bool}(undef, psp.nz)
    do_split .= false
    Threads.@threads for i in 1:psp.nz-1
        a = abs(log10(get_value(psp.eos_res[i].P)) -
                log10(get_value(psp.eos_res[i+1].P)))
        b = sm.opt.remesh.delta_log10P_split
        if a > b
            # if the condition is satisfied, we split the largest of the two cells
            if psp.dm[i] > psp.dm[i+1]
                do_split[i] = true
            else
                do_split[i+1] = true
            end
        end
    end
    extra_cells = sum(do_split)

    # next we determine cells to merge
    # do_merge[i] == true implies cell i is to be merged with cell i+1
    # vector is created with enough size to account for split cells, as we will be doing cell splitting
    # first and moving the indeces for cell merging accordingly
    do_merge = Vector{Bool}(undef, psp.nz + extra_cells)
    do_merge .= false
    if psp.model_number > 1
        Threads.@threads for i in 1:psp.nz-1
            a = abs(log10(get_value(psp.eos_res[i].P)) -
                    log10(get_value(psp.eos_res[i+1].P)))
            b = sm.opt.remesh.delta_log10P_split
            if a < 0.03*b && i > psp.nz/2
                # if the condition is satisfied, we merge these cells so long as neither is marked to split
                if !do_split[i] && !do_split[i+1]
                    do_merge[i] = true
                end
            end
        end
    end
    # remove cascading mergers
    for i in 1:psp.nz-1
        if do_merge[i] && do_merge[i+1]
            do_merge[i+1] = false
        end
    end

    # if allocated space is not enough, we need to reallocate everything
    if psp.nz + extra_cells > length(psp.dm)
        adjust_props_size!(sm, psp.nz + extra_cells, sm.nextra)
        # need to reassign psp and ssp since they point to the preallocated ones
        psp = sm.prv_step_props
        ssp = sm.start_step_props
        # have the garbage collector do its job right away
        GC.gc()
    end

    # keep track of original indeces before cell splitting
    # these are needed when doing cell merge
    original_index = zeros(Int64,psp.nz + extra_cells)

    # do splits
    for i=psp.nz:-1:1
        if !do_split[i]
            # move everything upwards by extra_cells
            for j in 1:sm.nvars
                ssp.ind_vars[(i+extra_cells-1)*sm.nvars+j] = psp.ind_vars[(i-1)*sm.nvars+j]
            end
            do_merge[i+extra_cells] = do_merge[i]
            ssp.m[i+extra_cells] = psp.m[i]
            ssp.dm[i+extra_cells] = psp.dm[i]
            original_index[i+extra_cells] = i
        else
            # split the cell and lower extra_cells by 1
            # get old values first
            dm_00 = psp.dm[i]
            var_00 = view(psp.ind_vars, ((i-1)*sm.nvars + 1):((i-1)*sm.nvars + sm.nvars))
            if i > 1
                dm_m1 = psp.dm[i-1]
                var_m1 = view(psp.ind_vars, ((i-2)*sm.nvars + 1):((i-2)*sm.nvars + sm.nvars))
            else
                dm_m1 = NaN
                var_m1 = []
            end
            if i < psp.nz
                dm_p1 = psp.dm[i+1]
                var_p1 = view(psp.ind_vars, ((i)*sm.nvars + 1):((i)*sm.nvars + sm.nvars))
            else
                dm_p1 = NaN
                var_p1 = []
            end
            # views on the new values
            varnew_low = view(ssp.ind_vars,
                              ((i+extra_cells-2)*sm.nvars+1):((i+extra_cells-2)*sm.nvars+sm.nvars))
            varnew_up = view(ssp.ind_vars,
                              ((i+extra_cells-1)*sm.nvars+1):((i+extra_cells-1)*sm.nvars+sm.nvars))
            # populate the new values based on the old ones
            remesh_splitting(sm.equation_set, sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)

            ssp.m[i+extra_cells] = psp.m[i]
            ssp.m[i+extra_cells-1] = psp.m[i]-0.5*psp.dm[i]
            ssp.dm[i+extra_cells] = 0.5*psp.dm[i]
            ssp.dm[i+extra_cells-1] = 0.5*psp.dm[i]

            original_index[i+extra_cells] = i
            original_index[i+extra_cells-1] = i

            extra_cells -= 1
        end
    end

    if (extra_cells != 0)
        throw(AssertionError("extra_cells is not zero: $(extra_cells), cell splitting has gone wrong"))
    end
    ssp.nz = psp.nz + sum(do_split)  # set the new nz

    # do merges
    # we move from center to surface 
    minus_cells = 0
    i = 1

    while i < ssp.nz
        if do_merge[i]
            i_o = original_index[i]
            dm_00 = psp.dm[i_o]
            var_00 = view(psp.ind_vars, ((i_o-1)*sm.nvars + 1):((i_o-1)*sm.nvars + sm.nvars))
            if i > 1
                dm_m1 = psp.dm[i_o-1]
                var_m1 = view(psp.ind_vars, ((i_o-2)*sm.nvars + 1):((i_o-2)*sm.nvars + sm.nvars))
            else
                dm_m1 = NaN
                var_m1 = []
            end
            dm_p1 = psp.dm[i_o+1]
            var_p1 = view(psp.ind_vars, ((i_o)*sm.nvars + 1):((i_o)*sm.nvars + sm.nvars))

            # view on the new values
            varnew = view(ssp.ind_vars,
                              ((i-1)*sm.nvars+1):((i-1)*sm.nvars+sm.nvars))
            remesh_merging(sm.equation_set, sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew)

            ssp.m[i] = ssp.m[i+1]
            ssp.dm[i] = dm_00 + dm_p1
        end
        if minus_cells > 0
            # move everything downwards (including cells that were merged above)
            for j in 1:sm.nvars
                ssp.ind_vars[(i-minus_cells-1)*sm.nvars+j] = ssp.ind_vars[(i-1)*sm.nvars+j]
            end
            ssp.m[i-minus_cells] = ssp.m[i]
            ssp.dm[i-minus_cells] = ssp.dm[i]
        end

        # update counter of merged cells to correctly displace following cells
        if do_merge[i]
            minus_cells += 1
            i += 1 # this will skip the cell above that was merged
        end
        i += 1
    end

    if (minus_cells != sum(do_merge))
        throw(AssertionError("minus_cells is not equal to sum(do_merge): $(minus_cells), $(sum(do_merge)), cell merging has gone wrong"))
    end
    ssp.nz = ssp.nz - sum(do_merge)  # set the new nz

    @show ssp.nz, sum(do_merge), sum(do_split)
end

function split_lnr_lnρ(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)

    lnρ_old = var_00[sm.vari[:lnρ]]  # we keep the same density on both new cells
    lnr_old = var_00[sm.vari[:lnr]]

    # radius of the lower cell is computed with the continuity equation
    # for i=1 we use the r=0 boundary condition as radius
    if i == 1
        r_face_below = 0
    else
        r_face_below = exp(var_m1[sm.vari[:lnr]])
    end
    r_low = (0.5*dm_00*3/(4π*exp(lnρ_old)) + r_face_below^3)^(1/3)

    varnew_low[sm.vari[:lnρ]] = lnρ_old
    varnew_low[sm.vari[:lnr]] = log(r_low)

    varnew_up[sm.vari[:lnρ]] = lnρ_old
    varnew_up[sm.vari[:lnr]] = lnr_old
end

function split_lum(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)
    # for the center face, we use the surface boundary condition L=0
    if i==1
        L_face_below = 0
    else
        L_face_below = var_m1[sm.vari[:lum]]
    end
    L_above = var_00[sm.vari[:lum]]

    varnew_low[sm.vari[:lum]] = 0.5*(L_face_below + L_above)
    varnew_up[sm.vari[:lum]] = L_above
end

function split_lnT(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)
    if i==1
        lnT_low = var_00[sm.vari[:lnT]]  # Central cell remains at the center
        lnT_cell_above = var_p1[sm.vari[:lnT]]

        mnew_up = 0.75*dm_00  # mass from the core
        mcell_above = dm_00 + 0.5*dm_p1

        lnT_up = lnT_low + (lnT_cell_above - lnT_low)*mnew_up/mcell_above
    elseif i==sm.prv_step_props.nz
        lnT_up = var_00[sm.vari[:lnT]]  # Surface remians at same temperature
        lnT_cell_below = var_m1[sm.vari[:lnT]]

        mnew_low = 0.5*dm_m1 + 0.25*dm_00  # mass of new lower cell from center of lower cell
        mup = 0.5*dm_m1 + dm_00 # mass of the surface from center of lower cell

        lnT_low = lnT_cell_below + (lnT_up - lnT_cell_below)*mnew_low/mup
    else
        lnT_cell_above = var_p1[sm.vari[:lnT]]
        lnT_cell_below = var_m1[sm.vari[:lnT]]
        lnT_old = var_00[sm.vari[:lnT]]

        mnew_low = 0.5*dm_m1+0.25*dm_00  # mass at cell center of new lower cell, from center of cell below
        mold = 0.5*dm_m1+0.5*dm_00  # old mass at cell center, from center of cell below
        lnT_low = lnT_cell_below + (lnT_old - lnT_cell_below)*mnew_low/mold

        mnew_up = 0.25*dm_00  # mass from center of cell before splitting
        mcell_above = 0.5*dm_00+0.5*dm_p1
        lnT_up = lnT_old + (lnT_cell_above - lnT_old)*mnew_up/mcell_above
    end
    varnew_low[sm.vari[:lnT]] = lnT_low
    varnew_up[sm.vari[:lnT]] = lnT_up
end

function split_xa(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)
    # use same composition on both cells to preserve species
    for j in 1:sm.network.nspecies
        varnew_low[sm.nvars+1-j] = var_00[sm.nvars+1-j]
        varnew_up[sm.nvars+1-j] = var_00[sm.nvars+1-j]
    end
end

function merge_lnr_lnρ(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew)

    # new density corresponds to mass enclosed within the corresponding radius
    # radius of the lower cell is computed with the continuity equation
    # for i=1 we use the r=0 boundary condition as radius
    if i == 1
        r_face_below = 0
    else
        r_face_below = exp(var_m1[sm.vari[:lnr]])
    end
    r_face_above = exp(var_p1[sm.vari[:lnr]])
    dV = 4pi/3*(r_face_above^3 - r_face_below^3)
    ρ_new = (dm_00 + dm_p1)/dV

    varnew[sm.vari[:lnρ]] = log(ρ_new)
    varnew[sm.vari[:lnr]] = var_p1[sm.vari[:lnr]]
end

function merge_lum(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew)
    L_above = var_p1[sm.vari[:lum]]
    varnew[sm.vari[:lum]] = L_above
end

function merge_lnT(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew)
    lnT_00 = var_00[sm.vari[:lnT]]
    lnT_p1 = var_p1[sm.vari[:lnT]]

    # linear interpolation to new cell center (TODO: verify this)
    new_lnT = (dm_p1*lnT_00 + dm_00*lnT_p1)/(dm_p1 + dm_00)
    varnew[sm.vari[:lnT]] = new_lnT
end

function merge_xa(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew)
    # preserve total mass of each isotope when merging abundance
    new_dm = dm_00 + dm_p1
    for j in 1:sm.network.nspecies
        varnew[sm.nvars+1-j] = (dm_00*var_00[sm.nvars+1-j] + dm_p1*var_p1[sm.nvars+1-j])/new_dm
    end
end

"""
    adjust_props_size(sm, new_nz::Int, nextra::Int)

Returns a new StellarModelProperties object with an adjusted size. Copies over the following from the currect active
properties:

  - nz
  - dt
  - time
  - ind_vars
  - mstar
  - m
  - dm
"""
function adjust_props_size!(sm::StellarModel, new_nz::Int, nextra::Int)
    # verify that new size can contain old sm
    if sm.prv_step_props.nz > new_nz + nextra
        throw(ArgumentError("Can't fit model of size nz=$(sm.prv_step_props.nz) using new_nz=$(new_nz) and nextra=$(nextra)."))
    end

    # Get the type of number being used
    number_type = eltype(sm.props.ind_vars)
    # get dual tag
    # parameters[2] of typeof(sm.props) contains the Dual type used Internally
    # With that, parameters[1] of a ForwardDiff.Dual will give the Tag
    internal_tag = typeof(sm.props).parameters[2].parameters[1]

    # new properties object
    adj_props = build_properties_for_equation_set(sm.equation_set, sm.nvars, new_nz, nextra, sm.network,
                                       sm.vari, number_type, internal_tag)
    # backup scalar quantities
    StellarModels.copy_scalar_properties!(adj_props, sm.start_step_props)
    # copy the mesh quantities (other properties are updated later)
    StellarModels.copy_mesh_properties!(sm, adj_props, sm.start_step_props)
    sm.start_step_props = adj_props

    # also the props used later need new size:
    adj_props = build_properties_for_equation_set(sm.equation_set, sm.nvars, new_nz, nextra, sm.network,
                                       sm.vari, number_type, internal_tag)
    StellarModels.copy_scalar_properties!(adj_props, sm.props)
    StellarModels.copy_mesh_properties!(sm, adj_props, sm.props)
    sm.props = adj_props

    # prv_step_props is also reallocated to be sure all have the same size
    adj_props = build_properties_for_equation_set(sm.equation_set, sm.nvars, new_nz, nextra, sm.network,
                                       sm.vari, number_type, internal_tag)
    StellarModels.copy_scalar_properties!(adj_props, sm.prv_step_props)
    StellarModels.copy_mesh_properties!(sm, adj_props, sm.prv_step_props)
    sm.prv_step_props = adj_props

    # also the solver needs new arrays!
    sm.solver_data = build_solver_data_for_equation_set(sm.equation_set, sm.nvars, new_nz, nextra,
                                sm.solver_data.use_static_arrays, number_type, internal_tag)
end
