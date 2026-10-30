struct StandardRemeshPolicy <: AbstractRemeshPolicy
end

struct NeverRemeshPolicy <: AbstractRemeshPolicy
end

function remesh_model!( remesh ::StandardRemeshPolicy, sm::StellarModel, state::StellarEvolutionState)
    should_remesh = sm.opt.remesh.do_remesh && state.use_remeshing

    if should_remesh
        StellarModels.remesher!(sm)
    else
        StellarModels.copy_mesh_properties!( sm, sm.start_step_props, sm.prv_step_props)
    end
end 