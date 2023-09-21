"admittance model"
const _mc_admittance_asset_types = [
    "line_admittance", "voltage_source_admittance", "load_admittance", "transformer_admittance", "shunt_admittance"
]

"custom version of 'transform_data_model' to build admittance model and deal with transformers"
function transform_admittance_data_model(
    data::Dict{String,<:Any};
    global_keys::Set{String}=Set{String}(),
    eng2math_passthrough::Dict{String,<:Vector{<:String}}=Dict{String,Vector{String}}(),
    eng2math_extensions::Vector{<:Function}=Function[],
    make_pu::Bool=true,
    make_pu_extensions::Vector{<:Function}=Function[],
    build_model::Bool=false,
    correct_network_data::Bool=true,
    kwargs...,
    )::Dict{String,Any}

    data_math = _map_eng2math_mc_admittance(
        data;
        eng2math_extensions = eng2math_extensions,
        eng2math_passthrough = eng2math_passthrough,
        make_pu_extensions = make_pu_extensions,
        global_keys = global_keys,
        build_model = build_model,
        kwargs...
    )

    correct_network_data && correct_network_data!(data_math; make_pu=false, make_pu_extensions=make_pu_extensions)

    correct_grounds!(data_math)
    populate_bus_voltages!(data_math)

    return data_math

end


"base function for converting mc engineering model to mathematical with admittances"
function _map_eng2math_mc_admittance(
    data_eng::Dict{String,<:Any};
    eng2math_passthrough::Dict{String,<:Vector{<:String}}=Dict{String,Vector{String}}(),
    eng2math_extensions::Vector{<:Function}=Function[],
    make_pu::Bool=true,
    make_pu_extentions::Vector{<:Function}=Function[],
    global_keys::Set{String}=Set{String}(),
    build_model::Bool=false,
    kwargs...,
    )::Dict{String,Any}

    _data_eng = deepcopy(data_eng)

    # any pre-processing of data here

    # TODO kron

    # TODO phase projection

    if ismultinetwork(data_eng)
        #  TODO multi network
    else
        data_math = Dict{String,Any}(
            "name" => get(_data_eng, "name", ""),
            "per_unit" => get(_data_eng, "per_unit", false),
            "data_model" => MATHEMATICAL,
            "is_projected" => get(_data_eng, "is_projected", false),
            "is_kron_reduced" => get(_data_eng, "is_kron_reduced", false),
            "settings" => deepcopy(_data_eng["settings"]),
            "time_elapsed" => get(_data_eng, "time_elapsed", 1.0),
        )
    end

    _map_eng2math_nw!(data_math, data_eng, eng2math_passthrough=eng2math_passthrough, eng2math_extensions=eng2math_extensions)

    _apply_mc_admittance!(_map_eng2math_mc_admittance_nw!, data_math, _data_eng; eng2math_passthrough=eng2math_passthrough, eng2math_extensions=eng2math_extensions)

    # admittance_bus_order!(data_math)
    return data_math
end


function _apply_mc_admittance!(func!::Function, data1::Dict{String,<:Any}, data2::Dict{String,<:Any}; kwargs...)
    func!(data1, data2; kwargs...)
end


function _map_eng2math_mc_admittance_nw!(data_math::Dict{String,<:Any}, data_eng::Dict{String,<:Any}; eng2math_passthrough::Dict{String,<:Vector{<:String}}=Dict{String,Vector{String}}(), eng2math_extensions::Vector{<:Function}=Function[])

    for type in _mc_admittance_asset_types # --> anything from missing from the model needed for the solve or admittance matrix maybe per unit to actual
        getfield(PowerModelsDistribution, Symbol("_map_eng2math_mc_$(type)!"))(data_math, data_eng; pass_props=get(eng2math_passthrough, type, String[]))
    end

    for eng2math_func! in eng2math_extensions
        eng2math_func!(data_math, data_eng)
    end
end


"mod with out per unit corrections see: common.jl in io PowerModelsDistribution"
function correct_network_data_admittance!(data::Dict{String,Any})
    if iseng(data)
        check_eng_data_model(data)
    elseif ismath(data)

        check_connectivity(data)
        correct_branch_directions!(data)
        check_branch_loops(data)
        correct_bus_types!(data)
        propagate_network_topology!(data)
    end
    return data
end

function populate_bus_voltages!(data::Dict{String,Any})
    for (i, transformer) in data["transformer"]
        f_bus = transformer["f_bus"]
        t_bus = transformer["t_bus"]
        if haskey(transformer, "tm_nom")
            transformer["dss"]["phases"] == 3 ? multi = 1/sqrt(3) : multi = 1
            if !haskey(data["bus"][string(f_bus)], "vbase")
                data["bus"][string(f_bus)]["vbase"] = transformer["tm_nom"][1]*multi
            end
            if !haskey(data["bus"][string(t_bus)], "vbase")
                data["bus"][string(t_bus)]["vbase"] = transformer["tm_nom"][2]*multi
            end
        end
    end

    propagate_voltages!(data)

    for (i, gen) in data["gen"]
        if !haskey(data["bus"][string(gen["gen_bus"])], "vbase")
            if occursin("voltage_source", gen["source_id"])
                data["bus"][string(gen["gen_bus"])]["vbase"] = gen["vg"][1]
            end
        end
    end
    propagate_voltages!(data)

end


function propagate_voltages!(data::Dict{String,Any})
    buses = collect(keys(data["bus"]))
    for (i,bus) in data["bus"]
        !haskey(bus, "vbase") ? filter!(n->n != string(i), buses) : nothing
    end

    found = true
    while found
        found = false
        for (i, branch) in data["branch"]
            f_bus = string(branch["f_bus"])
            t_bus = string(branch["t_bus"])
            if f_bus in buses
                if t_bus in buses
                    nothing
                else
                    data["bus"][t_bus]["vbase"] = data["bus"][f_bus]["vbase"]
                    found = true
                    push!(buses, t_bus)
                end
            elseif t_bus in buses
                if f_bus in buses
                    nothing
                else
                    data["bus"][f_bus]["vbase"] = data["bus"][t_bus]["vbase"]
                    found = true
                    push!(buses, f_bus)
                end
            end
        end
    end
end


function correct_grounds!(data::Dict{String,Any})
    for (i, transformer) in data["transformer"]
        for (i,config) in enumerate(transformer["configuration"])
            if config == WYE
                if occursin(".1.2.3.0", transformer["dss"]["buses"][i])
                    if i == 2
                        transformer["t_connections"] = [1,2,3,4]
                        data["bus"][string(transformer["t_bus"])]["terminals"] = [1,2,3,4]
                        data["bus"][string(transformer["t_bus"])]["grounded"] = Bool[0,0,0,1]
                    end
                end
            end
        end
    end
end
