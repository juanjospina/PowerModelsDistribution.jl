function transform_solution_ravens(
    solution_math::Dict{String,<:Any},
    data_math::Dict{String,<:Any};
    map::Union{Vector{<:Dict{String,<:Any}},Missing}=missing,
    make_si::Bool=true,
    convert_rad2deg::Bool=true,
     map_math2eng_extensions::Dict{String,<:Function}=Dict{String,Function}(),
    make_si_extensions::Vector{<:Function}=Function[],
    dimensionalize_math_extensions::Dict{String,<:Dict{String,<:Vector{<:String}}}=Dict{String,Dict{String,Vector{String}}}()
)::Dict{String,Any}

    @assert ismath(data_math) "cannot be converted. Not a MATH model."

    # convert solution to si?
    solution_math = solution_make_si(
        solution_math,
        data_math;
        mult_vbase=make_si,
        mult_sbase=make_si,
        convert_rad2deg=convert_rad2deg,
        make_si_extensions=make_si_extensions,
        dimensionalize_math_extensions=dimensionalize_math_extensions
    )

    # TODO: multinetwork/multiperiod support
    # if ismultinetwork(data_math)
    #     nws_math_sol = get(solution_math, "nw", Dict{String,Any}())
    #     nws_math_data = data_math["nw"]
    # else
    #     nws_math_sol = Dict("0" => solution_math)
    #     nws_math_data = Dict("0" => data_math)
    # end

    # Create OptimalPowerFlow AnalysisResult Dictionary
    # TODO: read original JSON file, and add result to that JSON file (check if AnalysisResult exists before overwriting it)
    solution_ravens = Dict()
    solution_ravens["AnalysisResult"] = Dict()
    solution_ravens["AnalysisResult"]["OptimalPowerFlow"] = Dict()
    solution_ravens["AnalysisResult"]["OptimalPowerFlow"]["Ravens.cimObjectType"] = "OperationsResult"
    solution_ravens["AnalysisResult"]["OptimalPowerFlow"]["IdentifiedObject.name"] = "OptimalPowerFlow"
    solution_ravens["AnalysisResult"]["OptimalPowerFlow"]["IdentifiedObject.mRID"] = "#_$(uppercase(string(UUIDs.uuid4())))"

    # Create a mapping from integers to strings
    phase_mapping = Dict(1 => "SinglePhaseKind.A", 2 => "SinglePhaseKind.B", 3 => "SinglePhaseKind.C")

    # Voltages
    solution_ravens["AnalysisResult"]["OptimalPowerFlow"]["OperationsResult.Voltages"] = []
    for (node_number, node_data) in solution_math["bus"]
        # Extract the phases for the node
        node_terminals = data_math["bus"][node_number]["terminals"]
        phase_kinds = [phase_mapping[x] for x in node_terminals]

        for (i, result_phase) in enumerate(phase_kinds)
            conn_node = split(data_math["bus"][node_number]["source_id"], '.')[2]
            voltage_info = Dict(
                "AnalysisResultData.phase" => result_phase,
                "ArVoltage.ConnectivityNode" => "ConnectivityNode::'$(conn_node)'",
                "AnalysisResultData.DataValues" => Dict(
                    "AvVoltage.v" => node_data["vm"][i]*solution_math["settings"]["voltage_scale_factor"],
                    "AvVoltage.angle" => node_data["va"][i],
                    "Ravens.cimObjectType" => "AvVoltage",
                ),
            )
            push!(solution_ravens["AnalysisResult"]["OptimalPowerFlow"]["OperationsResult.Voltages"], voltage_info)
        end
    end

    # PowerFlows
    solution_ravens["AnalysisResult"]["OptimalPowerFlow"]["OperationsResult.PowerFlows"] = []

    # Statuses
    solution_ravens["AnalysisResult"]["OptimalPowerFlow"]["OperationsResult.Statuses"] = []

    # PowerFlow solutions for Transformers elements
    for (xfrmr_number, xfrmr_data) in get(solution_math, "transformer", Dict{Any,Dict{String,Any}}())

        source_id_vect = split(data_math["transformer"][xfrmr_number]["source_id"], '.')
        cond_eq_type = source_id_vect[2]
        cond_eq_name = source_id_vect[3]
        end_num = parse(Int, source_id_vect[4])

        # Extract the phases for the branch
        terminals = data_math["transformer"][xfrmr_number]["f_connections"]
        phase_kinds = [phase_mapping[x] for x in terminals]

        for (i, result_phase) in enumerate(phase_kinds)

            pf_info = Dict(
                "AnalysisResultData.phase" => result_phase,
                "ArPowerFlow.ConductingEquipment" => "$(cond_eq_type)::'$(cond_eq_name)'",
                "AnalysisResultData.DataValues" => Dict(
                    "AvPowerFlow.p" => xfrmr_data["pf"][i]*solution_math["settings"]["power_scale_factor"],
                    "AvPowerFlow.q" => xfrmr_data["qf"][i]*solution_math["settings"]["power_scale_factor"],
                    "AvPowerFlow.endNumber" => end_num,
                    "Ravens.cimObjectType" => "ArPowerFlow",
                ),
            )
            push!(solution_ravens["AnalysisResult"]["OptimalPowerFlow"]["OperationsResult.PowerFlows"], pf_info)

        end
    end


    edge_elements = ["branch", "switch"]
    for edge_elmnt in edge_elements

        for (edge_number, edge_data) in get(solution_math, edge_elmnt, Dict{Any,Dict{String,Any}}())

            # Filter virtual elements that exist in the MATH model
            if !occursin("virtual", data_math[edge_elmnt][edge_number]["name"])

                cond_eq_type = split(data_math[edge_elmnt][edge_number]["source_id"], '.')[1]
                cond_eq_name = split(data_math[edge_elmnt][edge_number]["source_id"], '.')[2]

                num_ends = 2
                # # OPTIONAL opt out of edge elements beside transformers to write from and to flows
                # if edge_elmnt != "transformer"
                #     num_ends = 1
                # end

                for end_num in 1:num_ends # loop through ends
                    if end_num == 1
                        connection_flow = "f_connections"
                        p_flow_direction = "pf"
                        q_flow_direction = "qf"
                    else end_num == 2
                        connection_flow = "t_connections"
                        p_flow_direction = "pt"
                        q_flow_direction = "qt"
                    end

                    # Extract the phases for the branch
                    terminals = data_math[edge_elmnt][edge_number][connection_flow]
                    phase_kinds = [phase_mapping[x] for x in terminals]

                    for (i, result_phase) in enumerate(phase_kinds)

                        pf_info = Dict(
                            "AnalysisResultData.phase" => result_phase,
                            "ArPowerFlow.ConductingEquipment" => "$(cond_eq_type)::'$(cond_eq_name)'",
                            "AnalysisResultData.DataValues" => Dict(
                                "AvPowerFlow.p" => edge_data[p_flow_direction][i]*solution_math["settings"]["power_scale_factor"],
                                "AvPowerFlow.q" => edge_data[q_flow_direction][i]*solution_math["settings"]["power_scale_factor"],
                                "AvPowerFlow.endNumber" => end_num,
                                "Ravens.cimObjectType" => "ArPowerFlow",
                            ),
                        )
                        push!(solution_ravens["AnalysisResult"]["OptimalPowerFlow"]["OperationsResult.PowerFlows"], pf_info)

                    end
                end

                # Statuses
                object_prefix = ""
                if edge_elmnt == "branch"
                    object_prefix = "br_"
                end

                elemtn_status = data_math[edge_elmnt][edge_number]["$(object_prefix)status"] == 1 ? true : false
                status_info = Dict(
                    "ArStatus.ConductingEquipment" => "$(cond_eq_type)::'$(cond_eq_name)'",
                    "AnalysisResultData.DataValues" => Dict(
                        "AvStatus.inService" => elemtn_status,
                    )
                )
                push!(solution_ravens["AnalysisResult"]["OptimalPowerFlow"]["OperationsResult.Statuses"], status_info)

            end

        end
    end

    #  PowerFlow solutions for Node elements
    node_elements = ["load", "gen"]
    for node_elmnt in node_elements

        for (node_number, node_data) in get(solution_math, node_elmnt, Dict{Any,Dict{String,Any}}())

            # Filter virtual elements that exist in the MATH model
            if !occursin("virtual", data_math[node_elmnt][node_number]["name"])

                cond_eq_type = split(data_math[node_elmnt][node_number]["source_id"], '.')[1]
                cond_eq_name = split(data_math[node_elmnt][node_number]["source_id"], '.')[2]

                if node_elmnt == "load"
                    p_key = "pd"
                    q_key = "qd"
                elseif node_elmnt == "gen"
                    p_key = "pg"
                    q_key = "qg"
                else
                    p_key = "p"
                    q_key = "q"
                end

                # Extract the phases for the node element
                terminals = data_math[node_elmnt][node_number]["connections"]
                phase_kinds = [phase_mapping[x] for x in terminals]

                for (i, result_phase) in enumerate(phase_kinds)

                    pf_info = Dict(
                        "AnalysisResultData.phase" => result_phase,
                        "ArPowerFlow.ConductingEquipment" => "$(cond_eq_type)::'$(cond_eq_name)'",
                        "AnalysisResultData.DataValues" => Dict(
                            "AvPowerFlow.p" => node_data[p_key][i]*solution_math["settings"]["power_scale_factor"],
                            "AvPowerFlow.q" => node_data[q_key][i]*solution_math["settings"]["power_scale_factor"],
                            "Ravens.cimObjectType" => "ArPowerFlow",
                        ),
                    )
                    push!(solution_ravens["AnalysisResult"]["OptimalPowerFlow"]["OperationsResult.PowerFlows"], pf_info)
                end


                # Statuses
                object_prefix = ""
                if node_elmnt == "gen"
                    object_prefix = "gen_"
                end

                elemtn_status = data_math[node_elmnt][node_number]["$(object_prefix)status"] == 1 ? true : false
                status_info = Dict(
                    "ArStatus.ConductingEquipment" => "$(cond_eq_type)::'$(cond_eq_name)'",
                    "AnalysisResultData.DataValues" => Dict(
                        "AvStatus.inService" => elemtn_status,
                    )
                )
                push!(solution_ravens["AnalysisResult"]["OptimalPowerFlow"]["OperationsResult.Statuses"], status_info)
            end
        end
    end


    # @info "$(solution_ravens)"

    # open("./OPFTEST-mod.json","w") do f
    #     JSON.print(f, solution_ravens, 2)
    # end

    # asdasd

end
