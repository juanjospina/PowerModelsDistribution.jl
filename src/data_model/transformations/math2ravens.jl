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

    # multinetwork/multiperiod support
    mn_flag = false
    if ismultinetwork(data_math)
        # Assumes there is at least nw=1
        nws_math_data = data_math["nw"]["1"]
        solution = solution_math["nw"]["1"]
        mn_flag = true

        # sort nw keys (ensures vectors are organized in RAVENS)
        keys_array = collect(keys(solution_math["nw"]))
        sorted_nws = sort(keys_array, by=x -> parse(Int, x))

    else
        nws_math_data = data_math
        solution = solution_math
    end

    # Get time_elapsed to compute time steps -- default is 1.0 hour.
    time_elapsed = get(data_math, "time_elapsed", 1.0)

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

    # PowerFlows
    solution_ravens["AnalysisResult"]["OptimalPowerFlow"]["OperationsResult.PowerFlows"] = []

    # Statuses
    solution_ravens["AnalysisResult"]["OptimalPowerFlow"]["OperationsResult.Statuses"] = []

    # Switches States
    solution_ravens["AnalysisResult"]["OptimalPowerFlow"]["OperationsResult.Switches"] = []

    # PowerFlow solutions for Transformers elements
    seen_xfrmrs = Set{String}()     # Set to save xfrmr names

    # Buses/ConnectivityNodes
    for (node_number, node_data) in solution["bus"]

        # Extract the phases for the node
        node_terminals = nws_math_data["bus"][node_number]["terminals"]
        phase_kinds = [phase_mapping[x] for x in node_terminals]

        for (i, result_phase) in enumerate(phase_kinds)

            conn_node = split(nws_math_data["bus"][node_number]["source_id"], '.')[2]

            # Multinetwork support
            if mn_flag == true

                # Curve data
                mn_data = Dict()
                mn_data["AnalysisResultData.Curve"] = Dict()
                mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.xUnit"] = "UnitSymbol.h"
                mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.CurveDatas"] = []

                for nw in sorted_nws
                    mn_info = Dict(
                        "ArCurveData.xvalue" => parse(Float64, nw)*data_math["nw"][nw]["time_elapsed"],
                        "ArCurveData.DataValues" => Dict(
                            "AvVoltage.v" => solution_math["nw"][nw]["bus"][node_number]["vm"][i]*solution_math["nw"][nw]["settings"]["voltage_scale_factor"],
                            "Ravens.cimObjectType" => "AvVoltage",
                        ),
                    )

                    # Conditionally add the angle (some do not have it)
                    if haskey(solution_math["nw"][nw]["bus"][node_number], "va")
                        mn_info["ArCurveData.DataValues"]["AvVoltage.angle"] = solution_math["nw"][nw]["bus"][node_number]["va"][i]
                    end

                    push!(mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.CurveDatas"], mn_info)
                end

                voltage_info = Dict(
                    "AnalysisResultData.phase" => result_phase,
                    "ArVoltage.ConnectivityNode" => "ConnectivityNode::'$(conn_node)'",
                    "AnalysisResultData.Curve" => mn_data["AnalysisResultData.Curve"]
                )

            else
                voltage_info = Dict(
                    "AnalysisResultData.phase" => result_phase,
                    "ArVoltage.ConnectivityNode" => "ConnectivityNode::'$(conn_node)'",
                    "AnalysisResultData.DataValues" => Dict(
                        "AvVoltage.v" => node_data["vm"][i]*solution_math["settings"]["voltage_scale_factor"],
                        "AvVoltage.angle" => node_data["va"][i],
                        "Ravens.cimObjectType" => "AvVoltage",
                    ),
                )
            end

            # Push info to final dictionary
            push!(solution_ravens["AnalysisResult"]["OptimalPowerFlow"]["OperationsResult.Voltages"], voltage_info)

        end

    end

    # Transformers
    for (xfrmr_number, xfrmr_data) in get(solution, "transformer", Dict{Any,Dict{String,Any}}())

        # Extract data
        source_id_vect = split(nws_math_data["transformer"][xfrmr_number]["source_id"], '.')
        cond_eq_type = source_id_vect[2]
        cond_eq_name = source_id_vect[3]
        end_num = parse(Int, source_id_vect[4])

        # Extract the phases for the branch
        terminals = nws_math_data["transformer"][xfrmr_number]["f_connections"]
        phase_kinds = [phase_mapping[x] for x in terminals]

        for (i, result_phase) in enumerate(phase_kinds)

            # Multinetwork support
            if mn_flag == true

                # Curve data
                mn_data = Dict()
                mn_data["AnalysisResultData.Curve"] = Dict()
                mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.xUnit"] = "UnitSymbol.h"
                mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.CurveDatas"] = []

                for nw in sorted_nws
                    mn_info = Dict(
                        "ArCurveData.xvalue" => parse(Float64, nw)*data_math["nw"][nw]["time_elapsed"],
                        "ArCurveData.DataValues" => Dict(
                            "AvPowerFlow.p" => solution_math["nw"][nw]["transformer"][xfrmr_number]["pf"][i]*solution_math["nw"][nw]["settings"]["power_scale_factor"],
                            "AvPowerFlow.q" => solution_math["nw"][nw]["transformer"][xfrmr_number]["qf"][i]*solution_math["nw"][nw]["settings"]["power_scale_factor"],
                            "AvPowerFlow.endNumber" => end_num,
                            "Ravens.cimObjectType" => "ArPowerFlow",
                        ),
                    )
                    push!(mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.CurveDatas"], mn_info)
                end

                pf_info = Dict(
                    "AnalysisResultData.phase" => result_phase,
                    "ArPowerFlow.ConductingEquipment" => "$(cond_eq_type)::'$(cond_eq_name)'",
                    "AnalysisResultData.Curve" => mn_data["AnalysisResultData.Curve"]
                )

            else

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


            end

            push!(solution_ravens["AnalysisResult"]["OptimalPowerFlow"]["OperationsResult.PowerFlows"], pf_info)

        end

        # Status for Xfrmrs
        if !(cond_eq_name in seen_xfrmrs)

            xfrmr_status = nws_math_data["transformer"][xfrmr_number]["status"] == 1 ? true : false

            # Multinetwork support
            if mn_flag == true

                # Curve data
                mn_data = Dict()
                mn_data["AnalysisResultData.Curve"] = Dict()
                mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.xUnit"] = "UnitSymbol.h"
                mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.CurveDatas"] = []

                for nw in sorted_nws
                    mn_info = Dict(
                        "ArCurveData.xvalue" => parse(Float64, nw)*data_math["nw"][nw]["time_elapsed"],
                        "ArCurveData.DataValues" => Dict(
                            "AvStatus.inService" => xfrmr_status
                        ),
                    )
                    push!(mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.CurveDatas"], mn_info)
                end

                status_info = Dict(
                    "ArStatus.ConductingEquipment" => "$(cond_eq_type)::'$(cond_eq_name)'",
                    "AnalysisResultData.Curve" => mn_data["AnalysisResultData.Curve"]
                )

            else

                status_info = Dict(
                    "ArStatus.ConductingEquipment" => "$(cond_eq_type)::'$(cond_eq_name)'",
                    "AnalysisResultData.DataValues" => Dict(
                        "AvStatus.inService" => xfrmr_status,
                    )
                )

            end

            push!(solution_ravens["AnalysisResult"]["OptimalPowerFlow"]["OperationsResult.Statuses"], status_info)

        end

        # Store the xfrmr name
        push!(seen_xfrmrs, "$(cond_eq_name)")

    end

    # Edge elements (branches, switches)
    edge_elements = ["branch", "switch"]
    for edge_elmnt in edge_elements

        for (edge_number, edge_data) in get(solution, edge_elmnt, Dict{Any,Dict{String,Any}}())

            # Filter virtual elements that exist in the MATH model
            if !occursin("virtual", nws_math_data[edge_elmnt][edge_number]["name"])

                cond_eq_type = split(nws_math_data[edge_elmnt][edge_number]["source_id"], '.')[1]
                cond_eq_name = split(nws_math_data[edge_elmnt][edge_number]["source_id"], '.')[2]

                # Add Switch state (only switches)
                if edge_elmnt == "switch"

                    # initial state from math network data
                    sw_state = nws_math_data[edge_elmnt][edge_number]["state"] == 1 ? false : true

                    if mn_flag == true

                        # Curve data
                        mn_data = Dict()
                        mn_data["AnalysisResultData.Curve"] = Dict()
                        mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.xUnit"] = "UnitSymbol.h"
                        mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.CurveDatas"] = []

                        for nw in sorted_nws
                            if haskey(solution_math["nw"][nw][edge_elmnt][edge_number], "state")
                                sw_state = Int(solution_math["nw"][nw][edge_elmnt][edge_number]["state"]) == 0 ? true : false
                            end
                            mn_info = Dict(
                                "ArCurveData.xvalue" => parse(Float64, nw)*data_math["nw"][nw]["time_elapsed"],
                                "ArCurveData.DataValues" => Dict(
                                    "AvSwitch.open" => sw_state,
                                ),
                            )
                            push!(mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.CurveDatas"], mn_info)
                        end

                        state_info = Dict(
                            "ArSwitch.Switch" => "$(cond_eq_type)::'$(cond_eq_name)'",
                            "AnalysisResultData.Curve" => mn_data["AnalysisResultData.Curve"]
                        )

                    else
                        state_info = Dict(
                            "ArSwitch.Switch" => "$(cond_eq_type)::'$(cond_eq_name)'",
                            "AnalysisResultData.DataValues" => Dict(
                                "AvSwitch.open" => sw_state,
                            )
                        )

                    end

                    push!(solution_ravens["AnalysisResult"]["OptimalPowerFlow"]["OperationsResult.Switches"], state_info)

                end


                num_ends = 2
                # # OPTIONAL opt out of edge elements to write from and to flows
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
                    terminals = nws_math_data[edge_elmnt][edge_number][connection_flow]
                    phase_kinds = [phase_mapping[x] for x in terminals]

                    for (i, result_phase) in enumerate(phase_kinds)

                        # Multinetwork support
                        if mn_flag == true

                            # Curve data
                            mn_data = Dict()
                            mn_data["AnalysisResultData.Curve"] = Dict()
                            mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.xUnit"] = "UnitSymbol.h"
                            mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.CurveDatas"] = []

                            for nw in sorted_nws
                                mn_info = Dict(
                                    "ArCurveData.xvalue" => parse(Float64, nw)*data_math["nw"][nw]["time_elapsed"],
                                    "ArCurveData.DataValues" => Dict(
                                        "AvPowerFlow.p" => solution_math["nw"][nw][edge_elmnt][edge_number][p_flow_direction][i]*solution_math["nw"][nw]["settings"]["power_scale_factor"],
                                        "AvPowerFlow.q" => solution_math["nw"][nw][edge_elmnt][edge_number][q_flow_direction][i]*solution_math["nw"][nw]["settings"]["power_scale_factor"],
                                        "AvPowerFlow.endNumber" => end_num,
                                        "Ravens.cimObjectType" => "ArPowerFlow",
                                    ),
                                )
                                push!(mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.CurveDatas"], mn_info)
                            end

                            pf_info = Dict(
                                "AnalysisResultData.phase" => result_phase,
                                "ArPowerFlow.ConductingEquipment" => "$(cond_eq_type)::'$(cond_eq_name)'",
                                "AnalysisResultData.Curve" => mn_data["AnalysisResultData.Curve"]
                            )

                        else
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

                        end

                        push!(solution_ravens["AnalysisResult"]["OptimalPowerFlow"]["OperationsResult.PowerFlows"], pf_info)

                    end
                end

                # Statuses
                object_prefix = ""
                if edge_elmnt == "branch"
                    object_prefix = "br_"
                end

                elemtn_status = nws_math_data[edge_elmnt][edge_number]["$(object_prefix)status"] == 1 ? true : false

                # Multinetwork support
                if mn_flag == true

                    # Curve data
                    mn_data = Dict()
                    mn_data["AnalysisResultData.Curve"] = Dict()
                    mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.xUnit"] = "UnitSymbol.h"
                    mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.CurveDatas"] = []

                    for nw in sorted_nws
                        mn_info = Dict(
                            "ArCurveData.xvalue" => parse(Float64, nw)*data_math["nw"][nw]["time_elapsed"],
                            "ArCurveData.DataValues" => Dict(
                                "AvStatus.inService" => elemtn_status
                            ),
                        )
                        push!(mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.CurveDatas"], mn_info)
                    end

                    status_info = Dict(
                        "ArStatus.ConductingEquipment" => "$(cond_eq_type)::'$(cond_eq_name)'",
                        "AnalysisResultData.Curve" => mn_data["AnalysisResultData.Curve"]
                    )

                else

                    status_info = Dict(
                        "ArStatus.ConductingEquipment" => "$(cond_eq_type)::'$(cond_eq_name)'",
                        "AnalysisResultData.DataValues" => Dict(
                            "AvStatus.inService" => elemtn_status,
                        )
                    )

                end

                push!(solution_ravens["AnalysisResult"]["OptimalPowerFlow"]["OperationsResult.Statuses"], status_info)

            end

        end
    end

    # Nodal elements (loads, gens)
    node_elements = ["load", "gen", "storage"]
    for node_elmnt in node_elements

        for (node_number, node_data) in get(solution, node_elmnt, Dict{Any,Dict{String,Any}}())

            # Filter virtual elements that exist in the MATH model
            if !occursin("virtual", nws_math_data[node_elmnt][node_number]["name"])

                cond_eq_type = split(nws_math_data[node_elmnt][node_number]["source_id"], '.')[1]
                cond_eq_name = split(nws_math_data[node_elmnt][node_number]["source_id"], '.')[2]

                if node_elmnt == "load"
                    p_key = "pd"
                    q_key = "qd"
                elseif node_elmnt == "gen"
                    p_key = "pg"
                    q_key = "qg"
                elseif node_elmnt == "storage"
                    p_key = "ps"
                    q_key = "qs"
                else
                    p_key = "p"
                    q_key = "q"
                end

                # Extract the phases for the node element
                terminals = nws_math_data[node_elmnt][node_number]["connections"]
                phase_kinds = [phase_mapping[x] for x in terminals]

                for (i, result_phase) in enumerate(phase_kinds)

                    # Multinetwork support
                    if mn_flag == true

                        # Curve data
                        mn_data = Dict()
                        mn_data["AnalysisResultData.Curve"] = Dict()
                        mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.xUnit"] = "UnitSymbol.h"
                        mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.CurveDatas"] = []

                        for nw in sorted_nws
                            mn_info = Dict(
                                "ArCurveData.xvalue" => parse(Float64, nw)*data_math["nw"][nw]["time_elapsed"],
                                "ArCurveData.DataValues" => Dict(
                                    "AvPowerFlow.p" => solution_math["nw"][nw][node_elmnt][node_number][p_key][i]*solution_math["nw"][nw]["settings"]["power_scale_factor"],
                                    "AvPowerFlow.q" => solution_math["nw"][nw][node_elmnt][node_number][q_key][i]*solution_math["nw"][nw]["settings"]["power_scale_factor"],
                                    "Ravens.cimObjectType" => "ArPowerFlow",
                                ),
                            )
                            push!(mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.CurveDatas"], mn_info)
                        end

                        pf_info = Dict(
                            "AnalysisResultData.phase" => result_phase,
                            "ArPowerFlow.ConductingEquipment" => "$(cond_eq_type)::'$(cond_eq_name)'",
                            "AnalysisResultData.Curve" => mn_data["AnalysisResultData.Curve"]
                        )


                    else

                        pf_info = Dict(
                            "AnalysisResultData.phase" => result_phase,
                            "ArPowerFlow.ConductingEquipment" => "$(cond_eq_type)::'$(cond_eq_name)'",
                            "AnalysisResultData.DataValues" => Dict(
                                "AvPowerFlow.p" => node_data[p_key][i]*solution_math["settings"]["power_scale_factor"],
                                "AvPowerFlow.q" => node_data[q_key][i]*solution_math["settings"]["power_scale_factor"],
                                "Ravens.cimObjectType" => "ArPowerFlow",
                            ),
                        )

                    end

                    push!(solution_ravens["AnalysisResult"]["OptimalPowerFlow"]["OperationsResult.PowerFlows"], pf_info)

                end


                # Statuses
                object_prefix = ""
                if node_elmnt == "gen"
                    object_prefix = "gen_"
                end

                # original status from math data
                elemtn_status = nws_math_data[node_elmnt][node_number]["$(object_prefix)status"] == 1 ? true : false

                # Multinetwork support
                if mn_flag == true

                    # Curve data
                    mn_data = Dict()
                    mn_data["AnalysisResultData.Curve"] = Dict()
                    mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.xUnit"] = "UnitSymbol.h"
                    mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.CurveDatas"] = []

                    for nw in sorted_nws
                        if haskey(solution_math["nw"][nw][node_elmnt][node_number], "status")
                            elemtn_status = Int(solution_math["nw"][nw][node_elmnt][node_number]["status"]) == 1 ? true : false
                        end
                        mn_info = Dict(
                            "ArCurveData.xvalue" => parse(Float64, nw)*data_math["nw"][nw]["time_elapsed"],
                            "ArCurveData.DataValues" => Dict(
                                "AvStatus.inService" => elemtn_status
                            ),
                        )
                        push!(mn_data["AnalysisResultData.Curve"]["AnalysisResultCurve.CurveDatas"], mn_info)
                    end

                    status_info = Dict(
                        "ArStatus.ConductingEquipment" => "$(cond_eq_type)::'$(cond_eq_name)'",
                        "AnalysisResultData.Curve" => mn_data["AnalysisResultData.Curve"]
                    )

                else

                    status_info = Dict(
                        "ArStatus.ConductingEquipment" => "$(cond_eq_type)::'$(cond_eq_name)'",
                        "AnalysisResultData.DataValues" => Dict(
                            "AvStatus.inService" => elemtn_status,
                        )
                    )

                end

                push!(solution_ravens["AnalysisResult"]["OptimalPowerFlow"]["OperationsResult.Statuses"], status_info)

            end
        end
    end

    return solution_ravens

end
