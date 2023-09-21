function solution_mc_admittance_pf(v::_SP.SparseMatrixCSC{ComplexF64, Int64}, it::Int64, last_delta::Float64, i::_SP.SparseMatrixCSC{ComplexF64, Int64}, model)
    solution = Dict{String, Any}()
    solution["bus"] = Dict{String, Any}()
    for (indx,bus) in model.data["bus"]
        solution["bus"][indx] = Dict{String, Any}(
            "vm" => [0.0 for t in bus["terminals"]],
            "va" => [0.0 for t in bus["terminals"]],
            "name" => bus["source_id"],
            "vbase" => bus["vbase"],
            "im" => [0.0 for t in bus["terminals"]],
            "ia" => [0.0 for t in bus["terminals"]],
        )
        for (j, grounded) in enumerate(bus["grounded"])
            if grounded == 0
                t = bus["terminals"][j]
                solution["bus"][indx]["vm"][j] = abs(v[model.data["admittance_map"][(bus["index"], t)]])
                solution["bus"][indx]["va"][j] = angle(v[model.data["admittance_map"][(bus["index"], t)]]) * 180/pi
                solution["bus"][indx]["im"][j] = abs(i[model.data["admittance_map"][(bus["index"], t)]])
                solution["bus"][indx]["ia"][j] = angle(i[model.data["admittance_map"][(bus["index"], t)]]) * 180/pi
            end
        end
        solution["model"] = model
        solution["solver"] = Dict{String,Any}(
            "it" => it,
            "delta" => last_delta,
        )
    end

    return solution
end
