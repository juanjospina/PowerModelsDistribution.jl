"admittance model"

function _map_eng2math_mc_line_admittance!(data_math::Dict{String,<:Any}, data_eng::Dict{String,<:Any}; pass_props::Vector{String}=String[])
    if haskey(data_math, "branch")
        for (name, branch) in data_math["branch"]
            z = branch["br_r"] + 1im .* branch["br_x"]
            y_from = branch["g_fr"] + 1im .* branch["b_fr"]
            y_to = branch["g_to"] + 1im .* branch["b_to"]
            z1 = inv(z) + y_from
            z2 = -inv(z)
            z3 = z2
            z4 = inv(z) + y_to
            branch["p_matrix"] = [z1 z2;z3 z4]
        end
    end
end


function _map_eng2math_mc_shunt_admittance!(data_math::Dict{String,<:Any}, data_eng::Dict{String,<:Any}; pass_props::Vector{String}=String[])
    if haskey(data_math, "shunt")
        for (name, shunt) in data_math["shunt"]
            y = shunt["gs"] + 1im .* shunt["bs"]
            y1 = y
            y2 = -y
            y3 = y2
            y4 = y
            shunt["p_matrix"] = [y1 y2;y3 y4]
        end
    end
end


function _map_eng2math_mc_voltage_source_admittance!(data_math::Dict{String,<:Any}, data_eng::Dict{String,<:Any}; pass_props::Vector{String}=String[])
    if haskey(data_math, "gen")
        for (name, gen) in data_math["gen"]
            if occursin("voltage_source", gen["source_id"])
                vsource = data_eng["voltage_source"]["source"]
                z = vsource["rs"] + 1im .* vsource["xs"]
                z1 = inv(z[1:3,1:3])
                z2 = -inv(z[1:3,1:3])
                gen["p_matrix"] = [z1 z2;z2 z1]
            end
        end
    end
end


function _map_eng2math_mc_load_admittance!(data_math::Dict{String,<:Any}, data_eng::Dict{String,<:Any}; pass_props::Vector{String}=String[])
    if haskey(data_math, "load")
        for (name, load) in data_math["load"]
            n = length(load["connections"])
            y = zeros(Complex{Float64}, n, n)
            if load["configuration"] == WYE
                for (i,_i) in enumerate(load["connections"])
                    for (j,_j) in enumerate(load["connections"])
                        if _i != 4 && _j == 4
                            s = conj.(load["pd"][i] + 1im .* load["qd"][i])
                            _y = s / load["vnom_kv"]^2 / 1000
                            y[i,i] += _y
                            y[i,j] -= _y
                            y[j,i] -= _y
                            y[j,j] += _y
                        end
                    end
                end
                load["p_matrix"] = y
            elseif load["configuration"] == DELTA
                for (i,_i) in enumerate(load["connections"])
                    length(load["pd"]) == n ? s = conj.(load["pd"][i] + 1im .* load["qd"][i]) : s = conj.(load["pd"][1] + 1im .* load["qd"][1])
                    for (j,_j) in enumerate(load["connections"])
                        if i != j
                            _y = s / load["vnom_kv"]^2 / 1000
                            y[i,i] += _y
                            y[i,j] -= _y
                        end
                    end
                end
                load["p_matrix"] = y
            end
        end
    end
end


function _map_eng2math_mc_transformer_admittance!(data_math::Dict{String,<:Any}, data_eng::Dict{String,<:Any}; pass_props::Vector{String}=String[])
    lookup = Dict(
        (1,1) => [1,1],
        (1,2) => [5,3],
        (1,3) => [9,5],
        (2,1) => [3,2],
        (2,2) => [7,4],
        (2,3) => [11,6]
    )
    if haskey(data_math, "transformer")
        for (name, transformer) in data_math["transformer"]
        @info "Transformer name: $(name)"
        @info "Transformer MATH data: $(transformer)"
            if transformer["phases"] == 3
                z = sum(transformer["rw"]) + 1im .* transformer["xsc"][1]
                z_1volt= z * 3/transformer["sm_nom"][1]/1000
                z_b = [z_1volt 0 0;0 z_1volt 0;0 0 z_1volt]
                b = [1 0 0;-1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1]
                y1 = b*inv(z_b)*transpose(b)
                n = zeros(Float64, 12, 6)
                a = zeros(Int64,8,12)
                for w = 1:2
                    if transformer["configuration"][w] == WYE
                        w == 1 ? connections = transformer["f_connections"] : connections = transformer["t_connections"]
                        for (_,k) in enumerate(connections)
                            if haskey(lookup, (w,k))
                                i = lookup[(w,k)][1]
                                j = lookup[(w,k)][2]
                                n[i,j] = 1/(transformer["tm_nom"][w]/sqrt(3)*1000*transformer["tm_set"][w][k])
                                n[i+1,j] = - n[i,j]
                            end
                        end
                        if w == 1
                            a[1,1] = a[2,5] = a[3,9] = a[4,2] = a[4,6] = a[4,10] = 1
                        else
                            a[5,3] = a[6,7] = a[7,11] = a[8,4] = a[8,8] = a[8,12] = 1
                        end
                    elseif transformer["configuration"][w] == DELTA
                        w == 1 ? connections = transformer["f_connections"] : connections = transformer["t_connections"]
                        for (_,k) in enumerate(connections)
                            if haskey(lookup, (w,k))
                                i = lookup[(w,k)][1]
                                j = lookup[(w,k)][2]
                                n[i,j] = 1/(transformer["tm_nom"][w]*1000*transformer["tm_set"][w][k])
                                n[i+1,j] = - n[i,j]
                            end
                        end
                        if w == 1
                            a[1,1] = a[1,6] = a[2,5] = a[2,10] = a[3,9] = a[3,2] = 1
                        else
                            a[5,3] = a[6,7] = a[7,11] = a[8,4] = a[8,8] = a[8,12] = 1
                        end
                    end
                end
                y_w = n*y1*transpose(n)
                p_matrix = a*y_w*transpose(a)

                ybase = (transformer["sm_nom"][1]/3) / (transformer["tm_nom"][2]/sqrt(3))^2 /1000
                if haskey(transformer["dss"], "%noloadloss")
                    shunt = (transformer["dss"]["%noloadloss"] - 1im * transformer["dss"]["%imag"])/100*ybase
                    p_matrix[5,5] += shunt
                    p_matrix[5,8] -= shunt
                    p_matrix[6,6] += shunt
                    p_matrix[6,8] -= shunt
                    p_matrix[7,7] += shunt
                    p_matrix[7,8] -= shunt
                    p_matrix[8,5] -= shunt
                    p_matrix[8,6] -= shunt
                    p_matrix[8,7] -= shunt
                    p_matrix[8,8] += 3*shunt
                end

            elseif transformer["phases"] == 1
                z = sum(transformer["rw"]) + 1im .* transformer["xsc"][1]
                z_1volt= z * 1/transformer["sm_nom"][1]/1000
                b = [1 ;-1]
                y1 = b*1/z_1volt*transpose(b)
                n = zeros(Float64, 4, 2)
                a = zeros(Int64,4,4)
                for w = 1:2
                    if transformer["configuration"][w] == WYE
                        i = lookup[(w,1)][1]
                        j = lookup[(w,1)][2]
                        n[i,j] = 1/(transformer["tm_nom"][w]*1000*transformer["tm_set"][w][1])
                        n[i+1,j] = - n[i,j]
                        if w == 1
                            a[1,1] = a[2,2] = 1
                        else
                            a[3,3] = a[4,4] = 1
                        end
                    end
                end
                y_w = n*y1*transpose(n)
                p_matrix = a*y_w*transpose(a)
                ybase = (transformer["sm_nom"][1]) / (transformer["tm_nom"][1])^2 /1000
            end
            transformer["p_matrix"] = p_matrix
        end
    end
end
