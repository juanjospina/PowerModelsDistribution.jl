
function solve_mc_admittance_pf(model::AdmittanceModel)
    z = _SP.sparse(model.z)
    i = _SP.sparse(model.i)
    v = _SP.sparse(model.v)
    c = _SP.sparse(model.c)
    max_it = 10
    it = 1
    _i = i
    _v = deepcopy(v)
    last_v = deepcopy(v)
    delta_i = i
    while it != max_it
        _v = z*_i
        if maximum((abs.(_v-last_v))) < .0001
            break
        else
            _delta_i = c*(abs.(v).^2-abs.(_v).^2)./conj.(_v)
            _i = i + _delta_i
            last_v = _v
            it += 1
        end
    end

    return solution_mc_admittance_pf(_v, it, maximum((abs.(_v-last_v))), i + delta_i, model)
end
