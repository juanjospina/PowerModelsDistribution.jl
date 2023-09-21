@info "running power flow (pf) admittance-version tests"

@testset "test pf admittance" begin

    @testset "3-bus balanced pf admittance" begin
        pmd_admt = instantiate_mc_admmitance_model(case3_balanced)
        result = solve_mc_admittance_pf(pmd_admt)
    end

    @testset "3-bus unbalanced pf admittance" begin
        pmd_admt = instantiate_mc_admmitance_model(case3_unbalanced)
        result = solve_mc_admittance_pf(pmd_admt)
    end

end
