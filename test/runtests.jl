using SensitivityAnalysis, Test, Statistics, UnPack, Accessors

struct SimulateNormal{T}
    μ::T
    σ::T
end

function simulate(object::SimulateNormal)
    @unpack μ, σ = object
    N = 10000
    x = randn(N) .* σ .+ μ
    (mean = mean(x), std = std(x))
end

results = perturbation_results(SimulateNormal(0.0, 1.0), simulate)
push!(results,
      perturbation("μ", @optic(_.μ), ABSOLUTE),
      perturbation("σ", @optic(_.σ), RELATIVE))

moments = [moment("mean", x -> x.mean, ABSOLUTE),
           moment("std", x -> x.std, RELATIVE)]

sens_mu = moment_sensitivity.(results, "μ", moments)

@test sens_mu[1].label == "mean (absolute change)"
@test sens_mu[1].x == results.default_Δs
@test length(sens_mu[1].y) ==length( sens_mu[1].x)
