using SensitivityAnalysis, Test, Statistics, SimpleUnPack, Accessors

struct SimulateNormal{T}
    μ::T
    σ::T
    function SimulateNormal(μ::T, σ::T) where {T}
        σ > 1.001 && throw(DomainError("let's tigger catch"))
        new{T}(μ, σ)
    end
end

function simulate(object::SimulateNormal)
    @unpack μ, σ = object
    N = 10000
    x = randn(N) .* σ .+ μ
    (mean = mean(x), std = std(x))
end

results = perturbation_analysis(SimulateNormal(0.0, 1.0), simulate)
push!(results,
      perturbation("μ", @optic(_.μ), ABSOLUTE),
      perturbation("σ", @optic(_.σ), RELATIVE))
@test sort(keys(results)) == sort(["μ", "σ"])

moments = [moment("mean", x -> x.mean, ABSOLUTE; metadata = :meta),
           moment("std", x -> x.std, RELATIVE)]

sens_mu = moment_sensitivity.(results, "μ", moments)
let s = sens_mu[1]
    @test s.label == "mean (absolute change)"
    @test s.x == results.default_Δs
    @test length(s.y) ==length(s.x)
    @test s.metadata ≡ :meta
end

sens_sigma = moment_sensitivity.(results, "σ", moments)
@test isnan(sens_sigma[1].y[end])

@test_throws ArgumentError moment_sensitivity(results, "not found", moments[1])
