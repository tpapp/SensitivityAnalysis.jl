"""
"""
module SensitivityAnalysis

export RELATIVE, ABSOLUTE, perturbation, perturbation_results, moment, moment_sensitivity

import Accessors
using ArgCheck: @argcheck
import ThreadPools
using UnPack: @unpack

####
#### changes
####

###
### relative
###

struct Relative end

Base.show(io::IO, ::Relative) = print(io, "relative change")

const RELATIVE = Relative()

(::Relative)(x, r) = x .* (1 + r)

difference(::Relative, x, y) = y / x - 1

###
### absolute
###

struct Absolute end

Base.show(io::IO, ::Absolute) = print(io, "absolute change")

const ABSOLUTE = Absolute()

(::Absolute)(x, r) = x + r

difference(::Absolute, x, y) = y - x

####
#### perturbations
####

struct Perturbation
    label
    metadata
    optic
    change
    captured_errors
    Δs
end

function Base.show(io::IO, perturbation::Perturbation)
    @unpack label, change, captured_errors, Δs = perturbation
    domain = Δs ≡ nothing ? "default domain" : extrema(Δs)
    print(io, "perturb $(label) on $(domain), $(change) [capturing $(captured_errors)]")
end

function perturbation(label, optic, change; metadata = nothing, Δs = nothing, captured_errors = DomainError)
    Perturbation(label, metadata, optic, change, captured_errors, Δs)
end

function perturbation_with(object, perturbation::Perturbation, Δ)
    @unpack optic, change, captured_errors = perturbation
    try
        Accessors.modify(x -> change(x, Δ), object, optic)
    catch e
        if e isa captured_errors
            nothing
        else
            rethrow(e)
        end
    end
end

struct PerturbationResults
    label
    object
    simulate
    baseline_result
    default_Δs
    perturbations_and_results
end

function perturbation_results(object, simulate; baseline_result = simulate(object),
                              label = nothing, default_Δs = range(-0.1, 0.1; length = 11))
    PerturbationResults(label, object, simulate, baseline_result, default_Δs, Vector())
end

function Base.show(io::IO, perturbation_results::PerturbationResults)
    @unpack label, default_Δs, perturbations_and_results = perturbation_results
    print(io, "perturbation results")
    if label ≢ nothing
        print(io, "for ", label)
    end
    println(io, " with default domain $(extrema(default_Δs))")
    for (p, _) in perturbations_and_results
        println(io, "  ", p)
    end
end

Broadcast.broadcastable(x::PerturbationResults) = Ref(x)

function _effective_Δs(perturbation_results, perturbation)
    something(perturbation_results.default_Δs, perturbation.Δs)
end


function Base.push!(perturbation_results::PerturbationResults, perturbation::Perturbation)
    @unpack object, simulate, perturbations_and_results = perturbation_results
    results = ThreadPools.tmap(Δ -> simulate(perturbation_with(object, perturbation, Δ)),
                               _effective_Δs(perturbation_results, perturbation))
    push!(perturbations_and_results, perturbation => results)
    perturbation_results
end

struct Moment
    label
    metadata
    calculator
    change
end

function Base.show(io::IO, moment::Moment)
    @unpack label, change = moment
    print(io, label, " (", change, ")")
end

function moment(label, calculator, change; metadata = nothing)
    Moment(label, metadata, calculator, change)
end

function _moment_sensitivity(moment, baseline_result, results)
    @unpack calculator, change = moment
    m0 = calculator(baseline_result)
    map(results) do r
        if r ≡ nothing
            NaN
        else
            difference(change, calculator(r),  m0)
        end
    end
end

function _lookup_perturbation_by_label(perturbation_results,
                                       pattern::Union{AbstractPattern,AbstractString})
    @unpack perturbations_and_results = perturbation_results
    matches = findall(perturbations_and_results) do (p, _)
        contains(p.label, pattern)
    end
    if length(matches) > 1
        throw(ArgumentError("multiple labels match $(label)"))
    elseif isempty(matches)
        throw(ArgumentError("no labels match $(label)"))
    else
        first(matches)
    end
end

_lookup_perturbation_by_label(perturbation_results, index::Integer) = index

function moment_sensitivity(perturbation_results, index_or_pattern, moment)
    @unpack baseline_result, perturbations_and_results = perturbation_results
    index = _lookup_perturbation_by_label(perturbation_results, index_or_pattern)
    perturbation, results = perturbations_and_results[index]
    sensitivity = _moment_sensitivity(moment, baseline_result, results)
    Δs = _effective_Δs(perturbation_results, perturbation)
    (label = repr(moment), x = Δs, y = sensitivity)
end

end # module
