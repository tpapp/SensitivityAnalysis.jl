"""
Organize sensitivity / perturbation analysis.

See [`perturbation_analysis`](@ref) for a short example, and the documentation for details.
"""
module SensitivityAnalysis

export RELATIVE, ABSOLUTE, perturbation, perturbation_analysis, moment, moment_sensitivity

import Accessors
using ArgCheck: @argcheck
using DocStringExtensions: SIGNATURES
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

"Represents a relative change, for inputs (`y = x * (1 + Δ)`) or outputs (`Δ = y / x - 1`)."
const RELATIVE = Relative()

(::Relative)(x, Δ) = x .* (1 + Δ)

difference(::Relative, y, x) = y / x - 1

###
### absolute
###

struct Absolute end

Base.show(io::IO, ::Absolute) = print(io, "absolute change")

"Represents an absolute, for inputs (`y = x + Δ`) or outputs (`Δ = y - x`)."
const ABSOLUTE = Absolute()

(::Absolute)(x, Δ) = x + Δ

difference(::Absolute, y, x) = y - x

####
#### perturbations
####

"""
A container for a single perturbation. Internal, use [`perturbation`](@ref) to
create, which documents the fields.
"""
struct Perturbation
    label::AbstractString
    metadata
    lens
    change
    captured_errors
    Δs
end

function Base.show(io::IO, perturbation::Perturbation)
    @unpack label, change, captured_errors, Δs = perturbation
    domain = Δs ≡ nothing ? "default domain" : extrema(Δs)
    print(io, "perturb “$(label)” on $(domain), $(change) [capturing $(captured_errors)]")
end

"""
$(SIGNATURES)

# Arguments

- `label`: Label, used for lookup and printing.

- `lens`: For changing an object via the `Accessors.lens` API.

- `change`: A change calculator, eg [`ABSOLUTE`](@ref) or [`RELATIVE`](@ref)`, or a
  callable with the signature `(x, Δ) -> new_x`.

# Keyword arguments

- `metadata`: arbitrary metadata (eg for plot styles, etc). Passed through.

- `Δs`: an `AbstractVector` or iterable of `Δ` values to be used by `change` above. When
  `nothing`, a default will be used, see [`perturbation_analysis`](@ref).

- `captured_errors`: a type (eg a `Union{DomainError,ArgumentError}`) that is silently
  captured and converted to `nothing`. These results will show up as `NaN` in moment
  sensitivities.
"""
function perturbation(label::AbstractString, lens, change; metadata = nothing, Δs = nothing, captured_errors = DomainError)
    Perturbation(label, metadata, lens, change, captured_errors, Δs)
end

"""
$(SIGNATURES)

Internal function for computing a perturbation of `object` by `Δ`. Captures errors and
converts them to `nothing`.
"""
function _perturbation_with(object, perturbation::Perturbation, Δ)
    @unpack lens, change, captured_errors = perturbation
    try
        Accessors.modify(x -> change(x, Δ), object, lens)
    catch e
        if e isa captured_errors
            nothing
        else
            rethrow(e)
        end
    end
end

"""
A container for perturbation analysis. Internal, use [`perturbation_analysis`](@ref) to create.
"""
struct PerturbationAnalysis
    label::Union{Nothing,AbstractString}
    object
    simulate
    baseline_result
    default_Δs
    perturbations_and_results
end

"""
$(SIGNATURES)

# Arguments

- `object`: an opaque object that will be modified by perturbations

- `simulate`: a callable on (changed) `object`s that returns a value moments can be
  calculated from. Should be thread safe. Can return `nothing` in case the simulation fails
  for whatever reason, in which case moments will use a `NaN` for the corresponding change.

# Keyword arguments

- `baseline_result`: a simulation result from `object` as is.

- `label`: a global label for the whole analysis.

- `default_Δs`: when a perturbation has no `Δs` specified, these will be used instead.
  Defaults to 11 values on ±0.1.

# Usage

Once an analysis is created, you should `push!` perturbations into it. `label`s should be
unique, as they can be used for lookup.

```jldoctest
julia> using SensitivityAnalysis, Accessors

julia> anls = perturbation_analysis((a = 10.0, ), x -> (b = 2 * x.a, ))
perturbation results with default domain (-0.1, 0.1)

julia> push!(anls, perturbation("a", @optic(_.a), RELATIVE))
perturbation results with default domain (-0.1, 0.1)
  perturb “a” on default domain, relative change [capturing DomainError]

julia> moment_sensitivity(anls, "a", moment("b", x -> x.b, ABSOLUTE))
(label = "b (absolute change)", x = -0.1:0.02:0.1, y = [-2.0, -1.5999999999999979, -1.2000000000000028, -0.8000000000000007, -0.3999999999999986, 0.0, 0.3999999999999986, 0.8000000000000007, 1.2000000000000028, 1.6000000000000014, 2.0])
```

See also [`moment_sensitivity`](@ref).
"""
function perturbation_analysis(object, simulate; baseline_result = simulate(object),
                              label = nothing, default_Δs = range(-0.1, 0.1; length = 11))
    PerturbationAnalysis(label, object, simulate, baseline_result, default_Δs, Vector())
end

function Base.show(io::IO, perturbation_analysis::PerturbationAnalysis)
    @unpack label, default_Δs, perturbations_and_results = perturbation_analysis
    print(io, "perturbation results")
    if label ≢ nothing
        print(io, "for ", label)
    end
    print(io, " with default domain $(extrema(default_Δs))")
    for (p, _) in perturbations_and_results
        print(io, "\n  ", p)
    end
end

Broadcast.broadcastable(x::PerturbationAnalysis) = Ref(x)

"""
$(SIGNATURES)

Internal function for calculating `Δs`.
"""
function _effective_Δs(perturbation_analysis, perturbation)
    something(perturbation_analysis.default_Δs, perturbation.Δs)
end

function Base.push!(perturbation_analysis::PerturbationAnalysis, perturbation::Perturbation)
    @unpack object, simulate, perturbations_and_results = perturbation_analysis

    results = ThreadPools.tmap(_effective_Δs(perturbation_analysis, perturbation)) do Δ
        changed_object = _perturbation_with(object, perturbation, Δ)
        changed_object ≡ nothing && return nothing
        simulate(changed_object)
    end
    push!(perturbations_and_results, perturbation => results)
    perturbation_analysis
end

"""
$(SIGNATURES)

Internal function to look up perturbations by matching `pattern` on their labels. Checks for
unique matches.
"""
function _lookup_perturbation_by_label(perturbation_analysis,
                                       pattern::Union{AbstractPattern,AbstractString})
    @unpack perturbations_and_results = perturbation_analysis
    matches = findall(perturbations_and_results) do (p, _)
        contains(p.label, pattern)
    end
    if length(matches) > 1
        throw(ArgumentError("multiple labels match $(pattern)"))
    elseif isempty(matches)
        throw(ArgumentError("no labels match $(pattern)"))
    else
        first(matches)
    end
end

_lookup_perturbation_by_label(perturbation_analysis, index::Integer) = index

####
#### moments
####

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

"""
$(SIGNATURES)

Define a moment (a univariate statistic from simulation results).

# Arguments

- `label`: used for printing

- `calculator`: a callable on simulation results

- `change`: how changes should be interpreted, eg [`ABSOLUTE`](@ref)` or [`RELATIVE`](@ref)`.

# Keywoard arguments

- `metadata`: passed through unchanged, can be used for plot styles, etc.
"""
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
            difference(change, calculator(r), m0)
        end
    end
end

"""
$(SIGNATURES)

Return sensitivity to a particular moment.

# Arguments

- `perturbation_analysis`: a perturbation analysis, started by [`perturbation_analysis`](@ref)

- `index_or_pattern`: an integer that for a moment, or a pattern (eg a regex) that will be matched with labels

- `moment`: the moment to calculate, see [`moment`](@ref)
"""
function moment_sensitivity(perturbation_analysis, index_or_pattern, moment)
    @unpack baseline_result, perturbations_and_results = perturbation_analysis
    index = _lookup_perturbation_by_label(perturbation_analysis, index_or_pattern)
    perturbation, results = perturbations_and_results[index]
    sensitivity = _moment_sensitivity(moment, baseline_result, results)
    Δs = _effective_Δs(perturbation_analysis, perturbation)
    (label = repr(moment), x = Δs, y = sensitivity, metadata = moment.metadata)
end

end # module
