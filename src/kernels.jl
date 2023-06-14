"""
    KernelFunction

Abstract supertype to represent exponential relaxation/recovery functions for use in ILT.

Extended by singleton subtypes, each of which defines the functional form it represents as a functor. The current set of subtypes are:
- `InversionRecovery`, aliased as `t1ir`
- `SaturationRecovery`, aliased as `t1sr`
- `TransverseRelaxation`, aliased as `t2`
As these symbols are not exported, they must be referenced as `ILT.t1ir()` etc.
"""
abstract type KernelFunction end

struct InversionRecovery <: KernelFunction end
struct SaturationRecovery <: KernelFunction end
struct TransverseRelaxation <: KernelFunction end

const t1ir = InversionRecovery
const t1sr = SaturationRecovery
const t2 = TransverseRelaxation

function (k::t1ir)(t::Real, s::Real)
    1.0 - 2.0*exp(-t*s)
end

function (k::t1sr)(t::Real, s::Real)
    1.0 - exp(-t*s)
end

function (k::t2)(t::Real, s::Real)
    exp(-t*s)
end
