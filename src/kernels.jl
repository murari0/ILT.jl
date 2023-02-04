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