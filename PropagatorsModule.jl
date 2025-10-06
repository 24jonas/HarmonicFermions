# PropagatorsModule.jl

module PropagatorsModule

using LinearAlgebra

# Export the propagators constant so it can be used by other files.
export propagators

const propagators = (
    # ... (paste the exact same propagators definition from Method 1 here) ...
    PA = (
        zeta_1 = e -> 1 + e^2/2,
        lambda = e -> BigFloat(1.0),
        gamma = e -> sqrt(1 + e^2/4),
        k1 = e -> e
    ),
    TI = (
        zeta_1 = e -> begin
            omega_eff_sq = 1 + e^2 / 12
            omega_eff = sqrt(omega_eff_sq)
            return cosh(omega_eff * e)
        end,
        lambda = e -> BigFloat(1.0),
        gamma = e -> begin
            omega_eff_sq = 1 + e^2 / 12
            omega_eff = sqrt(omega_eff_sq)
            return iszero(e) ? BigFloat(1.0) : sqrt(sinh(omega_eff * e) / (omega_eff * e))
        end
    ),
    FA = (
        zeta_1 = e -> (1 + e^2/6) * cosh(e) - e^2/6,
        lambda = e -> 1 - e^2/12,
        gamma = e -> iszero(e) ? BigFloat(1.0) : sqrt( ((1+e^2/6)*sinh(e) - (e^3/6)*cosh(e)) / e )
    )
)

end # end of module