# PropagatorsModule.jl

module PropagatorsModule

using LinearAlgebra

# Export the propagators constant so it can be used by other files.
export propagators

const propagators = (
    # ... (paste the exact same propagators definition from Method 1 here) ...
    PA = (
        zeta_1 = e -> 1 + e^2/2,
        d_zeta_1 = e-> e,
        lambda = e -> BigFloat(1.0),
        gamma = e -> sqrt(1 + e^2/4),
        k1 = e -> e
    ),
    TI = (
        zeta_1 = e -> 1 + e^2/2 + e^4/24,
        d_zeta_1 = e -> e + e^3/6,
        lambda = e -> 1 + e^2/6, 
        gamma = e -> sqrt(1 + e^2/3 + e^4/24 + e^6/576),
        k1 = e -> e
    ),
    FA = (
        zeta_1 = e -> 1 + e^2/2 + e^4/24 + e^6/864,
        d_zeta_1 = e -> e + e^3/6 + e^5/144,
        lambda = e -> BigFloat(1.0),
        gamma = e -> sqrt(1+e^4/(432+36*e^2)),
        k1 = e -> e*(1+(e^2)/12)^2
    )
)

end # end of module