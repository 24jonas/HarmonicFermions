module EnergySol
using ArbNumerics

export harmEnergy, get_factor, get_energy_calc

# --- Optimized Iterative Core Logic ---

function degeneracy(m::Int, d::Int)
    if d == 1
        return 1
    end
    if d == 2
        return m + 1
    end
    if d == 3
        return div((m + 1) * (m + 2), 2)
    end
    return Int(binomial(BigInt(m + d - 1), BigInt(d - 1)))
end

function harmEnergy(num_fermions::Int, dimensions::Int, b::ArbFloat) ################# GRAND Z APPROACH
    # Z[n+1] stores the coefficient for xi^n
    Z = [ArbFloat(0) for _ = 1:(num_fermions+1)]
    Z_prime = [ArbFloat(0) for _ = 1:(num_fermions+1)]

    Z[1] = ArbFloat(1)

    m = 0
    capacity = 0

    while true
        g_m = degeneracy(m, dimensions)

        max_k = min(g_m, num_fermions)
        C = [ArbFloat(0) for _ = 1:(max_k+1)]
        C_prime = [ArbFloat(0) for _ = 1:(max_k+1)]

        C[1] = ArbFloat(1)
        b_to_m = b^m

        # Precompute the binomial coefficients for energy level m
        for k = 1:max_k
            factor = ArbFloat(g_m - k + 1) / ArbFloat(k)
            C[k+1] = C[k] * factor * b_to_m

            if b > 0
                C_prime[k+1] = C[k+1] * ArbFloat(m * k) / b
            elseif m * k == 1
                C_prime[k+1] = ArbFloat(g_m)
            else
                C_prime[k+1] = ArbFloat(0)
            end
        end

        # Convolve the polynomial
        Z_new = [ArbFloat(0) for _ = 1:(num_fermions+1)]
        Z_prime_new = [ArbFloat(0) for _ = 1:(num_fermions+1)]

        for n = 0:num_fermions
            for k = 0:min(n, max_k)
                Z_new[n+1] += C[k+1] * Z[n-k+1]
                Z_prime_new[n+1] += C_prime[k+1] * Z[n-k+1] + C[k+1] * Z_prime[n-k+1]
            end
        end

        capacity += g_m

        # TRUNCATION: Only allowed IF the Fermi sea is fully built
        if capacity >= num_fermions
            if Z[num_fermions+1] > 0
                # Break when the thermal addition is negligibly small
                diff = abs(Z_new[num_fermions+1] - Z[num_fermions+1])
                if diff / Z[num_fermions+1] < ArbFloat("1e-40")
                    Z = Z_new
                    Z_prime = Z_prime_new
                    break
                end
            end

            # Hard limit to prevent infinite loops at very high temperatures
            if m > num_fermions + 500
                Z = Z_new
                Z_prime = Z_prime_new
                break
            end
        end

        Z = Z_new
        Z_prime = Z_prime_new
        m += 1
    end

    Z_N = Z[num_fermions+1]
    Z_prime_N = Z_prime[num_fermions+1]

    if isnan(Z_prime_N) || Z_N <= 0
        return NaN
    end

    log_Z_prime = Z_prime_N / Z_N

    # Add back the analytical zero-point energy
    E_shift = (num_fermions * dimensions) / 2.0
    return Float64(E_shift + Float64(b * log_Z_prime))
end



function get_factor(mode::String)
    if mode == "thermo"
        return (lambda, gamma, w, lambda_s, gamma_s) -> begin
            factor_regular = lambda / gamma
            factor_star = (w) * lambda_s / gamma_s #No longer accounts for extra w from k1(e)
            return factor_regular, factor_star
        end
    else
        return (lambda, gamma, w, lambda_s, gamma_s) -> begin
            factor_regular = 0.5 * (gamma + 1 / gamma)
            factor_star = (w / 2) * (gamma_s + 1 / gamma_s) #No longer accounts for extra w from k1(e)
            return factor_regular, factor_star
        end
    end
end

function get_energy_calc(balanced::Bool)
    if balanced
        return (num_fermions, dimensions, b_s, ignored_cache) ->
            harmEnergy(floor(Int, num_fermions / 2), dimensions, b_s) +
            harmEnergy(ceil(Int, num_fermions / 2), dimensions, b_s)
    else
        return (num_fermions, dimensions, b_s, ignored_cache) ->
            harmEnergy(num_fermions, dimensions, b_s)
    end
end

end # module
