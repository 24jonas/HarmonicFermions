module EnergySol
using ArbNumerics

export harmEnergy, get_factor, get_energy_calc

# --- Optimized Iterative Core Logic ---

function harmEnergy(num_fermions::Int, dimensions::Int, b::ArbFloat)
    # 1. Pre-calculate the small 'z' terms iteratively
    # Using arrays prevents re-calculating powers of 'b' repeatedly
    z_vals = Vector{ArbFloat}(undef, num_fermions)
    z_prime_vals = Vector{ArbFloat}(undef, num_fermions)

    b_pow = b # Tracks b^i

    # We use a pure loop to populate z values. 
    # This avoids "pow" overhead and deep recursion.
    for i = 1:num_fermions
        # z_i formula
        one_minus_b_i = 1 - b_pow

        # Optimization: b^(i*d/2) calculation
        if dimensions == 2
            numerator = b_pow
        else
            numerator = b_pow^(dimensions/2)
        end

        z = numerator / (one_minus_b_i^dimensions)

        # z'_i formula
        # We need b^(i-1). Since b_pow is b^i, b^(i-1) = b_pow / b
        # For i=1, b^(0) is 1.0
        b_pow_minus_1 = (i == 1) ? one(ArbFloat) : (b_pow / b)

        log_deriv =
            (ArbFloat(i * dimensions) / (2 * b)) +
            (dimensions * i * b_pow_minus_1) / one_minus_b_i
        z_prime = z * log_deriv

        z_vals[i] = z
        z_prime_vals[i] = z_prime

        # Advance b_pow for next iteration (b^(i+1) = b^i * b)
        if i < num_fermions
            b_pow *= b
        end
    end

    # 2. Compute Big Z iteratively (Newton-Girard Formula)
    # Replaces the recursive function Z_recursive!
    # Z[k] stores Z_{k-1} (because Julia arrays are 1-indexed)
    Z = Vector{ArbFloat}(undef, num_fermions + 1)
    Z_prime = Vector{ArbFloat}(undef, num_fermions + 1)

    Z[1] = one(ArbFloat)      # Z_0
    Z_prime[1] = zero(ArbFloat) # Z'_0

    for k = 1:num_fermions
        sum_Z = zero(ArbFloat)
        sum_Z_prime = zero(ArbFloat)

        for i = 1:k
            sign = iseven(i - 1) ? 1 : -1

            # Z[k-i + 1] accesses Z_{k-i}
            term_Z = sign * z_vals[i] * Z[k-i+1]
            term_Z_prime = sign * (z_prime_vals[i] * Z[k-i+1] + z_vals[i] * Z_prime[k-i+1])

            sum_Z += term_Z
            sum_Z_prime += term_Z_prime
        end

        Z[k+1] = sum_Z / k
        Z_prime[k+1] = sum_Z_prime / k
    end

    # 3. Final Calculation
    Z_final = Z[num_fermions+1]
    Z_prime_final = Z_prime[num_fermions+1]

    # Safety check for log
    if isnan(Z_prime_final) || Z_final <= 0
        return ArbFloat(NaN)
    else
        log_Z_prime = Z_prime_final / Z_final
        return Float64(b * log_Z_prime)
    end
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
