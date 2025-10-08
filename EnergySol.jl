module EnergySol

# Make the main recursive function available for import.
export harmEnergy, get_factor, get_energy_calc


function z_and_derivative_high_precision(i::Int, d::Int, b::BigFloat)
    if !(0 < b < 1)
        return (zero(BigFloat), zero(BigFloat))
    end
    
    b_i = b^i
    one_minus_b_i = 1 - b_i
    
    # Direct formula for z_i = b^(i*d/2) / (1 - b^i)^d
    z = b^(BigFloat(i * d) / 2) / (one_minus_b_i)^d
    
    # We still need the logarithmic derivative to find z_i'
    # (log(z))' = (i*d)/(2*b) + (d*i*b^(i-1))/(1-b^i)
    log_deriv = (BigFloat(i * d) / (2 * b)) + (d * i * b^(i - 1)) / one_minus_b_i
    
    # Direct formula for the derivative: z_i' = z_i * (log(z_i))'
    z_prime = z * log_deriv
    
    return z, z_prime
end

function Z_and_derivative_recursive!(n::Int, d::Int, b::BigFloat, cache::Dict)
    # The cache now stores (Z_n, Z_n') instead of their logs.
    if haskey(cache, n)
        return cache[n]
    end
    # Base case: Z_0 = 1, Z_0' = 0
    if n == 0
        return (one(BigFloat), zero(BigFloat))
    end

    total_Z = zero(BigFloat)
    total_Z_prime = zero(BigFloat)

    for i in 1:n
        # --- Get z_i and its derivative z_i' ---
        z_i, z_i_prime = z_and_derivative_high_precision(i, d, b)

        # --- Recursive call for the remaining n-i particles ---
        Z_ni, Z_ni_prime = Z_and_derivative_recursive!(n - i, d, b, cache)
        
        sign = iseven(i - 1) ? 1 : -1

        # --- Sum the terms directly ---
        # For Z_n
        total_Z += sign * z_i * Z_ni
        
        # For Z_n', using the product rule: (f*g)' = f'*g + f*g'
        total_Z_prime += sign * (z_i_prime * Z_ni + z_i * Z_ni_prime)
    end

    # Final result is (1/n) * sum
    final_Z = total_Z / n
    final_Z_prime = total_Z_prime / n
    
    cache[n] = (final_Z, final_Z_prime)
    return final_Z, final_Z_prime
end 

function harmEnergy(num_fermions::Int, dimensions::Int, b::BigFloat, cache::Dict)

    Z_val, Z_prime = Z_and_derivative_recursive!(
                num_fermions, dimensions, b, cache
            )

    # The logarithmic derivative is Z'/Z.
    energy = if isnan(Z_prime) || Z_val <= 0 # Avoid division by zero or log of negative
        NaN
    else
        log_Z_prime = Z_prime / Z_val # This is the log derivative
        # The rest of the formula is the same
        Float64((b * log_Z_prime))
    end

    return energy

end



function get_factor(mode::String)
    if mode == "thermo"
        # Return a new function specifically for the "thermo" calculation
        return (lambda, gamma, w, lambda_s, gamma_s) -> begin
            factor_regular = lambda / gamma
            factor_star = (w^2) * lambda_s / gamma_s
            return factor_regular, factor_star
        end
    else
        # Return a new function for the default calculation
        return (lambda, gamma, w, lambda_s, gamma_s) -> begin
            factor_regular = 0.5 * (gamma + 1 / gamma)
            factor_star = (w / 2) * (gamma_s / w + w / gamma_s)
            return factor_regular, factor_star
        end
    end
end



function get_energy_calc(balanced::Bool)
    if balanced
        # Return a function for the balanced calculation
        return (num_fermions, dimensions, b_s, cache2) -> 
            harmEnergy(floor(Int, num_fermions / 2), dimensions, b_s, cache2) + harmEnergy(ceil(Int, num_fermions / 2), dimensions, b_s, cache2)
    else
        # Return a function for the unbalanced calculation.
        # It still accepts cache2 for a consistent function signature, even if unused.
        return (num_fermions, dimensions, b_s, cache2) -> 
            harmEnergy(num_fermions, dimensions, b_s, cache2)
    end
end



end