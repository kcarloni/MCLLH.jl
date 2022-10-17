
function LLHgeneral(k, weight_sum, weight_sq_sum; llh_funct)

    # handle ill-formed likelihoods
    if (weight_sum <= 0) | (weight_sq_sum < 0)
        (k == 0) && (return 0)
        return -Inf
    end
    μ, σsq = weight_sum, weight_sq_sum

    # in the limiting case <w²> = 0, return the poisson likelihood
    (σsq == 0) && (return logpmf_Poisson(k; λ=μ))

    return llh_funct(k)
end

"""
    LLHeff(k, weight_sum, weight_sq_sum)
    LLHeff(k, weights)

Computes the log of the L_Eff likelihood.

L_eff is the poisson likelihood, using a poisson distribution with
rescaled rate parameter to describe the Monte Carlo expectation, and
assuming a uniform prior on the rate parameter of the Monte Carlo.
This is the main result of the paper arXiv:1901.04645
"""
function LLHeff(k, weight_sum, weight_sq_sum; a=1, b=0)
    μ, σsq = weight_sum, weight_sq_sum
    α = μ^2/σsq + a
    β = μ/σsq
    llh_funct(k) = llh_gammaPriorPoisson(k; α, β)
    return LLHgeneral(k, μ, σsq; llh_funct)
end
LLHeff(k, weights; a=1, b=0) = LLHeff(k, sum(w), sum(@. w^2); a, b)


"""
    LLHmean(k, weight_sum, weight_sq_sum)
    LLHmean(k, weights)

Computes the log of the L_mean likelihood. 
L_mean is the poisson likelihood with gamma distribution prior, 
where the mean and variance are fixed to that of the weight distribution.
"""
LLHmean(k, weight_sum, weight_sq_sum) = LLHeff(k, weight_sum, weight_sq_sum; a=0, b=0)
LLHmean(k, weights) = LLHeff(k, weights; a=0, b=0)

"""
    LLHmode(k, weight_sum, weight_sq_sum)
    LLHmode(k, weights)

Computes the log of the L_mode likelihood. 
L_mode is the poisson likelihood with gamma distribution prior, 
where the mode and variance are fixed to that of the weight distribution.
"""
function LLHmode(k, weight_sum, weight_sq_sum)
    μ, σsq = weight_sum, weight_sq_sum
    s = sqrt(1 + 4σsq/μ^2)
    α = μ^2/σsq * 1/2 * (1 + s) + 1
    β = μ/σsq * 1/2 * (1 + s) 
    llh_funct(k) = llh_gammaPriorPoisson(k; α, β)
    return LLHgeneral(k, μ, σsq; llh_funct)
end 
LLHmode(k, weights) = LLHmode(k, sum(w), sum(@. w^2))