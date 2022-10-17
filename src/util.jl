
using SpecialFunctions: loggamma

# utility functions ==============================================

"""
    llh_gammaPriorPoisson(k; α, β)

Poisson distribution marginalized over the rate parameter, 
priored with a gamma distribution that has 
shape parameter α and inverse rate parameter β. 
Takes the number of observed events k and returns the log-likelihood.

"""
function llh_gammaPriorPoisson(k; α, β)
    values = [
        α * log(β), 
        -(k+α)*log1p(β),
        -loggamma(k+1), 
        loggamma(k+α),
        -loggamma(α)
    ]    
    return sum(values)
end

"""
    logpmf_Poisson(k; λ)

Return the logPMF of the Poisson distribution Pois(k, λ).
"""
logpmf_Poisson(k; λ) = k*log(λ) - loggamma(k+1) - λ