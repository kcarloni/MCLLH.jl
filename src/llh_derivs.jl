
# using SpecialFunctions: loggamma, digamma

# # for points at which the derivative is not defined? 
# const undef_val = nothing


# d_logpmf_Poisson_dk(k; λ) = log(λ) - digamma(k+1)
# d_logpmf_Poisson_dλ(k; λ) = k/λ - 1


# d_llh_gammaPriorPoisson_dk(k; α, β) = -log1p(β) - digamma(k) + digamma(k+α)
# d_llh_gammaPriorPoisson_dα(k; α, β) = log(β) - log1p(β) - digamma(α) + digamma(k+α)
# d_llh_gammaPriorPoisson_dβ(k; α, β) = α/β - (k+α)/(1+β)


# function d_LLHeff_dk(k, μ, σsq; a=1, b=0)
#     α = μ^2/σsq + a
#     β = μ/σsq + b
    
#     if (μ <= 0) || (σsq < 0); return (k==0) ? undef_val : 0.
#     elseif (σsq == 0); return d_logpmf_Poisson_dk(k; λ=μ)
#     else return d_llh_gammaPriorPoisson_dk(k; α, β)
#     end
# end
# function d_LLHeff_dμ(k, μ, σsq; a=1, b=0)
#     α = μ^2/σsq + a
#     β = μ/σsq + b
    
#     if (σsq < 0); return 0.
#     elseif (σsq == 0); return d_logpmf_Poisson_dλ(k; λ=μ)
#     elseif (μ < 0); return 0.
#     elseif (μ == 0); return undef_val
#     else
#         term_α = 2μ/σsq * d_llh_gammaPriorPoisson_dα(k; α, β)
#         term_β = 1/σsq * d_llh_gammaPriorPoisson_dβ(k; α, β)
#         return term_α + term_β
#     end
# end
# function d_LLHeff_dσsq(k, μ, σsq; a=1, b=0)
#     α = μ^2/σsq + a
#     β = μ/σsq + b
    
#     if (μ <= 0); return 0.
#     elseif (σsq <= 0); return 0.
#     else
#         term_α = -μ^2/σsq^2 * d_llh_gammaPriorPoisson_dα(k; α, β)
#         term_β = -μ/σsq^2 * d_llh_gammaPriorPoisson_dβ(k; α, β)
#         return term_α + term_β
# end


