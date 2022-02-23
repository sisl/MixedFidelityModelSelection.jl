initialize_betas(results; α=1, β=1) = Dict(k=>Beta(α,β) for k in keys(results))


function update_betas(results; α=1, β=1)
    betas = initialize_betas(results; α=α, β=β)
    for k in keys(results_regret)
        betas[k] = update_beta(betas[k], results_regret[k])
    end
    betas
end


function update_beta(beta, results)
    decisions = results[:last_action]
    true_decisions = get_true_decisions(results)
    for (i,decision) in enumerate(decisions)
        r = (decision == true_decisions[i])
        beta = Beta(beta.α + r, beta.β + 1 - r)
    end
    return beta
end


function betavalue(key; bayesian=true)
    if bayesian
        # Bayesian
        return mean(betas[key])
    else
        # Frequentist
        return (betas[key].α - 1) / (betas[key].α + betas[key].β - 2)
    end
end
