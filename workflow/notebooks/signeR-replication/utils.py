#!/usr/bin/env python
import aesara
import aesara.tensor as at
import arviz as az
import xarray as xr
import numpy as np
from collections import OrderedDict

srng = at.random.RandomStream()

def set_srng(seed):
    global srng
    srng = at.random.RandomStream(seed)

def update_params(exposure, signature, alpha_exposure, alpha_signature, counts):
    # probabilidade de uma mutação da amostra i pertencer à assinatura k E classe j (conjunta)
    joint = signature[None,:,:] * exposure[:,:,None]

    # probabilidade de uma mutação da amostra i E classe j pertencer à assinatura k
    conditional = joint.swapaxes(1,2)/joint.sum(axis=1)[:,:,None]

    # Gibbs updates
    Z = srng.multinomial(counts, conditional)
    new_exposure = srng.dirichlet(alpha_exposure + Z.sum(axis=1))
    new_signature = srng.dirichlet(alpha_signature + Z.sum(axis=0).T)

    return new_exposure, new_signature

def update_hyperparams_Dirichlet(alpha, proportions, prior_lambda):
    # Metropolis proposal
    metropolis_sigma = 10
    a, b = (alpha/metropolis_sigma)**2, alpha/metropolis_sigma**2
    alpha_proposal = srng.gamma(a, b)

    # evaluating priors and conditionals
    dirichlet_logp = at.sum((alpha_proposal-1)*at.log(proportions)) + at.gammaln(alpha_proposal.sum(axis=-1)).sum() - at.gammaln(alpha_proposal).sum()
    prev_dirichlet_logp = at.sum((alpha-1)*at.log(proportions)) + at.gammaln(alpha.sum(axis=-1)).sum() - at.gammaln(alpha).sum()
    exponential_logp = - prior_lambda * alpha_proposal.sum()
    prev_exponential_logp = - prior_lambda * alpha.sum()

    # Hastings adjustment
    a_rev, b_rev = (alpha_proposal/metropolis_sigma)**2, alpha_proposal/metropolis_sigma**2
    forward_logp = a*at.log(b) - at.gammaln(a) + (a-1)*at.log(alpha_proposal) - b*alpha_proposal
    backward_logp = a_rev*at.log(b_rev) - at.gammaln(a_rev) + (a_rev-1)*at.log(alpha) - b_rev*alpha

    log_acceptance_prob = (
        (dirichlet_logp - prev_dirichlet_logp) + 
        (exponential_logp - prev_exponential_logp) +
        (forward_logp - backward_logp)
        )
    return at.switch(at.lt(srng.uniform(), np.exp(log_acceptance_prob)), alpha_proposal, alpha)

def update_Exponential_hyperparams(alpha):
    return srng.gamma(alpha.shape.prod()+1, alpha.sum()+1).astype("float32")

def Sampler(counts_df, step, mcmc_steps=500, k=5):
    """
    Parameters
    ----------
    counts_df : pd.DataFrame
        DataFrame with the counts for each mutation type (columns) in each 
        sample (rows).
    step : callable
        Function receiving the parameters from the last iteration and returning
        an updated version for the next iteration.
    mcmc_steps : int
        Number of MCMC iterations to be performed.
    k : int
        Number of signatures in this analysis.
    """
    n, p = counts_df.shape
    counts = at.as_tensor(counts_df.values)

    # Initialization
    tensors = OrderedDict(
        exposure = at.ones((n,k))/k,
        signature = at.ones((k,p))/p,
        alpha_exposure = at.ones((n,k))/2,
        alpha_signature = at.ones((k,p))/2,
        lambda_exposure = at.as_tensor_variable(1.0),
        lambda_signature = at.as_tensor_variable(1.0)
    )

    # Loop definition
    outputs, updates = aesara.scan(
        step,
        outputs_info=list(tensors.values()),
        non_sequences=[counts],
        n_steps=mcmc_steps,
        strict=True,
    )
    sampler = aesara.function(inputs=[], outputs=outputs, updates=updates)
    mcmc_results = OrderedDict(zip(tensors.keys(), [np.expand_dims(var, 0) for var in sampler()]))

    # Numpy -> Arviz
    coords = dict(
        draw = ("draw", range(mcmc_steps)),
        specimen = ("specimen", counts_df.index),
        sig = ("sig", [f"S{num}" for num in range(1,k+1)]),
        substitution = ("substitution", counts_df.columns),
        chain = ("chain", [1]))

    dims = OrderedDict(
        exposure = ["specimen","sig"],
        signature = ["sig","substitution"],
        alpha_exposure = ["specimen","sig"],
        alpha_signature = ["sig","substitution"],
        lambda_exposure = [],
        lambda_signature = []
    )

    idata = az.InferenceData(
        posterior = xr.Dataset(
            data_vars = {var_name:(["chain","draw"]+var_coords, mcmc_results[var_name]) for var_name, var_coords in dims.items()},
            coords = coords))

    return idata

def sample_full_model(counts_df, mcmc_steps=500, k=5):
    def step(exposure, signature, alpha_exposure, alpha_signature, lambda_exposure, lambda_signature, counts):
        exposure, signature = update_params(exposure, signature, alpha_exposure, alpha_signature, counts)
        alpha_exposure = update_hyperparams_Dirichlet(alpha_exposure, exposure, lambda_exposure)
        alpha_signature = update_hyperparams_Dirichlet(alpha_signature, signature, lambda_signature)
        lambda_exposure = update_Exponential_hyperparams(alpha_exposure)
        lambda_signature = update_Exponential_hyperparams(alpha_signature)
        return exposure, signature, alpha_exposure, alpha_signature, lambda_exposure, lambda_signature
    
    return Sampler(counts_df, step, mcmc_steps, k)

def sample_fixed_hyperparams(counts_df, mcmc_steps=500, k=5):
    def step(exposure, signature, alpha_exposure, alpha_signature, lambda_exposure, lambda_signature, counts):
        exposure, signature = update_params(exposure, signature, alpha_exposure, alpha_signature, counts)
        return exposure, signature, alpha_exposure, alpha_signature, lambda_exposure, lambda_signature
    
    return Sampler(counts_df, step, mcmc_steps, k)