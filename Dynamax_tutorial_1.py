# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 13:48:10 2023

@author: neuro
"""
import jax.numpy as jnp
import jax.random as jr
import matplotlib.pyplot as plt
from dynamax.hidden_markov_model import GaussianHMM
from scipy import io
import pathlib as pl

import numpy as np
from itertools import count

fpath = pl.Path(r"Z:\Jesus\Jittering\FULLYCurated\05_190701_Jesus_Emili_Jittering_3800_1600_1500 yes VPM good\KS2_newChanMap\LFP_4_R.mat")
dat = io.loadmat(fpath)

y =  jnp.array( dat['lfp_z'] )
N = np.squeeze( dat['N'] )

keys = map(jr.PRNGKey, count())
num_states = 2
emission_dim = 1
num_timesteps = N

hmm = GaussianHMM(num_states, emission_dim)
params, props = hmm.initialize(next(keys), emission_means=(-3,2))
params, lps = hmm.fit_em(params, props, y )


# Make a Gaussian HMM and sample data from it
hmm = GaussianHMM(num_states, emission_dim)
true_params, _ = hmm.initialize(next(keys))
true_states, emissions = hmm.sample(true_params, next(keys), num_timesteps)

# Make a new Gaussian HMM and fit it with EM
params, props = hmm.initialize(key3, method="kmeans", emissions=emissions)
params, lls = hmm.fit_em(params, props, emissions, num_iters=20)

# Plot the marginal log probs across EM iterations
plt.plot(lls)
plt.xlabel("EM iterations")
plt.ylabel("marginal log prob.")

# Use fitted model for posterior inference
post = hmm.smoother(params, emissions)
print(post.smoothed_probs.shape) # (1000, 3)

from functools import partial
from jax import vmap

num_seq = 200
batch_true_states, batch_emissions = \
   vmap(partial(hmm.sample, true_params, num_timesteps=num_timesteps))(
      jr.split(key2, num_seq))
print(batch_true_states.shape, batch_emissions.shape) # (200,1000) and (200,1000,2)

# Make a new Gaussian HMM and fit it with EM
params, props = hmm.initialize(key3, method="kmeans", emissions=batch_emissions)
params, lls = hmm.fit_em(params, props, batch_emissions, num_iters=20)