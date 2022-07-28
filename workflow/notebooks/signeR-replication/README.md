# signeR replication
Files in this directory document an attempt to replicating the signature inference algorithm implemented in `signeR`, with some modifications. The main modifications are using plain MCMC instead of EM-MCMC and using Dirichlet priors instead of Gamma.

## Original paper
```
Rosales, Rafael A., Rodrigo D. Drummond, Renan Valieris, Emmanuel Dias-Neto, e Israel T. Da Silva. “signeR: An empirical Bayesian approach to mutational signature discovery”. Bioinformatics 33, nº 1 (2017): 8–16. https://doi.org/10.1093/bioinformatics/btw572.
```

## To-do
Due to a coding mistake, the alpha parameters are generated for each sample (in case of exposures) or for each signature (in case of signature spectra). This means there is effectively no information-sharing between signatures/samples, and is probably why the original signeR signatures are cleaner than ours. This issue will be solved as soon as possible.