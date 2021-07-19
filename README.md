# SparseCard

Code for approximate graph reduction techniques for cardinality-based DSFM, from paper 

"Augmented Sparsifiers for Generalized Hypergraph Cuts with Applications to Decomposable Submodular Function Minimization"

https://arxiv.org/abs/2007.08075

The `include` folder contains implementations for outside code needed for experimental comparisons.

The `src` folder includes implementations of our main methods.

There is a folder for each of the main experiments (one for image segmentation, one for hypergraph clustering).

### For image segmentation experiments

Code for running image segmentation experiments with competing continuous optimization techniques is given in 

`include/DSFM-with-incidence-relations-v2`

In order to reproduce these experiments, see 

`Run_all_image_exps.m`

### For hypergraph clustering experiments

You will need to place the stackoverflow-answers dataset in the `data` folder in order to run experiments

[https://www.cs.cornell.edu/~arb/data/stackoverflow-answers/](https://www.cs.cornell.edu/~arb/data/stackoverflow-answers/)

To reproduce experiments, run

`stackoverflow_runall.jl`

This will take a long time, as this file runs 4500 individual local clustering experiments.
