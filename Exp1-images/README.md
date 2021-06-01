## README


There are four image segmentations DSFM problems: 2 pictures (octopus and smallplant) times 2 different superpixel segmentations (200 and 500). 

**image\_graph\_utils:** useful functions for processing the data and running some of the experiments.

**Solve\_Exact\_Flow:** solves the 4 DSFM problems optimally via exact reduction to a maximum flow. Takes a long time. For large superpixel regions, optimal result (used only to compare approximation factor achieved by other methods) was obtained by running this on a server with much more memory. Output is saved already.


## To reproduce experiments

In order to reproduce experiments, you will need to run the following:

**save\_image\_noclique**: saves unary potentials and pairwise potentials for each of the 4 experiments

**SparseCard_images** Runs SparseCard on all instances.

**plot\_[dataset]:** code for plotting results plotting for each instance of DSFM.

You will also need to run experiments using competing methods, in the `include/DSFM...` folder.