Demo for DSFM with incidence relations 

reference: Revisiting Decomposable Submodular Function Minimization with Incidence Relations, NIP2018

main.m : main file to run the demo, MAP inference on the image smallplant
RCDM_cversion.cpp : RCD
RCDM_cversion_greedypar.cpp: RCD while using greedy partition
ACDM_cversion.cpp : ACD
APcompact_cversion.cpp: IAP
APfull_cversion.cpp: AP 


Updates to Readme in May, 2021
------------------------------
For the Neurips 2021 submission on "Approximate Decomposable Submodular Function Minimization for Cardinality-Based Components", this software has slightly been edited in order to allow make possible a proper comparison with other approaches implemented in Julia.

The following changes have been made to each of the main implementations (RCDM, ACDM, etc.)

* Added input parameter gaptol, which terminates implemented algorithms early if a certain tolerance is reached
* Added timer for keeping track of how long each iteration takes
* Removed calculation of smooth tolerance
* Returned minimum level set cost for each iteration
* Returned runtime at each iteration
* Returned pos_x at each iteration (for computing discrete gap at each iteration)

Other slight updates have also been added in order to compare against algorithms run separately using Julia code. 

In order to run image segmentation experiments, see Run_all_image_exps.m.

For post-hoc tests on alpha and c parameters, see grid_ADCM_wrapper.m and grid_RCDM_wrapper.m