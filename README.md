# Models to estimate genetic gain of soybean seed yield from annual multi‑environment feld trials

Krause, M.D., Piepho, HP., Dias, K.O.G. et al. Theor Appl Genet 136, 252 (2023). https://doi.org/10.1007/s00122-023-04470-3

This repositoty contains the following files:

1. dat.rds: a list with three elements: MET, parents, and G_matrix (see below).
2. RGG_Table6.R: codes of the models presented in Table 6, except for Model E9.

Future updates of this repositoty will include the code for Model E9, models presented in Table 7, and the simulator.

The elements saved in the data.rds object are the following:

**MET:** simulated data from simulation model B2 (s = 45)
G            : genotypes
L            : locations
Y            : years (factor)
Y_num        : years (numeric) 
trial        : trials (PYT and URT)
mappingL     : covariate mapping ($z_{jk}$)
first_Y_trial: first year of trial ($r_i$)
parent       : binary variable stating if the $i^{th}$ genotype was used as a parent in crossing blocks
check        : binary variable stating if the $i^{th}$ genotype was a check variety in MET
P1           : parent 1 of the $i^{th}$ genotype (assuming biparental crosses)
P2           : parent 2 of the $i^{th}$ genotype
eBLUE        : Empirical Best Linear Unbiased Estimates of genotype means from single-trial analysis (Model C3)
SE           : standard error of estimated eBLUE
weight       : weights computed as $\frac{1}{SE^2}$
sim_gv       : true simulated genetic value
