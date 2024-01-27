# Models to estimate genetic gain of soybean seed yield from annual multi‑environment feld trials

Krause, M.D., Piepho, HP., Dias, K.O.G. et al. Theor Appl Genet 136, 252 (2023). https://doi.org/10.1007/s00122-023-04470-3

This repositoty contains the following files:

1. dat.rds: a list with three elements: MET, parents, and G_matrix (see below).
2. RGG_Table6.R: codes of the models presented in Table 6, except for Model E9. Codes for the linearity metrics were also made available.

Future updates of this repositoty will include the code for Model E9, models presented in Table 7, and the simulator.

The elements saved in the data.rds object are the following:

**MET:** simulated data from simulation model B2 (s = 45) <br />
1. G            : genotypes <br />
2. L            : locations <br />
3. Y            : years (factor) <br />
4. Y_num        : years (numeric) <br />
5. trial        : trials (PYT and URT) <br />
6. mappingL     : covariate mapping ($z_{jk}$) <br />
7. first_Y_trial: first year of trial ($r_i$) <br />
8. parent       : binary variable stating if the $i^{th}$ genotype was used as a parent in crossing blocks <br />
9. check        : binary variable stating if the $i^{th}$ genotype was a check variety in MET <br />
10. P1           : parent 1 of the $i^{th}$ genotype (assuming biparental crosses) <br />
11. P2           : parent 2 of the $i^{th}$ genotype <br />
12. eBLUE        : Empirical Best Linear Unbiased Estimates of genotype means from single-trial analysis (Model C3) <br />
13. SE           : standard error of estimated eBLUE <br />
14. weight       : weights computed as $\frac{1}{SE^2}$ <br />
15. sim_gv       : true simulated genetic value <br />

**parents:** data from simulated crossing blocks <br />
1. G                    : genotypes used as parents <br />
2. number_crosses       : the number of hybridizations the $i^{th}$ genotype contributed to 
3. year_crossing_nursery: calendar year
4. sim_gv               : true simulated genetic value <br />

**G_matrix:** additive genomic relationship matrix computed with rrBLUP::A.mat() from simulated SNP dosages

