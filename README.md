# Models to estimate genetic gain of soybean seed yield from annual multi‑environment field trials

Krause, M.D., Piepho, HP., Dias, K.O.G. et al. Theor Appl Genet 136, 252 (2023). https://doi.org/10.1007/s00122-023-04470-3

This repository contains the following files:

1. dat_seed_45.rds: a list with four elements: MET, parents, G_table6 (G matrix for models in Table 6), and G_table7 (G matrix for models in Table 7). See details below.
2. RGG_Table6.R: models in Table 6, except E9, and linearity metrics.
3. RGG_Table6_E9.R: model E9 from Table 6.
4. RGG_Table7.R: models in Table 7.

The dat_seed_45.rds object contains the following datasets:

**MET:** simulated data from simulation scenario/model B2 (s = 45) <br />
1. G            : genotypes <br />
2. L            : locations <br />
3. Y            : years (factor) <br />
4. Y_num        : years (numeric) <br />
5. trial        : trials (BT3, PYT and URT) <br />
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

**G_table6 and G_table7:** additive genomic relationship matrices computed with rrBLUP::A.mat() from simulated SNP dosages

