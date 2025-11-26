# RiskCompensationNonTargetedHPV

Replication of the simulation studies in "Using NonTargeted HPV Infections in Observational Studies with Risk Compensation", by Etievant (2025). The estimation methods used in the article were previously proposed by Etievant et al. (Biometrics, 2023) and are available in file `EstimationFunctions.R`.

## Required packages 

```
ggplot2, parallel, openxlsx and xtable.
```

## Scripts

* Script `Simul_Observational.R` replicates the simulation proposed by Etievant in Section 3 in the Main Document.

* Script `Simul_RandomizedUnblinded.R` replicates the simulation proposed in Web Appendix I.

#### Remarks

* Each script relies on functions provided in `EstimationFunctions.R` and certain parameters provided in `Parameters.RData`.
  
* To replicate the additional simulation proposed in Web Appendix E.3, use script `Simul_Observational.R` with `zeta_T = -1` instead of `zeta_T = 1`.

* R code in scripts `Simul_Observational.R` and `Simul_RandomizedUnblinded.R` have been written to run the different scenarios in parallel and save computation time. Change the value of the `mc.cores` argument in the `mclapply` command if needed.


## Instructions to run each script

* Save the chosen script(s) and files `EstimationFunctions.R` and `Parameters.RData` in the same directory.

* Open and run the script(s).

* The results of the simulations are saved in .RData files, .pdf figures, and .csv tables. For example, when running script `Simul_Observational.R`, file `Comparison_Observational-beta1-0.73.pdf` will give Figure 4 in Section 4 in the Main Document. Table 1 in Section 4 in the Main Document can be directly obtained via the LaTeX table printed in the Console.


## Functions provided in `EstimationFunctions.R`

* **JointNC** - estimation with method Joint-NC. Relies on the joint estimation of the treatment effect on the targeted (i.e., primary) and untargeted (i.e., secondary) outcomes. Information on potential observed covariates or confounders is not used

* **JointMH** - estimation method Joint-MH. Relies on the joint estimation of the treatment effect on the targeted (i.e., primary) and untargeted (i.e., secondary) outcomes, using stratification on observed confounders with Mantel-Haenszel weights. 

* **JointReg** - estimation with method Joint-Reg. Relies on joint estimation of the treatment effect on the targeted (i.e., primary) and untargeted (i.e., secondary) outcomes, using regression models where observed confounders are included. For continuous confounders, quadratic polynomial functions are used.

* **JointRegCat** and **JointRegCont** - estimation with method Joint-Reg, just as in function **JointReg**, but when the observed confounders are all categorical or all continuous.
