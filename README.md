# RiskCompensationNonTargetedHPV

Replication of the simulation studies in "Using NonTargeted HPV Infections in Observational Studies with Risk Compensation", by Etievant (2025). The estimation methods used in the article were previously proposed by Etievant et al. (Biometrics, 2023) and are available in file `EstimationFunctions.R`.

### Required packages 

```
gee, ggplot2, grid, gtable, parallel, openxlsx and xtable.
```

### Scripts

* Script `Simul_Observational.R` replicates the simulations proposed by Etievant in Section 4 in the Main Document.

* Script `Simul_RandomizedUnblinded.R` replicates the simulations proposed in Web Appendix E. 

Each script relies on functions provided in `EstimationFunctions.R` and certain parameters provided in `Parameters.RData`.


### Instructions to run each script

* Save the chosen script(s) and files `EstimationFunctions.R` and `Parameters.RData` in the same directory.

* Open and run the script(s).

* The results of the simulations are saved in figures and tables. For example, when running script `Simul_Observational.R`, file `Comparison_Observational-beta1-0.73.pdf` will give Figure 4 in Section 4 in the Main Document.


### Functions provided in `EstimationFunctions.R`

* **JointNC** - estimation with method Joint-NC. Relies on the joint estimation of the treatment effect on the targeted (i.e., primary) and untargeted (i.e., secondary) outcomes. Information on potential observed covariates or confounders is not used

* **JointMH** - estimation method Joint-MH. Relies on the joint estimation of the treatment effect on the targeted (i.e., primary) and untargeted (i.e., secondary) outcomes, using stratification on observed confounders with Mantel-Haenszel weights. 

* **JointReg** - estimation with method Joint-Reg. Relies on joint estimation of the treatment effect on the targeted (i.e., primary) and untargeted (i.e., secondary) outcomes, using regression models where observed confounders are included. For continuous confounders, quadratic polynomial functions are used.

* **JointRegCat** and **JointRegCont** - estimation with method Joint-Reg, just as in function **JointReg**, but when the observed confounders are all categorical or all continuous.
