### Code package: M. Fischer, N. Hauzenberger, F. Huber, & M. Pfarrhofer (2023). General Bayesian time-varying parameter vector autoregressions for modeling government bond yields. Journal of Applied Econometrics, 38(1), 69–87.


[**Publication (open access).**](https://doi.org/10.1002/jae.2936)


We use monthly zero-coupon yields of US treasuries at 30 different yearly maturities. In addition, we consider three observed factors governing the joint dynamics of the TVPs: (1) a binary recession indicator (labeled REC) for the US, (2) the National Financial Conditions Index (NFCI), and (3) the Risk-Free (RF) interest rate from the Fama-French Portfolios and Factors database. We provide all data as an .rda file [US_ylds](./data/US_ylds.rda).


In the paper, we consider two versions of our model. The first models the yields in an unrestricted manner. The second model class includes specifications based on the three-factor Nelson–Siegel (NS) model, which imposes a factor structure on the yields. Replication files replicate the latter in an equation-by-equation manner. To update estimates for the Nelson-Siegel factor model, one needs to source the file [main-file.R](main-file.R), which calls the main MCMC sampler [tvp_flex_func.R](./aux_funcs/tvp_flex_func.R) and another file [aux_files.R](./aux_funcs/aux_files.R), including some auxiliary functions (e.g., for the HS priors). In addition, [kf.cpp](./aux_funcs/kf.cpp) is a Kalman filtering function for the FFBS step of the RW factors.


