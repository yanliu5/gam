# GAM

This repository contains R programs for the article, “Generalized additive regression for group testing data,” by Yan Liu, Christopher McMahan, Joshua Tebbs, Colin Gallagher and Christopher Bilder. This article was accepted by *Biostatistics*.

<pre>
1. <b>main_error_known.R</b>: Reproduces <b>Table 1</b>
2. <b>main_error_unknown.R</b>: Reproduces <b>Table 2</b>
3. <b>Simulated_dataset.csv</b>: Simulated data set that resembles the data set analyzed in Section 4 of the manuscript.
   It contains testing responses from Dorfman testing and covariates from individuals.
4. <b>GP_uninform_prior.R</b>: Fits a model to <b>Simulated_dataset.csv</b> using Gaussian Process (GP) method.
5. <b>GPP_uninform_prior.R</b>: Fits a model to <b>Simulated_dataset.csv</b> using Gaussian Predictive Process (GPP) method.
<pre>
