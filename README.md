# GAM

This repository contains R programs for the article, “Generalized additive regression for group testing data,” by Yan Liu, Christopher McMahan, Joshua Tebbs, Colin Gallagher and Christopher Bilder. This article has been submitted for publication.

<pre>
1. To reproduce <b>Table 1</b> in the manuscript; i.e., analyze outcome from four testing procedures
   (Individual Test, Master Pool Test, Dorfman Test and Array Test), run <b>main_error_known.R</b>
2. To reproduce <b>Table 2</b> in the manuscript; i.e., analyze outcome from two testing procedures 
   (Dorfman Test and Array Test), run <b>main_error_unknown.R</b>
3. <b>Simulated_dataset.csv</b> is the simulated data set that very closely resembles the Iowa data 
   set analyzed in Section 4 of the manuscript.
4. <b>GP_uninform.R</b> fits <b>Simulated_dataset.csv</b> using Gaussian Process (GP) method.
5. <b>GPP_uninform.R</b> fits <b>Simulated_dataset.csv</b> using Gaussian Predictive Process (GPP) method.
<pre>
