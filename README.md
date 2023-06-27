# BT-ITSP
This repository includes scripts used for pre-processing and analysis of sensory processing data as collected through the Infant-Toddler Sensory Profile (ITSP) in BATSS
as performed for the study: Individual variability in sensory processing: etiological structure and developmental links to autism (Bussu, Portugal, Falck-Ytter). 

The dataset consists of 622 same-sex twins tested at 5 months of age, and followed-up longitudinally until 36 months of age, when autistic traits were tested through the Quantitative Checklist for Autism in Toddlers, QCHAT questionnaire.

This repository includes a script for data preparation for the twin analysis and phenotypic GEE association (BT-ITSP_prepdata.R), a script for running univariate twin modelling on the ITSP scores (BT-ITSP_univariate.R), and a script for running multivariate twin modelling on the ITSP scores (BT-ITSP_multivariate.R). It also includes a script for running a factor analysis on the ITSP items, used to test structural validity of the instrument (BT-ITSP_factors.R).
Finally, it includes a script for running bivariate twin modelling with Cholesky decomposition, testing the etiological associations longitudinally between early sensory processing in infancy and later autistic traits (BT-ITSP_bivariate.R).

Of note, scripts require adaptation for file paths and for variable selection in relation to specific domain analysis (i.e., running on ITSP quadrant scores vs. section scores vs. factors).
