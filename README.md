# ImageJ

Find the paper published in [L&O methods](https://doi.org/10.1002/lom3.10569) 

Code and data for determining ImageJ accuracy to estimate hypsography from bathymetric maps.

  
Hypso_QC.R has the code for quality control of digitizer hypsographies. It takes in a directory of proportion hypsography stored as csv files and outputs which csv files had common errors (such as a comma instead of a decimal point) and what the errors are. There are potentially other but these cover the most common.
  
ImageJ_comparison.R has the code for changing the data into useable formats and calulating mean absolute difference, mean signed difference and other similar measures. It takes in the digitizer files (stored in Fourth_set_AVP, Fourth_set_CR) and the hypsography generated from Digital Elevation Models (stored in MN_bathy.rds) and outputs a file with the validation hypsographys (stored as all_validation_hypsos.csv), the maximum depth of each hypsography for each method (stored as max_depth_validation.csv).
