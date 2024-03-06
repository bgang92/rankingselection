This folder contains the code for the paper "Ranking and Selection in Large-Scale Inference of Heteroscedastic Units" by Fu et al.

Requirement:
platform       aarch64-apple-darwin20      
arch           aarch64                     
os             darwin20                    
system         aarch64, darwin20           
status                                     
major          4                           
minor          1.0                         
year           2021                        
month          05                          
day            18                          
svn rev        80317                       
language       R                           
version.string R version 4.1.0 (2021-05-18)
CVXR package version 1.0.9

Code abstract:
To address the challenge of over-representations of subpopulations when ranking and selecting top candidates from heteroscedastic units,
we propose a new multiple comparison framework that incorporates a modified power notion to prioritize the selection of important effects
and employs a novel ranking metric to assess the relative importance of units.

Code overview:
3M.func.R contains codes that implements the BH procedure and the SC procedure (Sun and Cai (2007)). The two procedures are used as benchmarks for comparison.

datagen.R contains codes that used to generate data in the simulation sections.

nest.func.R contains codes for bivaraite density estimation and deconvolution.

rih.R contains the main functions. rih.func implements algorithm 1 in the paper. rvalue.func and rvalue2.func computes the two r-values discussed in section 4.

To replicate the results from say section 5.3 run "section5_3.R"

realdata_CRSP.R and realdata_AYP.R are the code used for analyzing the two read datasets. The data used can be found in mean_se_5yrs.csv and AYP_05.csv.

cb.colors.R contains helper functions used for generating plots.

