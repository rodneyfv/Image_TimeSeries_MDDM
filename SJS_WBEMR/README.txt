TOOLBOX for wavelet estimation of mixture regression.

AUTHOR
Michel H. Montoril
Department of Statistics, Federal University of Juiz de Fora, Brazil

REFERENCES
"Wavelet-Based Estimators for Mixture Regression", 
submitted to Scandinavian Journal of Statistics.

WARNING
The function ksrlin.m presented below was developed by Yi Cao, and its 
license with copyright and terms of use is available in the file 
license_ksrlin.txt

LIST OF FILES:

SPECIFIC FUNCTIONS
  DWT.m              - Apply Mallat's algorithm to the data
  IDWT.m             - Apply inverse Mallat's algorithm to the data
  ksrlin.m           - Local Linear Smoother
  PHI.m              - Evaluates the scaling function at a vector of points
                       by Daubechies-Lagarias Algorithm
  WavCoefEsts.m      - Estimates wavelet coefficients of a function of probabilities 
                       in the context of the mixture regression problem
  WavShrink.m        - Regularizes wavelet coefficients by Hard and Soft 
                       threshold, and Efromovich's shrinkage method
  wd.m               - Applies wavelet decomposition based on DWT.m function
  wr.m               - Applies wavelet reconstruction based on IDWT.m function
	
SCRIPTS WITH APPLICATION
  Application.m      - Illustrates the methodology using a real data set

DATA SET
  Lai2005fig4.dat    - Real data set used in Section 5 of the paper. This data 
                   is available at
http://compbio.med.harvard.edu/Supplements/Bioinformatics05b/Profiles/Chrom_7_from40_to65Mb_GBM29.xls

LICENCE
  license_ksrlin.txt - This is the license of use of the function ksrlin.m with 
                       copyright and terms of use.


