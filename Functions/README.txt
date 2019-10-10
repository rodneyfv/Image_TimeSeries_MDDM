Description of the files in the folder 'Functions'

##########################################################

File: smooth2d.m
Description:
Function that applies to a matrix H the smoothing method 
described in Eilers and Goeman (2004, Bioinformatics).
First, columns are smoothed, then rows are smoothed.
lambda is a tunning parameter.


File: MatGramSchm.m
Description:
function to apply the Gram-Schmidt to the vectors of the matrix Y


file: Hellinger_PV97_interp_denoise.m
Description:
This function computes the total Hellinger distance of two images im1 and
im2. The density function of the wavelet coefficents in each sub-matrix
of the decomposition of these images, is estimated non-parametricaly
using the wavelet density estimator proposed by Pinheiro and Vidakovic
(1997). Here, the wavelet functions employed are computed based on pre-
computed values of the scaling and wavelet functions.


File: Hellinger_WT_interp.m
Description:
This function computes the total Hellinger distance of two images im1 and
im2. The density function of the wavelet coefficents in each sub-matrix
of the decomposition of these images, is estimated non-parametricaly
using the wavelet density estimator proposed by Pinheiro and Vidakovic
(1997). Here, the wavelet functions employed are computed based on pre-
computed values of the scaling and wavelet functions.


File: function [dist,distE]=Hellinger_ksdensity(im1,im2,wname,J)
Description:
This function computes the total Hellinger distance of two images im1 and
im2. The density function of the wavelet coefficents in each sub-matrix
of the decomposition of these images, is estimated non-parametricaly
using a kernel estimator.


File: Hellinger_DWT.m
Description:
This function computes the total Hellinger distance of two images im1 and
im2. A binning method is applied to the data, and an approximation of the
density function is computed at equispaced points. The values obtained
are decomposed with the DWT, whose wavelet coefficients generated are
used to compute the Hellinger distance.


File: Estim_Dim_Pval( Xdec, p, Nboot, alpha )
Description:
in this function we test the dimension of curves time series using the
method of Fonseca and Pinheiro (2019). The function receives a matrix of
wavelet coefficients of the curves, a lag p and performs bootstrap tests
with Nboot replications and significance level alpha. The function
returns the selected dimension vDimSel, the p-values obtained in the
bootstrap tests and a vector with eigenvalues associanted to the
eigenfunctions.


File: estimated_functions_Dim.m
Description:
this function receives a matrix Xdec with wavelet coefficients and using
a lag of p, it returns the coefficients mY of a estimated process of
dimension d.


File: DimEst_wavestrap.m
Description:
Nessa função é realizado o procedimento bootstrap para estimar a 
dimensão do processo, mas diferentemente do que é feito em DimEst_boot.m,
aqui aproveitamos as decomposições obtidas para realizar a reamostragem
nos termos formados com os coeficientes, o que deixa o código mais rápido.
A é a matriz de coeficientes de ondaletas dos funcionais observados, 
nas colunas, para diferentes dias, nas linhas.
NREP é o número de réplicas bootstrap. p é lag máximo. 


File: DecompForestImageTS_wavedec2_ThreshNormalize
Description:
Description: Here we take logs of the t-th image and of the mean image,
then we perform a smoothing stage and after that we normalize the t-th
image. Finally, we take a 2D-DWT, estimate the density function of the
coefficient's distribution for each detail subband separetely, and store
the wavelet coefficients of these densities on matrices that are returned
Input
idx and idy: indices of the x and y axis that are used in the Forest images
img_mean: mean image of the the Forest images
t: time of the image to be analyzed
J: level used in the 2D wavelet decomposition with wavedec2
wname: name of the wavelet used
npts: number of points used with the binning estimate of the coefficients densities
j1: level of the 1D wavelet decomposition of the wavelet densities
Output
xdec_H: matrix with wavelet coefficiets of the decomposition of the
estimated density of horizontal details, for each of the J detail
subbands
xdec_V and xdec_D: similar to xdec_H, but for vertical and diagonal details
L1: number of coefficients on each scale of the wavelet decomposition of
the coefficient's densities


File: DecompForestImageTS_wavedec2_NormalizeThresh.m
Description:
Description: Here we take logs of the t-th image and of the mean image,
then we normalize then and later perform a smoothing on their difference.
Finally, we take a 2D-DWT, estimate the density function of the
coefficient's distribution for each detail subband separetely, and store
the wavelet coefficients of these densities on matrices that are returned
Input
idx and idy: indices of the x and y axis that are used in the Forest images
img_mean: mean image of the the Forest images
t: time of the image to be analyzed
J: level used in the 2D wavelet decomposition with wavedec2
wname: name of the wavelet used
npts: number of points used with the binning estimate of the coefficients densities
j1: level of the 1D wavelet decomposition of the wavelet densities
Output
xdec_H: matrix with wavelet coefficiets of the decomposition of the
estimated density of horizontal details, for each of the J detail subbands
xdec_V and xdec_D: similar to xdec_H, but for vertical and diagonal details
L1: number of coefficients on each scale of the wavelet decomposition of
the coefficient's densities


File: DecompForestImageTS_wavedec2.m
Description:
With this function we access the Forest data, get the image corresponding
to time t and at the pixels idx and idy, then perform the MDDM process
using the wavedec2 function for the 2D wavelet decomposition of the
image minus the mean image, which is stored in img_mean. J is the level 
of wavedec2 transform, J_MDDM is the maximum level
of detail coefficients we analyze when computing the density function of 
the wavelet coefficients, wname is the wavelet type used, npts is
the number of points used in the estimated density with the binning and
j1 is the level used in the thresholding of the estimated density.


File: DecompForestImageTS_swt2.m
Description:
With this function we access the Forest data, get the image corresponding
to time t and at the pixels idx and idy, then perform the MDDM process
using the swt2 function for the 2D wavelet decomposition of the image. J
is the level of swt2 transform, wname is the wavelet type used, npts is
the number of points used in the estimated density with the binning and
j1 is the level used in the thresholding of the estimated density.


File: binning.m
Description:
Performs a binning (histogram estimator) of the data on vx, and
evaluates the estimated density function on the points on vpts. The
length of the window used is fixed on wind.


