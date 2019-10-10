Description of the files in the folder 'Images_matlab_87'

##########################################################

Folder: Wavelet_Density_Estimation_PV1997
Description:
Folder containing the functions used to estimate a density function
using the method proposed by Pinheiro and Vidakovic (1997)

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


File: Forest_01.m
Description:
In this code we perform the following analyses:
1 - We compute the wavelet decomposition of the images, then
we estimate the density of the coefficients for each subband
using binning, get this estimated density evaluated on
some points and the take the 1D wavelet transform of these 
points. For each subband, the time series of these points are
used with the method in Fonseca and Pinheiro (2019) to
estimate the dimension of the subspace gererating these densities.
2 - We plot the time series densities of the subband's coefficients.
3 - Use the estimated dimension of the subspace generating the
densities to obtain estimates of these densities with the 
estimated dimension, and then compute the MDDMs
4 - test using the function wavedec2 instead of the swt2 to
perform the 2D wavelet decomposition of the images. Also, in
this part the functions used are adapted to consider different values of J.
5 - in the last part, the mean image is computed, subtracted
of all images, and then the MDDMs are computed for these normalized images.


File: Forest_02.m
Description:
In this code we perform the following analyses:
1 - we use the function 'DecompForestImageTS_swt2.m' to compute
the wavelet coefficients of the 1D-DWT of the densities
estimated for subbands obtained from the 2D-DWT of the images
using the swt2 function of Matlab. Then we estimate the
dimension of the subspace generating these densities.
2 - a video is made with the plot of the estimated density
function corresponding to each subband of the 2D-DWT of the
images and each time point.
3 - MDDMs are tested to be plotted in gray scale
4 - using the estimated dimension of the subspace generating
the densities, estimates of these densities are obtained
and used to compute the MDDMs.
5 - the mean image is computed, subtractedof all images, and 
then the MDDMs are computed for these normalized images.


File: Forest_03.m
Description:
In this code we perform the following analyses:
1 - we plot the energy (sum of squared coefficients) on each level of a
2D-DWT of the mean image. We look for the level where we can observe a
expressive growth begins.
2 - we compute the mean image, then we take logs, and test the MDDMs
in the cases when normalization (subtract the mean) is performed in the
images before smoothing and the case when smoothing is performed before
normalization. The number of MDDMs depends on J, and is related for each
detail subband generated
3 - we identify the image that differs the most from all the others and the
two most different time points
4 - we make a video with the image time series


File: Forest_04.m
Description:
1 - Here we test the idea of estimating the dimension of the subspace
generating the density time series of wavelet coefficients, then using
the densities estimated with the eigenfunctions and the dimension chosen
to obtain residuals, which are employed to obtain a bootstrap
distribution of the Hellinger distance under H0.
2 - Here we make the same analysis as before, but also using the
approximation coefficients.


File: Rondonia_DimTS_01.m
Description:
In this code we analyze the time series of density functions obtained for
each sub-region of the 2D wavelet decomposition of the images of
Rondonia. After getting these time series, we estimate their dimension
using the technique proposed in the paper Fonseca and Pinheiro (2019).

