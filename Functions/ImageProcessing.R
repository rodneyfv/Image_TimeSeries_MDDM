require(sp)
require(gstat) 
require(fields)
require(spam64) #?large_matrix


##############################################################
##############################################################
##############################################################

# here we apply the kriking method for all the images and save
# them in a folder that will be used later with Matlab

for(i in 1:87){
  # reading the image
  wd <- getwd()
  w <- scan(file = sprintf(
    "/home/rodney/Documents/Images_TimeSeries_files/ImagesTStxt/image-%i.txt",
    i), sep = ",")
  tmp <- matrix(w, ncol = 512)
  I <- t(tmp)
  rm(w); rm(tmp)
  # assuming the pixels are in the 
  # grid (i,j) \in \mathbb{Z}^d
  coord <- cbind(x = rep(1:512, 512), y = rep(1:512, each = 512))
  #fullData <- as.data.frame(cbind(coord, w))
  
  idx <- seq(1, 512, by = 1) 
  idy <- seq(1, 512, by = 1) 
  Is <- I[idx, idy]
  ws <- as.numeric(Is)
  vecIndex <- matrix(1:(nrow(I)*ncol(I)), nrow = nrow(I), ncol = ncol(I))[idx, idy]
  coords <- coord[vecIndex,]
  
  I.hat <- mKrig(y = ws, x = coords,
                 cov.function = "stationary.taper.cov",
                 lambda = 0.005/0.01, # nugget/variance
                 Covariance = "Exponential",
                 theta = 10, # Spatial dependence
                 Taper = "Wendland",
                 Taper.args = list(theta = 4, # approx neighborhood
                                   k = 2, 
                                   dimension = 2),
                 chol.args = list(pivot = TRUE,
                                  memory = list(nnzR= 900000)))
  w.hat <- as.numeric(predict(I.hat))
  # Distance matrix is relatively dense at theta = 4
  Is.hat <- matrix(w.hat, ncol = ncol(Is))
  
  rm(I.hat); rm(w.hat); rm(ws); rm(vecIndex)
  rm(I); rm(coord); rm(coords)
  # saving the treated image
  write.table(Is.hat,sprintf("%s/ImagesHatTS/imagehs-%i.txt",wd,i),
            row.names = FALSE, col.names = FALSE)
}
  


