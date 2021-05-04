# load libraries
library(EBImage)
library(this.path)

#set working directory
setwd(this.dir())

# load functions
source("functions/makeNucMask.R")
source("functions/makeAggregateMask2.R")
source("functions/getIntensitySeries.R")

# list all folders and constructs to analyze (in this case only one)
folders <- list.dirs("example_data/", recursive = F)
constructs <- list.dirs("example_data/", full.names = F, recursive = F)
outfiles <- c("optodropletInduction.txt",
              "optodropletInduction.pdf")

# determine total number of images to analyze (sum of images in all folders)
n <- 0
for (j in 1:length(folders)){
  files <- list.files(folders[j], pattern=".tif", full.names=T)
  files <- files[!grepl("masks", files)]
  files_short <- list.files(path=folders[j], pattern=".tif", full.names=F)
  files_short <- files_short[!grepl("masks", files_short)]
  print(paste("Found", length(files), "images for", constructs[j]))
  n <- n+length(files)
}
print(paste("Analyzing a total of", n, "cells"))

# set segmentation parameters (to be adjusted for different constructs/microscopy settings)
a <- 0.95         #otsu scaling factor
b <- 5            #sigma for gblur for nucleus segmentation
b2 <- 2           #sigma for gblur for pre-existing aggregate segmentation
c <- 10           #intensity scaling factor for low-intensity images
c2 <- 2           #intensity scaling factor for all other images and for pre-existing aggregate segmentation
threshold <- 0.8  #threshold used for pre-existing aggregate segmentation
minsize <- 100    #minimum size of segmented pre-existing aggregates in pixels

# initialize result table
counter <- 0
results <- matrix(NA, nrow = n, ncol = 8)
results <- as.data.frame(results)
colnames(results) <- c("filename", "construct", "pre_mean", "pre_cv", "post_cv", "post2_cv", "A_lo", "A_hi")


##########################################
############# image analysis #############
##########################################


for (j in 1:length(folders)){
  print(paste("-------- Processing", constructs[j],"--------"))
  files <- list.files(folders[j], pattern=".tif", full.names=T)
  files <- files[!grepl("masks", files)]
  files_short <- list.files(folders[j], pattern=".tif", full.names=F)
  files_short <- files_short[!grepl("masks", files_short)]
  
  for(i in 1:length(files)) {
    
    # load images
    data <- channel(readImage(files[i]), "gray")
    pre <- data[,,1]
    image <- combine(pre, data[,,5], data[,,9])                          #combine the relevant frames #1, #5, #9, discard all others
    
    # enhance intensity, blur and remove aggregates for nucleus segmentation
    if (quantile(pre,0.8)<0.04){                                         #corresponds to a pixel intensity value of <11 in our images
      print(paste("Skipping very low intensity cell", files_short[i]))
      next
    }else if (quantile(pre,0.8)<0.09){                                   #corresponds to a pixel intensity value of <23 in our images
      test <- gblur(c*image, sigma = b)                                  #enhance the brightness and blur the image
      test[test>quantile(test,0.75)] <- quantile(test,0.75)              #remove very bright spots that hinder nucleus segmentation 
      test[test>1] <- 1                                                  #set the maximum intensity to 1 (i. e. a pixel intensity of 255)
    }else{
      test <- gblur(c2*image, sigma = b)
      test[test>quantile(test,0.75)] <- quantile(test,0.75)
      test[test>1] <- 1
    }
    
    # make nuclear masks
    nuc_mask <- makeNucMask(test, a*otsu(test))
    
    # use nucleus mask of image #1 in case segmentation failed for image #5 and/or #9
    if (abs(sum(nuc_mask[,,1]>0)-sum(nuc_mask[,,2]>0))>0.1*max(sum(nuc_mask[,,1]>0),sum(nuc_mask[,,2]>0))){
      nuc_mask[,,2] <- nuc_mask[,,1]
      print(paste("Image #5 segmentation failed for", files_short[i]))
    }
    if (abs(sum(nuc_mask[,,1]>0)-sum(nuc_mask[,,3]>0))>0.1*max(sum(nuc_mask[,,1]>0),sum(nuc_mask[,,3]>0))){
      nuc_mask[,,3] <- nuc_mask[,,1]
      print(paste("Image #9 segmentation failed for", files_short[i]))
    }
    
    # make image with masks (to visually check the segmentation results)
    image_with_masks <- toRGB(image)
    image_with_masks <- paintObjects(nuc_mask, image_with_masks, opac=c(1,0), col=c("yellow","yellow"))
    
    # segment out bright spots that are already present in the "pre" images if necessary
    if (quantile(pre,0.9999)>0.9){                                      #i. e. if the 26 brightest pixels (in a 512x512 px image) have an intensity value >=230, aggregates are segmented
      print(paste("Segmenting out large aggregates in", files_short[i]))
      test2 <- gblur(c2*image, sigma = b2)                              #enhance the brightness and blur the image
      test2[test2>1] <- 1                                               #set the maximum intensity to 1 (i. e. a pixel intensity of 255)
      aggregate_mask <- makeAggregateMask(test2, threshold, nuc_mask, minsize)
      final_mask <- nuc_mask-aggregate_mask
      image_with_masks <- paintObjects(aggregate_mask, image_with_masks, opac=c(1,0), col=c("magenta","magenta"))
    } else {
      final_mask <- nuc_mask
    }
    
    # optional: save images with masks for visual inspection of results
    #writeImage(image_with_masks, bits.per.sample=8, files=paste0(sub("\\..*", "", files[i]),"_with_masks.tif"))
    
    # measure CV
    image[final_mask==0] <- NA
    nucleus_series <- getIntensitySeries(image, final_mask)
    results[counter+i,] <- c(files_short[i], constructs[j], nucleus_series[[1]][1], nucleus_series[[2]]/nucleus_series[[1]], 0, 0)
    results$A_lo[counter+i] <- as.numeric(results$post_cv[counter+i])-as.numeric(results$pre_cv[counter+i])
    results$A_hi[counter+i] <- as.numeric(results$post2_cv[counter+i])-as.numeric(results$pre_cv[counter+i])
  }
  counter <- counter+length(files)
}


##########################################
############# result cleanup #############
##########################################


for (k in 3:8) {
  results[,k] <- as.numeric(results[,k])
}

#remove rows with NAs
results <- results[complete.cases(results),]

# save result table containing file name, constructs, expression levels, raw coefficients of variation and droplet abundance
write.table(x = results, quote = F, row.names = F, file=outfiles[1], sep = '\t')


##########################################
############## plot results ##############
##########################################


pdf(outfiles[2], width=6, height=5)

############
### A_lo ###
############

plot(c(0.03,0.6), c(0,1), las=1,type="n", main="Droplet formation", xlab="Expression level (a.u.)", ylab="Droplet abundance A_lo (a.u.)", log="x")
norm <- max(results$A_lo, na.rm = T)

for (j in 1:length(constructs)){
  cv <- results$A_lo[results$construct == constructs[j]]/norm
  mean <- results$pre_mean[results$construct == constructs[j]]
  
  points(x = mean, y = cv, pch=16, col=j)
  y <- cv; init <- list(k=0.1,n=1); up <- list(k=5,n=10); lw <- list(k=0,n=0)
  fitres <- nls(y ~ 1/(1+(k/mean)^n), start=init, lower=lw, upper=up, alg="port", control=list(maxiter = 50, tol = 1e-05, minFactor = 1e-07, printEval = FALSE, warnOnly = TRUE))
  lines(seq(0,1,by=0.001), 1/(1+(coefficients(fitres)[1]/seq(0,1,by=0.001))^coefficients(fitres)[2]), lwd=2, col=j)
  print(paste(c(paste(constructs[j],"c_crit ="), paste(constructs[j],"n =")), coefficients(fitres)))
}

legend("topleft",
       legend = constructs,
       col = seq.int(from = 1, to = length(constructs)),
       lwd = 2, pch = 16)


############
### A_hi ###
############

plot(c(0.03,0.6), c(0,1), las=1,type="n", main="Droplet formation", xlab="Expression level (a.u.)", ylab="Droplet abundance A_hi (a.u.)", log="x")
norm2 <- max(results$A_hi, na.rm = T)

for (j in 1:length(constructs)){
  cv <- results$A_hi[results$construct == constructs[j]]/norm2
  mean <- results$pre_mean[results$construct == constructs[j]]
  
  points(x = mean, y = cv, pch=16, col=j)
  y <- cv; init <- list(k=0.1,n=1); up <- list(k=5,n=10); lw <- list(k=0,n=0)
  fitres <- nls(y ~ 1/(1+(k/mean)^n), start=init, lower=lw, upper=up, alg="port", control=list(maxiter = 50, tol = 1e-05, minFactor = 1e-07, printEval = FALSE, warnOnly = TRUE))
  lines(seq(0,1,by=0.001), 1/(1+(coefficients(fitres)[1]/seq(0,1,by=0.001))^coefficients(fitres)[2]), lwd=2, col=j)
  print(paste(c(paste(constructs[j],"c_crit ="), paste(constructs[j],"n =")), coefficients(fitres)))
}

legend("topleft",
       legend = constructs,
       col = seq.int(from = 1, to = length(constructs)),
       lwd = 2, pch = 16)


dev.off()
