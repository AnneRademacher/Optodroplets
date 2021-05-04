#' This function retrieves the mean intensity within a given mask
#'
#' @param im Image series
#' @return Intensity series (mean, standard deviation, area)

getIntensitySeries <- function(im, mask) {
  # retrieve the mean intensity, its standard deviation, and the area of the mask
  mean_intensity <- rep(0, numberOfFrames(im))
  sd_intensity <- rep(0, numberOfFrames(im))
  area <- rep(0, numberOfFrames(im))
  
  for(i in 1:numberOfFrames(im)) {
    roi <- mask[,,i]
    frame <- im[,,i]
    
    mean_intensity[i] <- mean(frame[roi>0])
    sd_intensity[i] <- sd(frame[roi>0])
    area[i] <- sum(roi>0)
  }
  
  # return results
  return(list(mean_intensity, sd_intensity, area))
}