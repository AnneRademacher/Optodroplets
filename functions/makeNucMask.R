#' This function segments the area occupied by the cell nucleus
#'
#' @param im Image series
#' @return Segmentation result

makeNucMask <- function(im, threshold, dilate_disc=3, erode_disc=3) {
  ## segmentation in the first image
  # find regions of interest
  mask <- im[,,1]>threshold[1]
  mask <- dilate(mask, makeBrush(dilate_disc, shape='disc'))
  mask <- erode(mask, makeBrush(erode_disc, shape='disc'))
  mask <- fillHull(mask) # fill holes in the mask
  mask <- bwlabel(mask) # find connected sets of pixels
  
  # select the largest region of interest
  regions <- table(mask)
  largest_region <- which(regions[2:length(regions)]==max(regions[2:length(regions)]))
  mask[mask != largest_region] <- 0
  
  ## segmentation in the following images
  for(i in 2:numberOfFrames(im)) {
    # find regions of interest
    buffer <- im[,,i]>threshold[i]
    buffer <- dilate(buffer, makeBrush(dilate_disc, shape='disc'))
    buffer <- erode(buffer, makeBrush(erode_disc, shape='disc'))
    buffer <- fillHull(buffer) # fill holes in the mask
    buffer <- bwlabel(buffer) # find connected sets of pixels
    
    # select the largest region of interest
    regions <- table(buffer)
    largest_region <- which(regions[2:length(regions)]==max(regions[2:length(regions)]))
    buffer[buffer != largest_region] <- 0
    
    # append current mask
    mask <- abind(mask, buffer, along=3)
  }
  
  # binarize
  mask[mask>1] <- 1 
  
  # return segmentation result
  return(mask)
}