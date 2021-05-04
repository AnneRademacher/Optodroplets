#' This function segments the area occupied by aggregates
#'
#' @param im Image series
#' @return Segmentation result

makeAggregateMask <- function(im, threshold, nuc_mask, minsize,dilate_disc=3, erode_disc=3, ignore_boundary=11) {
  ## segmentation in the first image
  # find regions of interest
  mask <- im[,,1]>threshold
  mask <- dilate(mask, makeBrush(dilate_disc, shape='disc'))
  mask <- erode(mask, makeBrush(erode_disc, shape='disc'))
  mask[erode(nuc_mask[,,1], makeBrush(ignore_boundary, shape='disc'))==0] <- 0
  mask <- bwlabel(mask) # find connected sets of pixels
  mask <- fillHull(mask) # fill holes in the mask
  
  # only keep large regions of interest
  regions <- table(mask)
  small_regions <- which(regions<minsize)-1
  mask[mask %in% small_regions] <- 0
  mask[mask!=0] <- 1
  
  # registration for the following images
  m <- mask
  
  for(i in 2:numberOfFrames(im)) {
    mask <- im[,,i]>threshold
    mask <- dilate(mask, makeBrush(dilate_disc, shape='disc'))
    mask <- erode(mask, makeBrush(erode_disc, shape='disc'))
    mask[erode(nuc_mask[,,i], makeBrush(ignore_boundary, shape='disc'))==0] <- 0
    mask <- bwlabel(mask) # find connected sets of pixels
    mask <- fillHull(mask) # fill holes in the mask
    
    # only keep large regions of interest
    regions <- table(mask)
    small_regions <- which(regions<minsize)-1
    mask[mask %in% small_regions] <- 0
    mask[mask!=0] <- 1
    m <- abind(m, mask, along=3)
  }
  
  # return mask
  return(m)
}