add.wt=function(dat, ph1="ph1", ph2="ph2", Wstratum="Wstratum", wt="wt", verbose=FALSE) {
  
  tmp = dat[[ph1]]==1
  wts_table <- table(dat[tmp,Wstratum], dat[tmp,ph2])
  
  if (any(0==wts_table[, 2] & 0!=wts_table[, 1])) {
    print(wts_table)
    stop("there exist empty strata")
  }
  if(verbose) print(wts_table)
  
  wts_norm <- rowSums(wts_table) / wts_table[, 2]
  dat[[wt]] <- wts_norm[dat[[Wstratum]] %.% ""]
  dat[dat[[ph1]]!=1, wt] = NA   # the step above assigns weights for some subjects outside ph1. the next step makes them NA
  dat

}