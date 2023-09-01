draw.x.axis.cor=function(xlim, llox, llox.label, for.ggplot=FALSE){
  
  xx=seq(ceiling(xlim[1]), floor(xlim[2]))        
  if (is.na(llox)) {
    labels = sapply (xx, function(x) if (x>=3) bquote(10^.(x)) else 10^x )
    
  } else if (llox.label=="delta") {
    labels = sapply (xx, function(x) if (x>=3 | x<=-3) bquote(10^.(x)) else 10^x )
    
  } else {
    
    xx=xx[xx>log10(llox*1.8)]
    labels = sapply (xx, function(x) if(x>=3) bquote(10^.(x)) else 10^x)
    xx=c(log10(llox), xx)
    labels=c(llox.label, labels)
  }
  
  # add e.g. 30 between 10 and 100
  if (length(xx)<4) {
    for (i in 1:length(xx)) {
      x=xx[i]
      xx=c(xx, x+log10(3))
      labels=c(labels, if (x>=3) bquote(3%*%10^.(x)) else 3*10^(x) )
    }

    labels=c(if (min(xx)-1>=3) bquote(3%*%10^.(min(xx)-1)) else 3*10^(min(xx)-1), labels)    
    xx=c(min(xx)-1+log10(3), xx)
  }

  if (for.ggplot) {
    return(list(ticks = xx, labels = labels))
  } else {
    axis(1, at=xx, labels=sapply(labels, function (label) as.expression(label)))
  }
  
}


# get plotting range
get.xlim=function(dat, marker) {
  lloxs=with(assay_metadata, ifelse(llox_label=="lloq", lloq, lod))
  
  assay=marker.name.to.assay(marker)
  
  # the default
  ret=range(dat[[marker]], log10(lloxs[assay]/2), na.rm=T)
  
  # may be customized, e.g. to have the same xlim for different variants in the same type of assay
  # if (TRIAL=="moderna_boost") {
  #   if(assay %in% c("bindSpike", "bindRBD")) {
  #     ret=range(dat[[DayPrefix%.%time%.%"bindSpike"]], 
  #               dat[[DayPrefix%.%time%.%"bindRBD"]], 
  #               log10(lloxs[c("bindSpike","bindRBD")]/2), na.rm=T)
  #     
  #   } 
  # }

  delta=(ret[2]-ret[1])/20     
  c(ret[1]-delta, ret[2]+delta)
}


# get histogram object to add to VE plots etc
get.marker.histogram=function(marker, wt, trial, marker.break=marker) {
  # first call hist to get breaks, then call weighted.hist
  tmp.1=hist(marker.break,breaks=ifelse(trial=="moderna_real",25,15),plot=F)  # 15 is treated as a suggestion and the actual number of breaks is determined by pretty()
  tmp=weighted.hist(marker,wt, breaks=tmp.1$breaks, plot=F)
  attr(tmp,"class")="histogram" 
  tmp
}



get.range.cor=function(dat, assay, time) {
  lloxs=with(assay_metadata, ifelse(llox_label=="lloq", lloq, lod))
  
  if(assay %in% c("bindSpike", "bindRBD")) { # & all(c("pseudoneutid50", "pseudoneutid80") %in% assays)
    ret=range(dat[[DayPrefix%.%time%.%"bindSpike"]], 
              dat[[DayPrefix%.%time%.%"bindRBD"]], 
              log10(lloxs[c("bindSpike","bindRBD")]/2), na.rm=T)
    
  } else if(assay %in% c("pseudoneutid50", "pseudoneutid80")) { #  & all(c("pseudoneutid50", "pseudoneutid80") %in% assays)
    ret=range(dat[[DayPrefix%.%time%.%"pseudoneutid50"]], 
              dat[[DayPrefix%.%time%.%"pseudoneutid80"]], 
              #log10(uloqs[c("pseudoneutid50","pseudoneutid80")]),
              log10(lloxs[c("pseudoneutid50","pseudoneutid80")]/2), na.rm=T) 
  } else {
    ret=range(dat[[DayPrefix%.%time%.%assay]], 
              log10(lloxs[assay]/2), na.rm=T)        
  }
  delta=(ret[2]-ret[1])/20     
  c(ret[1]-delta, ret[2]+delta)
}




# x is the marker values
# assay is one of assays, e.g. pseudoneutid80
report.assay.values=function(x, assay){
  lars.quantiles=seq(0,1,length.out=30) [round(seq.int(1, 30, length.out = 10))]
  sens.quantiles=c(0.15, 0.85)
  # cannot have different lengths for different assays, otherwise downstream code may break
  fixed.values = log10(c("500"=500, "1000"=1000))
  # if we want to add "llox/2"=unname(lloxs[assay]/2))) to fixed.values, we have to get assay right, which will take some thought because marker.name.to.assay is hardcoded
  out=sort(c(quantile(x, c(lars.quantiles,sens.quantiles), na.rm=TRUE), fixed.values[fixed.values<max(x, na.rm=T) & fixed.values>min(x, na.rm=T)]))    
  out
  #out[!duplicated(out)] # unique strips away the names. But don't take out duplicates because 15% may be needed and because we may want the same number of values for each assay
}
#report.assay.values (dat.vac.seroneg[["Day57pseudoneutid80"]], "pseudoneutid80")

