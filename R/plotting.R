
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
get.xlim=function(dat, marker, lloxs) {
  assay=marker.name.to.assay(marker)
  
  # the default
  ret=range(dat[[marker]], log10(lloxs[assay]/2), na.rm=T)
  
  # may be customized, e.g. to have the same xlim for different variants in the same type of assay
  # if (TRIAL=="moderna_boost") {
  #   if(assay %in% c("bindSpike", "bindRBD")) {
  #     ret=range(dat[["Day"%.%time%.%"bindSpike"]], 
  #               dat[["Day"%.%time%.%"bindRBD"]], 
  #               log10(lloxs[c("bindSpike","bindRBD")]/2), na.rm=T)
  #     
  #   } 
  # }

  delta=(ret[2]-ret[1])/20     
  c(ret[1]-delta, ret[2]+delta)
}


