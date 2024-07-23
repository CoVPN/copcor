draw.x.axis.cor=function(xlim, llox, llox.label, for.ggplot=FALSE){
  
  xx=seq(ceiling(xlim[1]), floor(xlim[2]))
  
  has.llox=F
  if (is.na(llox)) {
    labels = sapply (xx, function(x) if (x>=3) bquote(10^.(x)) else 10^x )
    
  } else if (llox.label=="delta") {
    labels = sapply (xx, function(x) if (x>=3 | x<=-3) bquote(10^.(x)) else 10^x )
    
  } else {
    
    xx=xx[xx>log10(llox*1.8)]
    labels = sapply (xx, function(x) if(x>=3) bquote(10^.(x)) else 10^x)
    xx=c(log10(llox), xx)
    labels=c(llox.label, labels)
    has.llox=T
  }
  
  # add e.g. 30 between 10 and 100
  if (length(xx)<4) {
    for (i in ifelse(has.llox,2,1):length(xx)) {
      x=xx[i]
      xx=c(xx, x+log10(3))
      labels=c(labels, if (x>=3) bquote(3%*%10^.(x)) else 3*10^(x) )
    }
    
    # add another tickmark at the very left
    if (!has.llox) {
      labels=c(if (min(xx)-1>=3) bquote(3%*%10^.(min(xx)-1)) else 3*10^(min(xx)-1), labels)    
      xx=c(min(xx)-1+log10(3), xx)
    }
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
report.assay.values=function(x, assay, grid_size=10){
  lars.quantiles=seq(0,1,length.out=30) [round(seq.int(1, 30, length.out = grid_size))]
  sens.quantiles=c(0.15, 0.85)
  # cannot have different lengths for different assays, otherwise downstream code may break
  fixed.values = log10(c("500"=500, "1000"=1000))
  # if we want to add "llox/2"=unname(lloxs[assay]/2))) to fixed.values, we have to get assay right, which will take some thought because marker.name.to.assay is hardcoded
  out=sort(c(quantile(x, c(lars.quantiles,sens.quantiles), na.rm=TRUE), fixed.values[fixed.values<max(x, na.rm=T) & fixed.values>min(x, na.rm=T)]))    
  out
  #out[!duplicated(out)] # unique strips away the names. But don't take out duplicates because 15% may be needed and because we may want the same number of values for each assay
}
#report.assay.values (dat.vac.seroneg[["Day57pseudoneutid80"]], "pseudoneutid80")



# 26Oct2020      Erika Rudnicki
theforestplot <- function(cohort=NA,group,nEvents=NA,totFU=NA,rate=NA,point.estimates,lower.bounds,upper.bounds,p.values,
                          table.labels,zero.line=1.0,dashed.line=NA,
                          x.ticks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2),
                          decimal.places = 1,fontsize = 1,width.pdf=7,height.pdf=7,graphwidth="auto",...){
  
  plotdata <- structure(
    list(
      mean  = c(NA, point.estimates),
      lower = c(NA, lower.bounds),
      upper = c(NA, upper.bounds)
    ),
    .Names = c("mean", "lower", "upper"),
    row.names = c(NA,-(length(point.estimates) + 1)),
    class = "data.frame"
  )
  
  # remove.redundancy <- function(x){ifelse(!is.na(dplyr::lag(x)) & x==dplyr::lag(x), NA, x)}
  show.decimals <- function(x){format(round(x, decimal.places), nsmall=decimal.places)}
  
  group_edit <- sapply(group, function(x){
    if(grepl(">=", x)){
      ssplit <- strsplit(x, split = ">=")[[1]]
      eval(parse(text = paste0('expression("', ssplit[1], '" >= "', ssplit[2], '")')))
    }else if(grepl("<", x)){
      ssplit <- strsplit(x, split = "<")[[1]]
      eval(parse(text = paste0('expression("', ssplit[1], '" < "', ssplit[2], '")')))
    }else{
      x
    }
  }, simplify = FALSE)
  group <- group_edit
  
  if(all(is.na(p.values))){
    tabletext <- list(
      # c(table.labels[1], remove.redundancy(as.character(cohort))),
      c(table.labels[1], group),
      c(table.labels[3], nEvents),
      # c(table.labels[4], totFU),
      # c(table.labels[5], rate),
      c(paste0(table.labels[2]), 
        paste0(sapply(point.estimates, show.decimals), " (", sapply(lower.bounds, show.decimals), ", ", sapply(upper.bounds, show.decimals), ")")),
      c(" ", rep(NA, length(point.estimates)))
    )
  } else{
    tabletext <- list(
      # c(table.labels[1], remove.redundancy(as.character(cohort))),
      c(table.labels[1], group),
      c(table.labels[3], nEvents),
      # c(table.labels[4], totFU),
      # c(table.labels[5], rate),
      c(paste0(table.labels[2]), 
        paste0(sapply(point.estimates, show.decimals), " (", sapply(lower.bounds, show.decimals), ", ", sapply(upper.bounds, show.decimals), ")")),
      c("P-value", p.values)
    )}    
  
  replaceNA <- function(x){ifelse(grepl("NA", x), NA, x)}
  tabletext[[3]] <- sapply(tabletext[[3]], replaceNA)
  
  replaceDash <- function(x){gsub("-", "\u2013", x)}
  tabletext[[3]] <- sapply(tabletext[[3]], replaceDash)
  
  if(!is.na(dashed.line)){grid.line <- structure(dashed.line, gp = gpar(lty = 2, col = "red", lwd=0.5))} else{grid.line <- FALSE}
  
  forestplot(tabletext,plotdata,is.summary = FALSE,col = fpColors(box = "darkblue",line = "darkblue",summary = "royalblue",zero="black"),    
             graph.pos = 3,graphwidth = graphwidth,
             hrzl_lines = list("2" = gpar(lty=1)),
             zero = zero.line,lwd.zero = 0.5,lwd.ci = 0.5,lwd.xaxis = 0.5,xticks = x.ticks,boxsize = 0.1,grid=grid.line,txt_gp = fpTxtGp(
               ticks = gpar(fontfamily = "", cex = fontsize * 0.8),
               label = gpar(fontfamily = "", cex = fontsize * 0.9),
               summary = gpar(cex = fontsize)
             ),
             colgap = unit(2, "mm"),align = c("l", "l", "l"),mar = unit(c(4,1,9,1), "mm"), #bltr
             clip = c(min(x.ticks), max(x.ticks)), ...
  )
}
