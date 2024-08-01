utils::globalVariables(c("TRIAL", "config", "assay_metadata", "verbose", "DayPrefix", "COR", "variant"))


# convert ID50 between assays

PPD2DKE_pseudoneutid50_D614G = function(log_titers) (log_titers + 0.371)/1.073
DKE2IU_pseudoneutid50_D614G  = function(log_titers)  log_titers - 0.616185 #log10(0.242)
MNG2IU_pseudoneutid50_D614G  = function(log_titers)  log_titers - 1.185087 #log10(0.0653)

PPD2DKE_pseudoneutid50_BA1   = function(log_titers) (log_titers + 0.303)/1.25
DKE2imputedIU_pseudoneutid50_BA1  = function(log_titers)  log_titers - 0.616185 #log10(0.242)
MNG2imputedIU_pseudoneutid50_BA1  = function(log_titers)  log_titers - 1.185087 #log10(0.0653)



# _ causes trouble in captions, and that has to be taken care of by putting \protect{} around the word
escape=function(x) {
  for (i in c("_","^")) {
    x=gsub(i, "\\"%.%i, x, fixed = TRUE)
  }
  x
}



# extract assay name from marker names, which include Day, e.g.
# e.g. Day22pseudoneutid50 => pseudoneutid50, Delta22overBpseudoneutid50 => pseudoneutid50
marker.name.to.assay=function(a) {
  
  if (startsWith(a,"Day")) {
    # Day22pseudoneutid50 => pseudoneutid50 
    sub("Day[[0123456789]+", "", a)
    
  } else if (startsWith(a,"BD")) {
    # BD29pseudoneutid50 => pseudoneutid50
    sub("BD[[0123456789]+", "", a)
    
  } else if (startsWith(a,"B")) {
    # Bpseudoneutid50 => pseudoneutid50
    substring(a, 2)
    
  } else if (startsWith(a,"M")) {
    # M18bindL1L2_HPV6 => bindL1L2_HPV6   
    sub("M[[0123456789]+", "", a)
    
  } else if (startsWith(a,"Mon")) {
    # Mon18bindL1L2_HPV6 => bindL1L2_HPV6   
    sub("Mon[[0123456789]+", "", a)
    
  } else if (contain(a,"overBD1")) {
    # DeltaBD29overBD1pseudoneutid50 => pseudoneutid50
    sub("DeltaBD[[0123456789]+overBD1", "", a)    
    
  } else if (contain(a,"overB")) {
    # Delta22overBpseudoneutid50 => pseudoneutid50
    sub("Delta[[0123456789]+overB", "", a)
    
  } else if (contain(a,"over")) {
    sub("Delta[[0123456789]+over[[0123456789]+", "", a)
    
  } else stop("marker.name.to.assay: not sure what to do")
}




# get marginalized risk to the followup followup.day without marker
get.marginalized.risk.no.marker=function(formula, dat.ph1, followup.day){
  if (!is.list(formula)) {
    # model=T is required because the type of prediction requires it, see Note on ?predict.coxph
    fit.risk = coxph(formula, dat.ph1, model=T) 
    # dat.ph1$EventTimePrimary=followup.day # it cannot be like this
    dat.ph1[[as.character(formula[[2]][[2]])]] = followup.day
    risks = 1 - exp(-predict(fit.risk, newdata=dat.ph1, type="expected"))
    mean(risks)
  } else {
    # competing risk estimation
    out=pcr2(formula, dat.ph1, followup.day)
    mean(out)
  }
}


# return a character column not a factor column, because different subpopulations may be cut at different steps
add.trichotomized.markers=function(dat, markers, ph2.col.name="ph2", wt.col.name="wt", verbose=F) {
  
  if(is.null(dat[[wt.col.name]])) stop("col does not exist: "%.%wt.col.name) 
  if(is.null(dat[[ph2.col.name]])) stop("col does not exist: "%.%ph2.col.name) 
  
  # this allows adding to marker.cutpoints
  if (is.null(attr(dat, "marker.cutpoints"))) {
    marker.cutpoints <- list()    
  } else {
    marker.cutpoints <- attr(dat, "marker.cutpoints")    
  }
  
  for (a in markers) {
    
    if (verbose) myprint(a, newline=T)
    
    if (is.null(dat[[a %.% "cat"]])) {
      new.col = T
      # if a new column, will return a column of factor

    } else {
      new.col = F
      # if not a new column, will return a column of type character
      dat[[a %.% "cat"]] = as.character(dat[[a %.% "cat"]])
    }
    
    tmp.a=dat[[a]]
    
    # if we estimate cutpoints using all non-NA markers, it may have an issue when a lot of subjects outside ph2 have non-NA markers
    # since that leads to uneven distribution of markers between low/med/high among ph2
    # this issue did not affect earlier trials much, but it is a problem with vat08m. We are changing the code for trials after vat08m
    if (TRIAL %in% c("hvtn705","hvtn705V1V2","hvtn705second","hvtn705secondRSA","hvtn705secondNonRSA",
                     "moderna_real","moderna_mock","prevent19",
                     "janssen_pooled_EUA","janssen_na_EUA","janssen_la_EUA","janssen_sa_EUA")) {
      flag=rep(TRUE, length(tmp.a))
    } else {
      flag=as.logical(dat[[ph2.col.name]])
    }
    
    # set it to NA so that it won't interfere with table
    tmp.a[!flag]=NA
    
    if(startsWith(a, "Delta")) {
      # fold change
      q.a <- wtd.quantile(tmp.a[flag], weights = dat[[wt.col.name]][flag], probs = c(1/3, 2/3))
      
    } else {
      # not fold change
      uloq=assay_metadata$uloq[assay_metadata$assay==marker.name.to.assay(a)]
      uppercut=log10(uloq); uppercut=uppercut*ifelse(uppercut>0,.9999,1.0001)
      
      lowercut=min(tmp.a, na.rm=T)*1.0001; 
      lowercut=lowercut*ifelse(lowercut>0,1.0001,.9999)
      pos.rate=mean(tmp.a[flag]>=lowercut, na.rm=T)
      myprint(pos.rate)
      
      binary.cut=FALSE
      if (mean(tmp.a[flag]>uppercut, na.rm=T)>1/3) {
        # if more than 1/3 of vaccine recipients have value > ULOQ, let q.a be (median among those < ULOQ, ULOQ)
        if (verbose) cat("more than 1/3 of vaccine recipients have value > ULOQ\n")
        q.a=c(wtd.quantile(tmp.a[dat[[a]]<=uppercut & flag], weights = dat[[wt.col.name]][tmp.a<=uppercut & flag], probs = c(1/2)),  uppercut)
      
      } else if (2/3<=pos.rate) {
        # this implementation uses all non-NA markers, which include a lot of subjects outside ph2, and that leads to uneven distribution of markers between low/med/high among ph2
        #q.a <- wtd.quantile(tmp.a, weights = dat[[wt.col.name]], probs = c(1/3, 2/3))
        q.a <- wtd.quantile(tmp.a[flag], weights = dat[[wt.col.name]][flag], probs = c(1/3, 2/3))
      
      } else if (1/3<=pos.rate & pos.rate<2/3) {
        # if more than 1/3 of vaccine recipients have value at min, let q.a be (min, median among those > LLOQ)
        if (verbose) cat("more than 1/3 of vaccine recipients have at min\n")
        q.a=c(lowercut, wtd.quantile(tmp.a[dat[[a]]>=lowercut & flag], weights = dat[[wt.col.name]][tmp.a>=lowercut & flag], probs = c(1/2))  )
        
      } else {
        # if pos.rate less than 1/3, make a binary variable
        binary.cut=TRUE
        cat("binary cut\n")
        q.a=c(lowercut)
        
      }
      
    }
    
    tmp=try(factor(cut(tmp.a, breaks = c(-Inf, q.a, Inf))), silent=T)
    if (inherits(tmp, "try-error")) {
      # if there is a huge point mass, an error would occur
      failed = TRUE
    } else if(!binary.cut & length(table(tmp)) != 3) {
      # or it may not break into 3 groups
      failed = TRUE
    } else {
      failed = FALSE
    }
    
    if(failed) {
      cat("\nfirst cut fails, call cut again with breaks=3 \n")
      # cut is more robust but it does not incorporate weights
      tmp=cut(tmp.a, breaks=3)
      stopifnot(length(table(tmp))==3)
    }
    
    if(new.col) {
      dat[[a %.% "cat"]] = as.character(tmp)
    } else {
      # only touch values in ph2
      dat[[a %.% "cat"]][flag] = as.character(tmp[flag])
    }
    
    # extract cut points from factor level labels
    # marker.cutpoints[[a]] <- q.a
    tmpname = names(table(tmp))[2]
    if(is.na(tmpname)) {
      # this happens when all values are the same
      tmpname = names(table(tmp))[1]
    }
    tmpname = substr(tmpname, 2, nchar(tmpname)-1)
    marker.cutpoints[[a]] <- as.numeric(strsplit(tmpname, ",")[[1]])
    
    # rm inf
    if (marker.cutpoints[[a]][2]==Inf) {
      marker.cutpoints[[a]] = marker.cutpoints[[a]][1]
    } else if (marker.cutpoints[[a]][1]==-Inf) {
      marker.cutpoints[[a]] = marker.cutpoints[[a]][2]
    }
      
    stopifnot(binary.cut & length(marker.cutpoints[[a]]) == 1 | length(table(tmp)) == 3 & length(marker.cutpoints[[a]]) == 2)

    if(verbose) {
      print(table(dat[dat[[ph2.col.name]]==1, a %.% "cat"]))
      cat("\n")
    }
    
  }
  
  attr(dat, "marker.cutpoints")=marker.cutpoints
  dat
  
}



add.dichotomized.markers=function(dat, markers, ph2.col.name="ph2", wt.col.name="wt", verbose=F) {
  
  if(is.null(dat[[wt.col.name]])) stop("col does not exist: "%.%wt.col.name) 
  if(is.null(dat[[ph2.col.name]])) stop("col does not exist: "%.%ph2.col.name) 
  
  # this allows adding to dich.marker.cutpoints
  if (is.null(attr(dat, "dich.marker.cutpoints"))) {
    dich.marker.cutpoints <- list()    
  } else {
    dich.marker.cutpoints <- attr(dat, "dich.marker.cutpoints")    
  }
  
  for (a in markers) {
    
    if (verbose) myprint(a, newline=T)
    
    if (is.null(dat[[a %.% "dich"]])) {
      new.col = T
      # if a new column, will return a column of factor
      
    } else {
      new.col = F
      # if not a new column, will return a column of type character
      dat[[a %.% "dich"]] = as.character(dat[[a %.% "dich"]])
    }
    
    tmp.a=dat[[a]]
    
    flag=as.logical(dat[[ph2.col.name]])

    # set it to NA so that it won't interfere with table
    tmp.a[!flag]=NA
    
    if(startsWith(a, "Delta")) {
      # fold change
      q.a <- wtd.quantile(tmp.a[flag], weights = dat[[wt.col.name]][flag], probs = c(1/2))
      
    } else {
      # not fold change
      uloq=assay_metadata$uloq[assay_metadata$assay==marker.name.to.assay(a)]
      uppercut=log10(uloq); uppercut=uppercut*ifelse(uppercut>0,.9999,1.0001)
      
      lowercut=min(tmp.a, na.rm=T)*1.0001; 
      lowercut=lowercut*ifelse(lowercut>0,1.0001,.9999)
      pos.rate=mean(tmp.a[flag]>=lowercut, na.rm=T)
      myprint(pos.rate)
      
      if (mean(tmp.a[flag]>uppercut, na.rm=T)>0.4) {
        if (verbose) cat("more than 40% of vaccine recipients have value > ULOQ, cut at ULOQ\n")
        q.a=uppercut
        
      } else if (1-pos.rate>0.4) {
        if (verbose) cat("pos.rate less than 0.6, cut at pos threshold\n")
        q.a=c(lowercut)
        
      } else {
        q.a <- wtd.quantile(tmp.a[flag], weights = dat[[wt.col.name]][flag], probs = c(1/2))
        
      }
      
    }
    
    tmp=try(factor(cut(tmp.a, breaks = c(-Inf, q.a, Inf))), silent=T)
    if (inherits(tmp, "try-error")) {
      # if there is a huge point mass, an error would occur
      failed = TRUE
    } else if(length(table(tmp)) != 2) {
      # or it may not break into 3 groups
      failed = TRUE
    } else {
      failed = FALSE
    }
    
    if(failed) {
      cat("\nfirst cut fails, call cut again with breaks=2 \n")
      # cut is more robust but it does not incorporate weights
      tmp=cut(tmp.a, breaks=2)
      stopifnot(length(table(tmp))==2)
    }
    
    if(new.col) {
      dat[[a %.% "dich"]] = tmp
      
      # extract cut points from factor level labels
      # dich.marker.cutpoints[[a]] <- q.a
      print(names(table(tmp))[2])
      tmpname = names(table(tmp))[2]
      tmpname = substr(tmpname, 2, nchar(tmpname)-1)
      dich.marker.cutpoints[[a]] <- as.numeric(strsplit(tmpname, ",")[[1]][1])
      print(dich.marker.cutpoints[[a]])
      
    } else {
      # only update values in ph2
      dat[[a %.% "dich"]][flag] = as.character(tmp[flag])
      
    }
    
    if(verbose) {
      print(table(dat[dat[[ph2.col.name]]==1, a %.% "dich"]))
      cat("\n")
    }
    
  }
  
  attr(dat, "dich.marker.cutpoints")=dich.marker.cutpoints
  dat
  
}




# a function to print tables of cases counts with different marker availability
# note that D57 cases and intercurrent cases may add up to more than D29 cases because ph1.D57 requires EarlyendpointD57==0 while ph1.D29 requires EarlyendpointD29==0
make.case.count.marker.availability.table=function(dat) {
  if (config$study_name=="COVE" | config$study_name=="MockCOVE" ) {
    idx.trt=1:0
    names(idx.trt)=c("vacc","plac")
    cnts = sapply (idx.trt, simplify="array", function(trt) {
      idx=1:3
      names(idx)=c("Day 29 Cases", "Day 57 Cases", "Intercurrent Cases")
      tab=t(sapply (idx, function(i) {           
        tmp.1 = with(dat[with(dat, Trt==trt & Bserostatus==0 & (if(i==2) EventIndPrimaryD57 else EventIndPrimaryD29) &   (if(i==2) ph1.D57 else if(i==1) ph1.D29 else ph1.intercurrent.cases)),], is.na(BbindSpike)     | is.na(BbindRBD) )
        tmp.2 = with(dat[with(dat, Trt==trt & Bserostatus==0 & (if(i==2) EventIndPrimaryD57 else EventIndPrimaryD29) &   (if(i==2) ph1.D57 else if(i==1) ph1.D29 else ph1.intercurrent.cases)),], is.na(Day29bindSpike) | is.na(Day29bindRBD))
        tmp.3 = with(dat[with(dat, Trt==trt & Bserostatus==0 & (if(i==2) EventIndPrimaryD57 else EventIndPrimaryD29) &   (if(i==2) ph1.D57 else if(i==1) ph1.D29 else ph1.intercurrent.cases)),], is.na(Day57bindSpike) | is.na(Day57bindRBD))    
        
        c(sum(tmp.1 & tmp.2 & tmp.3), sum(tmp.1 & tmp.2 & !tmp.3), sum(tmp.1 & !tmp.2 & tmp.3), sum(tmp.1 & !tmp.2 & !tmp.3), 
          sum(!tmp.1 & tmp.2 & tmp.3), sum(!tmp.1 & tmp.2 & !tmp.3), sum(!tmp.1 & !tmp.2 & tmp.3), sum(!tmp.1 & !tmp.2 & !tmp.3))
      }))
      colnames(tab)=c("---", "--+", "-+-", "-++", "+--", "+-+", "++-", "+++")
      tab
    })
    cnts
  } else if (config$study_name=="ENSEMBLE" | config$study_name=="MockENSEMBLE" ) {
    idx.trt=1:0
    names(idx.trt)=c("vacc","plac")
    cnts = sapply (idx.trt, simplify="array", function(trt) {
      idx=1:1
      tab=t(sapply (idx, function(i) {           
        tmp.1 = with(dat[with(dat, Trt==trt & Bserostatus==0 & if(i==2) EventIndPrimaryD57 else EventIndPrimaryD29 &   if(i==2) ph1.D57 else if(i==1) ph1.D29 else ph1.intercurrent.cases),], is.na(BbindSpike)     | is.na(BbindRBD) )
        tmp.2 = with(dat[with(dat, Trt==trt & Bserostatus==0 & if(i==2) EventIndPrimaryD57 else EventIndPrimaryD29 &   if(i==2) ph1.D57 else if(i==1) ph1.D29 else ph1.intercurrent.cases),], is.na(Day29bindSpike) | is.na(Day29bindRBD))
        
        c(sum(tmp.1 & tmp.2), sum(!tmp.1 & tmp.2), sum(tmp.1 & !tmp.2), sum(!tmp.1 & !tmp.2))
      }))
      colnames(tab)=c("--", "+-", "-+", "++")
      tab
    })
    t(drop(cnts))
  } else {
    NA
  }
}
#make.case.count.marker.availability.table(dat.mock)


