utils::globalVariables(c("TRIAL", "config", "assay_metadata", "verbose"))


# extract assay name from marker names, which include Day, e.g.
# e.g. Day22pseudoneutid50 => pseudoneutid50, Delta22overBpseudoneutid50 => pseudoneutid50
marker.name.to.assay=function(a) {
  
  if (startsWith(a,"Day")) {
    # Day22pseudoneutid50 => pseudoneutid50 
    sub("Day[[0123456789]+", "", a)
    
  } else if (startsWith(a,"BD")) {
    # BD29pseudoneutid50 => pseudoneutid50
    sub("BD[[0123456789]+", "", a)
    
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
    dat.ph1$EventTimePrimary=followup.day
    risks = 1 - exp(-predict(fit.risk, newdata=dat.ph1, type="expected"))
    mean(risks)
  } else {
    # competing risk estimation
    out=pcr2(formula, dat.ph1, followup.day)
    mean(out)
  }
}



add.trichotomized.markers=function(dat, markers, ph2.col.name="ph2", wt.col.name="wt") {
  
  if(verbose) print("add.trichotomized.markers ...")
  
  marker.cutpoints <- list()    
  for (a in markers) {
    if (verbose) myprint(a, newline=F)
    tmp.a=dat[[a]]
    
    # if we estimate cutpoints using all non-NA markers, it may have an issue when a lot of subjects outside ph2 have non-NA markers
    # since that leads to uneven distribution of markers between low/med/high among ph2
    # this issue did not affect earlier trials much, but it is a problem with vat08m. We are changing the code for trials after vat08m
    if (TRIAL %in% c("hvtn705","hvtn705V1V2","hvtn705second","hvtn705secondRSA","hvtn705secondNonRSA","moderna_real","moderna_mock","prevent19",
                     "janssen_pooled_EUA","janssen_na_EUA","janssen_la_EUA","janssen_sa_EUA")) {
      flag=rep(TRUE, length(tmp.a))
    } else {
      flag=dat[[ph2.col.name]]
    }
    
    if(!startsWith(a, "Delta")) {
      # not fold change
      uloq=assay_metadata$uloq[assay_metadata$assay==marker.name.to.assay(a)]
      uppercut=log10(uloq); uppercut=uppercut*ifelse(uppercut>0,.9999,1.0001)
      lowercut=min(tmp.a, na.rm=T)*1.0001; lowercut=lowercut*ifelse(lowercut>0,1.0001,.9999)
      if (mean(tmp.a>uppercut, na.rm=T)>1/3) {
        # if more than 1/3 of vaccine recipients have value > ULOQ, let q.a be (median among those < ULOQ, ULOQ)
        if (verbose) cat("more than 1/3 of vaccine recipients have value > ULOQ\n")
        q.a=c(wtd.quantile(tmp.a[dat[[a]]<=uppercut & flag], weights = dat[[wt.col.name]][tmp.a<=uppercut & flag], probs = c(1/2)),  uppercut)
      } else if (mean(tmp.a<lowercut, na.rm=T)>1/3) {
        # if more than 1/3 of vaccine recipients have value at min, let q.a be (min, median among those > LLOQ)
        if (verbose) cat("more than 1/3 of vaccine recipients have at min\n")
        q.a=c(lowercut, wtd.quantile(tmp.a[dat[[a]]>=lowercut & flag], weights = dat[[wt.col.name]][tmp.a>=lowercut & flag], probs = c(1/2))  )
      } else {
        # this implementation uses all non-NA markers, which include a lot of subjects outside ph2, and that leads to uneven distribution of markers between low/med/high among ph2
        #q.a <- wtd.quantile(tmp.a, weights = dat[[wt.col.name]], probs = c(1/3, 2/3))
        q.a <- wtd.quantile(tmp.a[flag], weights = dat[[wt.col.name]][flag], probs = c(1/3, 2/3))
      }
    } else {
      # fold change
      q.a <- wtd.quantile(tmp.a[flag], weights = dat[[wt.col.name]][flag], probs = c(1/3, 2/3))
    }
    tmp=try(factor(cut(tmp.a, breaks = c(-Inf, q.a, Inf))), silent=T)
    
    do.cut=FALSE # if TRUE, use cut function which does not use weights
    # if there is a huge point mass, an error would occur, or it may not break into 3 groups
    if (inherits(tmp, "try-error")) do.cut=TRUE else if(length(table(tmp)) != 3) do.cut=TRUE
    
    if(!do.cut) {
      dat[[a %.% "cat"]] <- tmp
      marker.cutpoints[[a]] <- q.a
    } else {
      cat("\nfirst cut fails, call cut again with breaks=3 \n")
      # cut is more robust but it does not incorporate weights
      tmp=cut(tmp.a, breaks=3)
      stopifnot(length(table(tmp))==3)
      dat[[a %.% "cat"]] = tmp
      # extract cut points from factor level labels
      tmpname = names(table(tmp))[2]
      tmpname = substr(tmpname, 2, nchar(tmpname)-1)
      marker.cutpoints[[a]] <- as.numeric(strsplit(tmpname, ",")[[1]])
    }
    stopifnot(length(table(dat[[a %.% "cat"]])) == 3)
    if(verbose) {
      print(table(dat[[a %.% "cat"]]))
      cat("\n")
    }
  }
  
  attr(dat, "marker.cutpoints")=marker.cutpoints
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


