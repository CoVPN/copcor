###################################################################################################
# bootstrap marginalized risks
# type: 
#    1 for S=s
#    2 for S>=s
#    3 for categorical S
# data: ph1 data
# t: a time point near to the time of the last observed outcome will be defined
marginalized.risk.svycoxph.boot=function(form.0, marker.name, type, data, t, B, ci.type="quantile", numCores=1, additional.terms=NULL, variant=NULL) {
#marker.name=a; type=1; data=dat.vac.seroneg; t=tfinal.tpeak; B=B; ci.type="quantile"; numCores=1; additional.terms=NULL; variant=NULL
  
  # store the current rng state 
  save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
  if (inherits(save.seed, "try-error")) {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) } 
  
  # input check
  if (TRIAL %in% c("janssen_partA_VL") & is.null(variant)) stop("Variant needs to be specified for janssen_partA_VL")
    
  
  data.ph2=subset(data, data$ph2==1)
  
  # used in both point est and bootstrap
  # many variables are not passed but defined in the scope of marginalized.risk.svycoxph.boot
  fc.1=function(data.ph2, data, f1, categorical.s, n.dean=FALSE, in.boot=FALSE){
    # This is no longer necessary b/c pcr2 is updated to handle the situation when there are no competing events
    # if (comp.risk) 
    #   if(all(model.frame(f1[[2]], data.ph2)[[1]][,2]==0)) 
    #     # if there are no competing events, drop competing risk formula
    #     f1 = f1[[1]]
    # 
    if (comp.risk) {
      # competing risk implementation
      newdata=data.ph2
      out=sapply(ss, function(x) {
        newdata[[marker.name]]=x
        risks = try(pcr2(f1, data.ph2, t, weights=data.ph2$wt, newdata=newdata), silent=in.boot)
        ifelse (inherits(risks, "try-error"), NA, weighted.mean(risks, data.ph2$wt))
      })
      if (n.dean) c(NA,out) else out

    } else {        
      # non-competing risk implementation
      result <- tryCatch({
        fit.risk.1=svycoxph(f1, design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=data))
        out=marginalized.risk(fit.risk.1, marker.name, data.ph2, t=t, ss=ss, weights=data.ph2$wt, categorical.s=categorical.s)
        if (n.dean) {
          c(n.dean= last(coef(fit.risk.1)/sqrt(diag(fit.risk.1$var))) * sqrt(1/fit.risk.1$n + 1/fit.risk.1$nevent), out) 
        } else out
      }, 
      warning = function(w) {
        rep(NA, ifelse(n.dean,1,0)+length(ss))
      },
      error = function(e) {
        rep(NA, ifelse(n.dean,1,0)+length(ss))
      },
      finally = {
        # cat("This runs no matter what!\n")
      })
      result
      
    } 
  }
  
  fc.2=function(data.ph2, form.0, in.boot=FALSE){
    # This is no longer necessary b/c pcr2 is updated to handle the situation when there are no competing events
    # if (comp.risk) 
    ##    The following is not enough to ensure there are no competing events b/c subsets of datasets are used in pcr2
    #   if(all(model.frame(form.0[[2]], data.ph2)[[1]][,2]==0)) 
    #     # if there are no competing events, drop competing risk formula
    #     form.0 = form.0[[1]]
    # 
    if (comp.risk) {
      sapply(ss, function(x) {
        if (verbose>=2) myprint(x)
        newdata=data.ph2[data.ph2[[marker.name]]>=x, ]
        risks=try(pcr2(form.0, newdata, t, weights=newdata$wt), silent=in.boot)
        ifelse (inherits(risks, "try-error"), NA, weighted.mean(risks, newdata$wt))
      })
    } else {
      out = try(marginalized.risk.threshold (form.0, marker.name, data=data.ph2, weights=data.ph2$wt, t=t, ss=ss))
      if ( !inherits(out, "try-error" )) {
        out
      } else {
        rep(NA, length(ss)) 
      }
    }
  }    
  
  f2=as.formula(paste0("~.+",marker.name, if(!is.null(additional.terms)) "+"%.%additional.terms))
  
  comp.risk=is.list(form.0)
  if (comp.risk) {
    f1=lapply(form.0, function(x) update(x, f2))
  } else {
    f1=update(form.0, f2)        
  }
  
  
  if (type==1) {
    # conditional on S=s (quantitative)
    # don't sort ss or do ss=ss[!duplicated(ss)] because e.g. 15% will be lost and later code depends on that
    ss=sort(c(
      # Lars quantiles so that to be consistent with his analyses, also add every 5% to include s1 and s2 for sensitivity analyses
      report.assay.values(data[[marker.name]][data$EventIndPrimary==1], marker.name.to.assay(marker.name)), 
      # 2.5% and 97.5% as the leftmost and rightmost points 
      wtd.quantile(data[[marker.name]], data$wt, c(0.025,0.05,0.95,0.975)),
      # equally spaced values so that the curves look good  
      seq(min(data[[marker.name]], na.rm=TRUE), max(data[[marker.name]], na.rm=TRUE), length=100)[-c(1,100)],
      # useful for reports
      if (log10(100)>min(data[[marker.name]], na.rm=TRUE) & log10(100)<max(data[[marker.name]], na.rm=TRUE)) log10(100)
    ))
    
    prob = if (TRIAL %in% c("janssen_partA_VL")) {
      # do MI
      rowMeans(sapply(1:10, function(imp) {
        data.ph2$EventIndOfInterest = ifelse(data.ph2$EventIndPrimary==1 & data.ph2[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
        data.ph2$EventIndCompeting  = ifelse(data.ph2$EventIndPrimary==1 & data.ph2[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
        data$EventIndOfInterest = ifelse(data$EventIndPrimary==1 & data[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
        data$EventIndCompeting  = ifelse(data$EventIndPrimary==1 & data[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
        fc.1(data.ph2, data, f1, n.dean=TRUE, categorical.s=FALSE)
      }))
    } else {
      fc.1(data.ph2, data, f1, n.dean=TRUE, categorical.s=FALSE)
    }
    
    n.dean=prob[1]
    prob=prob[-1]

  } else if (type==2) {
    # conditional on S>=s
    ss=quantile(data[[marker.name]], seq(0,.9,by=0.05), na.rm=TRUE)
    if(verbose>=2) myprint(ss)
    prob = if (TRIAL %in% c("janssen_partA_VL")) {
      # do MI
      rowMeans(sapply(1:10, function(imp) {
        data.ph2$EventIndOfInterest = ifelse(data.ph2$EventIndPrimary==1 & data.ph2[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
        data.ph2$EventIndCompeting  = ifelse(data.ph2$EventIndPrimary==1 & data.ph2[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
        fc.2(data.ph2, form.0)
      }))
    } else {
        fc.2(data.ph2, form.0)
    }
    
    
  } else if (type==3) {
    # conditional on S=s (categorical)
    ss=unique(data[[marker.name]]); ss=sort(ss[!is.na(ss)])
    if(verbose>=2) myprint(ss)        
    prob = if (TRIAL %in% c("janssen_partA_VL")) {
      # do MI
      rowMeans(sapply(1:10, function(imp) {
        data.ph2$EventIndOfInterest = ifelse(data.ph2$EventIndPrimary==1 & data.ph2[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
        data.ph2$EventIndCompeting  = ifelse(data.ph2$EventIndPrimary==1 & data.ph2[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
        data$EventIndOfInterest = ifelse(data$EventIndPrimary==1 & data[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
        data$EventIndCompeting  = ifelse(data$EventIndPrimary==1 & data[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
        fc.1(data.ph2, data, f1, n.dean=FALSE, categorical.s=TRUE)
      }))
    } else {
      fc.1(data.ph2, data, f1, n.dean=FALSE, categorical.s=TRUE)
    }
    
    
  } else if (type==4) {
    # conditional on S=s (quantitative)
    if (comp.risk) {
      stop("need to implement this (like type 1 but coef only)") 
    } else {
      tmp.design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=data)
      fit.risk=try(svycoxph(f1, design=tmp.design)) # since we don't need se, we could use coxph, but the weights computed by svycoxph are a little different from the coxph due to fpc
    }
    
  } else stop("wrong type")
  
  
  ###############
  # bootstrap
  
  if(config$sampling_scheme=="case_cohort") ptids.by.stratum=get.ptids.by.stratum.for.bootstrap (data)     
  
  seeds=1:B; names(seeds)=seeds

  out=mclapply(seeds, mc.cores = numCores, FUN=function(seed) {   
    seed=seed+560
    if (verbose>=2) myprint(seed)
    
    if (TRIAL=="moderna_boost") {
      dat.b = bootstrap.cove.boost.2(data, seed)
      
    } else if(config$sampling_scheme=="case_cohort") {
      dat.b = get.bootstrap.data.cor (data, ptids.by.stratum, seed) 
      
    } else if(config$sampling_scheme=="case_control") {
      # includes hvtn705second
      # min.cell.size default to 1
      dat.b = bootstrap.case.control.samples(data, seed, delta.name="EventIndPrimary", strata.name="tps.stratum", ph2.name="ph2") 
      
    } else if(config$sampling_scheme=="cohort") {
      dat.b = bootstrap.cohort(data, seed) 
      
    } else stop("not sure which bootstrap function to use")
    
    dat.b.ph2=subset(dat.b, dat.b$ph2==1)  
    
    # if there is no missing variant info in a bootstrap dataset, only need to run the MI code once
    if (TRIAL %in% c("janssen_partA_VL")) {
      nImp = ifelse(any(with(subset(dat.b.ph2, dat.b.ph2$EventIndPrimary==1), is.na(seq1.variant))), 10, 1)
    }
    
    if(type==1) {
      # conditional on s
      if (TRIAL %in% c("janssen_partA_VL")) {
        rowMeans(sapply(1:nImp, function(imp) {
          dat.b.ph2$EventIndOfInterest = ifelse(dat.b.ph2$EventIndPrimary==1 & dat.b.ph2[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
          dat.b.ph2$EventIndCompeting  = ifelse(dat.b.ph2$EventIndPrimary==1 & dat.b.ph2[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
          dat.b$EventIndOfInterest = ifelse(dat.b$EventIndPrimary==1 & dat.b[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
          dat.b$EventIndCompeting  = ifelse(dat.b$EventIndPrimary==1 & dat.b[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
          fc.1(dat.b.ph2, dat.b, f1, n.dean=TRUE, categorical.s=FALSE, in.boot=T)
        }))
      } else {
        fc.1(dat.b.ph2, dat.b, f1, n.dean=TRUE, categorical.s=FALSE, in.boot=T)
      }
      
    } else if (type==2) {
      # conditional on S>=s
      if (TRIAL %in% c("janssen_partA_VL")) {
        rowMeans(sapply(1:nImp, function(imp) {
          dat.b.ph2$EventIndOfInterest = ifelse(dat.b.ph2$EventIndPrimary==1 & dat.b.ph2[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
          dat.b.ph2$EventIndCompeting  = ifelse(dat.b.ph2$EventIndPrimary==1 & dat.b.ph2[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
          fc.2(dat.b.ph2, form.0, in.boot=T)
        }))
      } else {
        fc.2(dat.b.ph2, form.0, in.boot=T)
      }
      
    } else if (type==3) {
      # conditional on a categorical S
      if (TRIAL %in% c("janssen_partA_VL")) {
        rowMeans(sapply(1:nImp, function(imp) {
          dat.b.ph2$EventIndOfInterest = ifelse(dat.b.ph2$EventIndPrimary==1 & dat.b.ph2[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
          dat.b.ph2$EventIndCompeting  = ifelse(dat.b.ph2$EventIndPrimary==1 & dat.b.ph2[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
          dat.b$EventIndOfInterest = ifelse(dat.b$EventIndPrimary==1 & dat.b[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
          dat.b$EventIndCompeting  = ifelse(dat.b$EventIndPrimary==1 & dat.b[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
          fc.1(dat.b.ph2, dat.b, f1, n.dean=FALSE, categorical.s=TRUE, in.boot=T)
        }))
      } else {
        fc.1(dat.b.ph2, dat.b, f1, n.dean=FALSE, categorical.s=TRUE, in.boot=T)
      }
      
    } else if (type==4) {
      # conditional on S=s (quantitative)
      fit.risk.b=try(svycoxph(f1, design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.b)))
      if ( !inherits(fit.risk.b, "try-error")) {
      } else {
        NA
      }
      
    } else stop("wrong type")
    
  })
  
  res=do.call(cbind, out)
  if (type==1) {
    # the first row is n.dean
    boot.n.dean=res[1,]
    res=res[-1,]
  }
  res=res[,!is.na(res[1,])] # remove NA's
  if (verbose) str(res)
  
  # restore rng state 
  assign(".Random.seed", save.seed, .GlobalEnv)    
  
  if (ci.type=="quantile") {
    ci.band=t(apply(res, 1, function(x) quantile(x, c(.025,.975), na.rm=T)))
  } else {
    stop("only quantile bootstrap CI supported for now")
  }
  
  ret = list(marker=if(type==3) names(prob) else ss, prob=prob, boot=res, lb=ci.band[,1], ub=ci.band[,2], if(type==1) n.dean=c(n.dean, boot.n.dean))   
  if (type==1 & !comp.risk) names(ret)[length(ret)]="n.dean" # this is necessary because when using if, that element won't have a name
  ret  
}    
