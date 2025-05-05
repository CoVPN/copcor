# Output:
# 2 tables for continuous markers 
# 1 table for discrete markers
# 4 forestplot

# possible errors:
# Error in .local(x, i, j, ..., value) : not-yet-implemented 'Matrix[<-' method
#   likely caused by missing values


cor_coxph_coef_1 = function(
  form.0,
  design_or_dat, # either design when tps is T, or data frame otherwise
  fname.suffix, #used in the file names to save results
  save.results.to,
  config,
  config.cor,
  markers,
  markers.names.short,
  
  dat.plac = NULL,
  show.q=TRUE, # whether to show fwer and q values in tables
  
  forestplot.markers=NULL, # make forestplot only for a list of subsets of markers, each member of the list is an array
  for.title="",
  
  run.trichtom=TRUE,
  
  cmp.label = "Placebo",
  verbose=FALSE
  
) {
  
  if(verbose) print("Running cor_coxph_coef_1")
  
  if (is.null(forestplot.markers)) forestplot.markers = list(1:length(markers))
  
  has.plac=!is.null(dat.plac)
  
  # twophase sampling 
  tps=inherits(design_or_dat, "survey.design")
  if (tps) dat = design_or_dat$phase1$full$variables else dat=design_or_dat
  
  tpeak=config.cor$tpeak
  
  myprint(fname.suffix, show.q, tps, has.plac, verbose)
  
  
  ###################################################################################################
  if(verbose) cat("Regression for continuous markers\n")
  
  myboxplot(Day31T4_IFNg_OR_IL2_N_BA.4.5~EventIndPrimary, dat)
  
  fits=list()
  for (a in markers) {
    f = update(form.0, as.formula(paste0("~.+", a)))
    if (tps) {
      fits[[a]]=svycoxph(f, design=design_or_dat) 
    } else {
      fits[[a]]=coxph(f, design_or_dat) 
    }
  }
  
  # scaled marker
  fits.scaled=list()
  for (a in markers) {
    f= update(form.0, as.formula(paste0("~.+scale(", a, ")")))
    if (tps) {
      fits.scaled[[a]]=svycoxph(f, design=design_or_dat) 
    } else {
      fits.scaled[[a]]=coxph(f, design_or_dat) 
    }
  }
  
  # put coxph model coef together to save
  fits.cont.coef.ls = lapply(fits, function (fit) getFixedEf(fit, robust=tps))
  
  natrisk=nrow(dat)
  dat$yy = dat[[as.character(form.0[[2]][[3]])]]
  nevents = sum(dat$yy==1)
  
  # make pretty table
  rows=length(coef(fits[[1]]))
  est=getFormattedSummary(fits, exp=T, robust=tps, rows=rows, type=1)
  ci= getFormattedSummary(fits, exp=T, robust=tps, rows=rows, type=7)
  p=  getFormattedSummary(fits, exp=T, robust=tps, rows=rows, type=10)
  est.scaled=getFormattedSummary(fits.scaled, exp=T, robust=tps, rows=rows, type=1)
  ci.scaled= getFormattedSummary(fits.scaled, exp=T, robust=tps, rows=rows, type=7)
  
  pvals.cont = sapply(fits, function(x) {
    tmp=getFixedEf(x)
    p.val.col=which(startsWith(tolower(colnames(tmp)),"p"))
    tmp[nrow(tmp),p.val.col]
  })
  
  
  # make forest plots
  
  if (!is.list(forestplot.markers)) forestplot.markers=list(forestplot.markers)
  
  for (i in 1:2) {# not scaled and scaled
    
    res=getFormattedSummary(if (i==1) fits else fits.scaled, exp=F, robust=tps, rows=rows, type=0)
    res=t(res)
    # res: est, se, lb, ub, pvalue 
    p.val.col=which(startsWith(tolower(colnames(res)),"p"))
    
    for (iM in 1:length(forestplot.markers)) {
      
      est.ci = rbind(exp(res[forestplot.markers[[iM]], 1]), 
                     exp(res[forestplot.markers[[iM]], 3]), 
                     exp(res[forestplot.markers[[iM]], 4]), 
                     res[forestplot.markers[[iM]], p.val.col]
      )
      colnames(est.ci)=markers.names.short[forestplot.markers[[iM]]]
      
      # inf values break theforestplot
      est.ci[abs(est.ci)>100]=sign(est.ci[abs(est.ci)>100])*100
      
      # make sure point lb < ub and lb is not too close to 0, which, when log transformed, lead to errors
      kp = est.ci[2,]>est.ci[3,] | est.ci[2,]<1e-10
      est.ci=est.ci[,!kp,drop=F]

      # make two versions, one log and one antilog
      
      fig.height = 4*ncol(est.ci)/13
      if(ncol(est.ci)<=8) fig.height=2*fig.height
      if(ncol(est.ci)<=4) fig.height=2*fig.height
      
      mypdf(onefile=F, width=10,height=fig.height, file=paste0(save.results.to, "hr_forest_", ifelse(i==1,"","scaled_"), fname.suffix, if (iM>1) iM)) 
      theforestplot(point.estimates=est.ci[1,], lower.bounds=est.ci[2,], upper.bounds=est.ci[3,], group=colnames(est.ci), 
                    nEvents=rep(NA, ncol(est.ci)), # as table.labels below shows, we are not showing nevents
                    p.values=formatDouble(est.ci[4,], 3, remove.leading0=F), 
                    decimal.places=2, graphwidth=unit(120, "mm"), fontsize=1.2, 
                    table.labels = c("", "  HR (95% CI)",""), 
                    title=for.title, 
                    xlog=F,
                    x.ticks = get.forestplot.ticks(est.ci, forestplot.xlog=F)  # controls the limit
      )
      dev.off()
      
      mypdf(onefile=F, width=10,height=fig.height, file=paste0(save.results.to, "hr_forest_log_", ifelse(i==1,"","scaled_"), fname.suffix, if (iM>1) iM)) 
      theforestplot(point.estimates=est.ci[1,], lower.bounds=est.ci[2,], upper.bounds=est.ci[3,], group=colnames(est.ci), 
                    nEvents=rep(NA, ncol(est.ci)), # as table.labels below shows, we are not showing nevents, 
                    p.values=formatDouble(est.ci[4,], 3, remove.leading0=F), 
                    decimal.places=2, graphwidth=unit(120, "mm"), fontsize=1.2, 
                    table.labels = c("", "  HR (95% CI)",""), 
                    title=for.title, 
                    xlog=T,
                    x.ticks = get.forestplot.ticks(est.ci, forestplot.xlog=T)  # controls the limit
      )
      dev.off()
      
    }
  }
  
  ###################################################################################################
  if (run.trichtom) {
      
    if(verbose) cat("regression for trichotomized markers\n")

    marker.levels = sapply(markers, function(a) length(table(dat[[a%.%"cat"]]))); marker.levels
    
    fits.tri=list()
    for (a in markers) {
      if(verbose>=2) myprint(a)
      f = update(form.0, as.formula(paste0("~.+", a, "cat")))
      if (tps) {
        fits.tri[[a]]=svycoxph(f, design=design_or_dat) 
      } else {
        fits.tri[[a]]=coxph(f, design_or_dat) 
      }
    }
    fits.tri=fits.tri
    
    fits.tri.coef.ls= lapply(fits.tri, function (fit) getFixedEf(fit, robust=tps))
    
    
    # get generalized Wald p values
    p.cov = length(strsplit(config$covariates, "\\+")[[1]])-1
    overall.p.tri=sapply(fits.tri, function(fit) {
      rows = (1+p.cov):length(coef(fit))
      if (length(fit)==1) NA else {
        stat=coef(fit)[rows] %*% solve(vcov(fit,robust=tps)[rows,rows]) %*% coef(fit)[rows]
        pchisq(stat, length(rows), lower.tail = FALSE)
      }
    })
    #
    overall.p.0=formatDouble(unlist(lapply(markers, function (a) c(overall.p.tri[a], rep(NA,marker.levels[a]-1)))), digits=3, remove.leading0 = F);   
    overall.p.0=sub("0.000","<0.001",overall.p.0)
    
  }
  
  
  
  ###################################################################################################
  if(show.q & verbose) cat("multitesting adjustment for continuous and trichotomized markers together\n")
  
  # If primary_assays is not defined in config, multitesting adjustment is over all assays. 
  # If primary_assays defined, multitesting adjustment is only over this subset. If this set is empty, then no multitesting adjustment is done
  
  if (show.q) {
  
  p.unadj.1 = c(cont=pvals.cont, tri=overall.p.tri) 
  # the markers to perform multitestign are defined in config$primary_assays
  p.unadj=c()
  if (!is.null(config$primary_assays)) {
    if (length(config$primary_assays)>0) {
      p.unadj = p.unadj.1[c("cont.Day"%.%tpeak%.%config$primary_assays, "tri.Day"%.%tpeak%.%config$primary_assays)]
    } 
  }
  if (TRIAL=="prevent19") {
    # bindSpike tertiary has no cases in the upper tertile, cannot do P value
    p.unadj = p.unadj[startsWith(names(p.unadj), "cont."), drop=F]
  }

  
  if (length(p.unadj)>1) {
    if (verbose) cat("doing Westfall and Young no select markers\n")
    
    #### Westfall and Young permutation-based adjustment
    if(!file.exists(paste0(save.results.to, "pvals.perm.",fname.suffix,".Rdata"))) {
      
      dat.ph2 = design_or_dat$phase1$sample$variables
      design_or_dat.perm=design_or_dat
      #design_or_dat.perm$phase1$full$variables
      
      #    # if want to only do multitesting when liveneutmn50 is included
      #    if (!"liveneutmn50" %in% assays) numPerm=5
      
      # TODO: there is no need to permute markers
      numPerm = config$num_perm_replicates
      numCores <- unname(ifelse(Sys.info()["sysname"] == "Windows",
                                1, 
                                min(20, config$num_boot_replicates, future::availableCores())))
      
      out=mclapply(1:numPerm, mc.cores = numCores, FUN=function(seed) {   
        # store the current rng state 
        save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
        if (inherits(save.seed,"try-error")) {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }          
        set.seed(seed)        
        
        # permute markers in design_or_dat.perm
        new.idx=sample(1:nrow(dat.ph2))
        tmp=dat.ph2
        for (a in markers) {
          tmp[[a]]=tmp[[a]][new.idx]
          tmp[[a%.%"cat"]]=tmp[[a%.%"cat"]][new.idx]
        }
        design_or_dat.perm$phase1$sample$variables = tmp
        
        # rename markers so that out has the proper names. this is only effective within permutation
        names(markers)=markers
        out=c(
          cont=sapply (markers, function(a) {
            f= update(form.0, as.formula(paste0("~.+", a)))
            if (tps) {
              fit=try(svycoxph(f, design=design_or_dat.perm), silent=T)
              if (inherits(fit,"try-error")) fit=NA
              
            } else {
              # TODO
            }
            if (length(fit)==1) NA else mylast(c(getFixedEf(fit)))
          })        
          ,    
          tri=sapply (markers, function(a) {
            f= update(form.0, as.formula(paste0("~.+", a, "cat")))
            if (tps) {
              fit=try(svycoxph(f, design=design_or_dat.perm), silent=T)
              if (inherits(fit,"try-error")) fit=NA
              
            } else {
              # TODO
            }
            if (length(fit)==1) NA else mylast(c(getFixedEf(fit)))
          })
        )
        
        # restore rng state 
        assign(".Random.seed", save.seed, .GlobalEnv)    
        
        out
      })
      pvals.perm=do.call(rbind, out)
      save(pvals.perm, file=paste0(save.results.to, "pvals.perm."%.%fname.suffix%.%".Rdata"))
      
    } else {
      load(file=paste0(save.results.to, "pvals.perm."%.%fname.suffix%.%".Rdata"))
    }
    # save number of permutation replicates
    write(nrow(pvals.perm), file=paste0(save.results.to, "permutation_replicates_"%.%fname.suffix))
    
    
    if(any(is.na(p.unadj))) {
      pvals.adj = cbind(p.unadj=p.unadj, p.FWER=NA, p.FDR=NA)
    } else {
      pvals.adj = p.adj.perm (p.unadj, pvals.perm[,names(p.unadj)], alpha=1)  
    }
    if(verbose) print(pvals.adj)
    
    
  } else {
    if (verbose) cat("doing Holm and FDR adjustment on all markers\n")
    
    pvals.adj.fdr=p.adjust(p.unadj.1, method="fdr")
    pvals.adj.hol=p.adjust(p.unadj.1, method="holm")
    
    pvals.adj=cbind(p.unadj.1, p.FWER=pvals.adj.hol, p.FDR=pvals.adj.fdr)
    write(NA, file=paste0(save.results.to, "permutation_replicates_"%.%fname.suffix))     # so the rmd file can compile
  }
  
  # not all markers were multitesting adjusted  
  pvals.adj = cbind(p.unadj=p.unadj.1, pvals.adj[match(names(p.unadj.1), rownames(pvals.adj)),2:3, drop=F])
  
  if (TRIAL=="prevent19" | TRIAL=='prevent19_stage2') {
    # bindSpike tertiary has no cases in the upper tertile, cannot do P value
    # this code somehow works even though there are also RBD and ID50
    pvals.adj=rbind(pvals.adj, tri.Day35bindSpike=c(NA,NA,NA))
  }
  
  p.1=formatDouble(pvals.adj["cont."%.%names(pvals.cont),"p.FWER"], 3, remove.leading0=F); p.1=sub("0.000","<0.001",p.1)
  p.2=formatDouble(pvals.adj["cont."%.%names(pvals.cont),"p.FDR" ], 3, remove.leading0=F); p.2=sub("0.000","<0.001",p.2)
  
  }
  
  
  ###################################################################################################
  # make continuous markers table
  
  tab.1=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), t(est), t(ci), t(p), if(show.q) p.2, if(show.q) p.1)
  rownames(tab.1)=markers.names.short
  tab.1
  
  if (show.q) {
    header=paste0("\\hline\n 
         \\multicolumn{1}{l}{", '', "} & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{HR per 10-fold incr.}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{q-value}   & \\multicolumn{1}{c}{FWER} \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{} & \\multicolumn{1}{c}{***} & \\multicolumn{1}{c}{} \\\\ 
         \\hline\n 
    ")
  } else {
    header=paste0("\\hline\n 
         \\multicolumn{1}{l}{", '', "} & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{HR per 10-fold incr.}                     & \\multicolumn{1}{c}{P-value}    \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{}  \\\\ 
         \\hline\n 
    ")
  }
  
  mytex(tab.1, file.name="CoR_univariable_svycoxph_pretty_"%.%fname.suffix, align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
        col.headers=header,
        longtable=T, 
        label=paste0("tab:CoR_univariable_svycoxph_pretty"), 
        caption.placement = "top", 
        caption=paste0("Inference for Day ", tpeak, " antibody marker covariate-adjusted correlates of risk of ", config.cor$txt.endpoint, " in the ", escape(fname.suffix), " group: Hazard ratios per 10-fold increment in the marker*")
  )
  tab.cont=tab.1
  
  tab.1.nop12=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), t(est), t(ci), t(p))
  rownames(tab.1.nop12)=markers.names.short
  
  # scaled markers
  tab.1.scaled=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), t(est.scaled), t(ci.scaled), t(p), if(show.q) p.2, if(show.q) p.1)
  rownames(tab.1.scaled)=markers.names.short
  tab.1.scaled
  
  if (show.q) {
    header=paste0("\\hline\n 
         \\multicolumn{1}{l}{} & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{HR per SD incr.}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{q-value}   & \\multicolumn{1}{c}{FWER} \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{} & \\multicolumn{1}{c}{***} & \\multicolumn{1}{c}{} \\\\ 
         \\hline\n 
    ")
  } else {
    header=paste0("\\hline\n 
         \\multicolumn{1}{l}{} & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{HR per SD incr.}                     & \\multicolumn{1}{c}{P-value}   \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{}  \\\\ 
         \\hline\n 
    ")
  }
  
  mytex(tab.1.scaled, file.name="CoR_univariable_svycoxph_pretty_scaled_"%.%fname.suffix, align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
        col.headers=header,
        longtable=T, 
        label=paste0("tab:CoR_univariable_svycoxph_pretty_scaled"), 
        caption.placement = "top", 
        caption=paste0("Inference for Day ", tpeak, " antibody marker covariate-adjusted correlates of risk of ", config.cor$txt.endpoint, " in the ", escape(fname.suffix), " group: Hazard ratios per SD increment in the marker*")
  )
  tab.cont.scaled=tab.1.scaled
  
  
  
  ###################################################################################################
  # make trichotomized markers table
  
  if (run.trichtom) {
  
  if(show.q) {
    overall.p.1=formatDouble(pvals.adj["tri."%.%names(pvals.cont),"p.FWER"], 3, remove.leading0=F);   overall.p.1=sub("0.000","<0.001",overall.p.1)
    overall.p.2=formatDouble(pvals.adj["tri."%.%names(pvals.cont),"p.FDR" ], 3, remove.leading0=F);   overall.p.2=sub("0.000","<0.001",overall.p.2)
    
    # add space
    overall.p.1=formatDouble(unlist(lapply(markers, function (a) c(pvals.adj["tri."%.%a,"p.FWER"], rep(NA,marker.levels[a]-1)))), digits=3, remove.leading0 = F);   
    overall.p.1=sub("0.000","<0.001",overall.p.1)
    overall.p.2=formatDouble(unlist(lapply(markers, function (a) c(pvals.adj["tri."%.%a,"p.FDR"],  rep(NA,marker.levels[a]-1)))), digits=3, remove.leading0 = F);   
    overall.p.2=sub("0.000","<0.001",overall.p.2)
  }
  
  # if "Delta"%.%tpeak%.%"overB" is included, nevents have a problem because some markers may have only two category in the cases
  
  # n cases and n at risk
  natrisk = round(unlist(lapply (markers%.%"cat", function(a) aggregate(dat[dat$ph2==1,] [["wt"]], dat[dat$ph2==1,] [a], sum, na.rm=T, drop=F)[,2] )))
  nevents = round(unlist(lapply (markers%.%"cat", function(a) {
    aggregate(dat[dat$ph2==1 & dat$yy==1,] [["wt"]], dat[dat$ph2==1 & dat$yy==1,] [a], sum, na.rm=T, drop=F)[,2] 
  }
  )))
  natrisk[is.na(natrisk)]=0
  nevents[is.na(nevents)]=0
  colSums(matrix(natrisk, nrow=3))
  
  # regression parameters
  est=unlist(lapply(markers, function (a) 
    if(length(fits.tri[[a]])==1) rep(NA,marker.levels[a]) else {
      c(1,     getFormattedSummary(list(fits.tri[[a]]), exp=T, robust=tps, rows=p.cov+1:(marker.levels[[a]]-1), type=1))
    }
  ))
  ci= unlist(lapply(markers, function (a) 
    if(length(fits.tri[[a]])==1) rep(NA,marker.levels[[a]]) else {
      c("N/A", getFormattedSummary(list(fits.tri[[a]]), exp=T, robust=tps, rows=p.cov+1:(marker.levels[[a]]-1), type=7)) 
    }
  ))
  p=  unlist(lapply(markers, function (a) 
    if(length(fits.tri[[a]])==1) rep(NA,marker.levels[[a]]) else {
      c("N/A", getFormattedSummary(list(fits.tri[[a]]), exp=T, robust=tps, rows=p.cov+1:(marker.levels[[a]]-1), type=10)) 
    }
  ))
  
  
  tab=cbind(
    unlist(lapply(markers, function (a) c("Lower", if (marker.levels[a]==3) "Middle", "Upper"))),
    paste0(nevents, "/", format(natrisk, big.mark=",",digit=0, scientific=F)), 
    formatDouble(nevents/natrisk, digits=4, remove.leading0=F),
    est, ci, p, 
    overall.p.0, if(show.q) overall.p.2, if(show.q) overall.p.1
  )
  # str(nevents); str(natrisk)
  rownames(tab)=unlist(lapply(markers, function (a) c(markers.names.short[a],rep("",marker.levels[a]-1)))) 
  tab.cat=tab[1:(nrow(tab)),]
  #cond.plac=dat.plac[[config.cor$EventTimePrimary]]<=tfinal.tpeak # not used anymore
  
  
  if(show.q) {
    header=paste0("\\hline\n 
         \\multicolumn{1}{l}{", '', "} & \\multicolumn{1}{c}{Tertile}   & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{1}{c}{Attack}   & \\multicolumn{2}{c}{Haz. Ratio}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{Overall P-}      & \\multicolumn{1}{c}{Overall q-}   & \\multicolumn{1}{c}{Overall} \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{}          & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{rate}   & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{} & \\multicolumn{1}{c}{value***} & \\multicolumn{1}{c}{value $\\dagger$} & \\multicolumn{1}{c}{FWER} \\\\ 
         \\hline\n 
    ")
  } else {
    header=paste0("\\hline\n 
         \\multicolumn{1}{l}{", '', "} & \\multicolumn{1}{c}{Tertile}   & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{1}{c}{Attack}   & \\multicolumn{2}{c}{Haz. Ratio}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{Overall P-}     \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{}          & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{rate}   & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{} & \\multicolumn{1}{c}{value***} \\\\ 
         \\hline\n 
    ")
  }
  
  if (has.plac) {
    add.to.row=list(list(nrow(tab)), # insert at the beginning of table, and at the end of, say, the first table
                    c(paste0(" \n \\multicolumn{8}{l}{} \\\\ \n", 
                             "\n \\multicolumn{2}{l}{", cmp.label, "} & ", 
                             paste0(sum(dat.plac$yy), "/", format(nrow(dat.plac), big.mark=",")), "&",  
                             formatDouble(sum(dat.plac$yy)/nrow(dat.plac), digits=4, remove.leading0=F), "&",  
                             "\\multicolumn{4}{l}{}  \\\\ \n")
                      #"\\hline\n \\multicolumn{4}{l}{Standard Deviation 1 mcg/mL}\\\\ \n"
                    )
    )
  } else {
    add.to.row = NULL
  }
  
  # use longtable because this table could be long, e.g. in hvtn705second
  mytex(tab[1:(nrow(tab)),], file.name="CoR_univariable_svycoxph_cat_pretty_"%.%fname.suffix, align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
        col.headers=header,        
        add.to.row=add.to.row,
        longtable=T, 
        label=paste0("tab:CoR_univariable_svycoxph_cat_pretty_", fname.suffix), 
        caption.placement = "top", 
        caption=paste0("Inference for Day ", tpeak, " antibody marker covariate-adjusted correlates of risk of ", config.cor$txt.endpoint, " in the ", escape(fname.suffix), " group: Hazard ratios for Middle vs. Upper tertile vs. Lower tertile*")
  )
  
  
  # save two subjects for collate
  if (has.plac) {
    save.s.1=paste0(sum(dat.plac$yy), "/", format(nrow(dat.plac), big.mark=","))
    save.s.2=formatDouble(sum(dat.plac$yy)/nrow(dat.plac), digits=4, remove.leading0=F)
  }
  
  
  tab.nop12=cbind(
    unlist(lapply(markers, function (a) c("Lower", if (marker.levels[a]==3) "Middle", "Upper"))),
    paste0(nevents, "/", format(natrisk, big.mark=",",digit=0, scientific=F)), 
    formatDouble(nevents/natrisk, digits=4, remove.leading0=F),
    est, ci, p, overall.p.0
  )
  rownames(tab.nop12)=unlist(lapply(markers, function (a) c(markers.names.short[a],rep("",marker.levels[a]-1)))) 
  
  }
  
  
  ###################################################################################################
  # multivariate_assays models
  
  if (!is.null(config$multivariate_assays)) {
    if(verbose) cat("Multiple regression\n")
    
    for (a in config$multivariate_assays) {
      for (i in 1:2) {
        # 1: per SD; 2: per 10-fold
        a.tmp=a
        aa=trim(strsplit(a, "\\+")[[1]])
        for (x in aa[!contain(aa, "\\*")]) {
          # replace every variable with Day210x, with or without scale
          a.tmp=gsub(x, paste0(if(i==1) "scale","(Day",tpeak,x,")"), a.tmp) 
        }
        f= update(form.0, as.formula(paste0("~.+", a.tmp)))
        if (tps) {
          fit=svycoxph(f, design=design_or_dat) 
        } else {
          fit=coxph(f, dat) 
        }
        var.ind=length(coef(fit)) - length(aa):1 + 1
        
        fits=list(fit)
        est=getFormattedSummary(fits, exp=T, robust=tps, rows=var.ind, type=1)
        ci= getFormattedSummary(fits, exp=T, robust=tps, rows=var.ind, type=7)
        est = paste0(est, " ", ci)
        p=  getFormattedSummary(fits, exp=T, robust=tps, rows=var.ind, type=10)
        
        #generalized Wald test for whether the set of markers has any correlation (rejecting the complete null)
        stat=coef(fit)[var.ind] %*% solve(vcov(fit)[var.ind,var.ind]) %*% coef(fit)[var.ind] 
        p.gwald=pchisq(stat, length(var.ind), lower.tail = FALSE)
        
        tab=cbind(est, p)
        # tmp=match(aa, colnames(labels.axis))
        # tmp[is.na(tmp)]=1 # otherwise, labels.axis["Day"%.%tpeak, tmp] would throw an error when tmp is NA
        # rownames(tab)=ifelse(aa %in% colnames(labels.axis), labels.axis["Day"%.%tpeak, tmp], aa)
        colnames(tab)=c(paste0("HR per ",ifelse(i==1,"sd","10 fold")," incr."), "P value")
        tab
        tab=rbind(tab, "Generalized Wald Test"=c("", formatDouble(p.gwald,3, remove.leading0 = F)))
        
        mytex(tab, file.name=paste0("CoR_multivariable_svycoxph_pretty", match(a, config$multivariate_assays), if(i==2) "_per10fold", fname.suffix), align="c", include.colnames = T, save2input.only=T, 
              input.foldername=save.results.to)
      }
    }
    
  }
  
  
  
  ###################################################################################################
  # additional_models
  
  if (!is.null(config$additional_models)) {
    if(verbose) cat("Additional models\n")
    
    for (a in config$additional_models) {
      tmp=gsub("tpeak",tpeak,a)
      f= update(Surv(EventTimePrimary, EventIndPrimary) ~1, as.formula(paste0("~.+", tmp)))
      if (tps) {
        fits[[a]]=svycoxph(f, design=design_or_dat) 
      } else {
        fits[[a]]=coxph(f, dat) 
      }
      
      fits=list(fit)
      est=getFormattedSummary(fits, exp=T, robust=tps, type=1)
      ci= getFormattedSummary(fits, exp=T, robust=tps, type=7)
      est = paste0(est, " ", ci)
      p=  getFormattedSummary(fits, exp=T, robust=tps, type=10)
      
      tab=cbind(est, p)
      colnames(tab)=c("HR", "P value")
      tab
      
      mytex(tab, file.name=paste0("CoR_add_models", match(a, config$additional_models)), align="c", include.colnames = T, save2input.only=T, 
            input.foldername=save.results.to)
    }
    
  }
  
  
  if (attr(config,"config")=="janssen_pooled_EUA") {
    f=Surv(EventTimePrimary, EventIndPrimary) ~ risk_score + as.factor(Region) * Day29pseudoneutid50    
    if (tps) {
      fits[[a]]=svycoxph(f, design=design_or_dat) 
    } else {
      fits[[a]]=coxph(f, dat) 
    }
    var.ind=5:6
    
    fits=list(fit)
    est=getFormattedSummary(fits, exp=T, robust=tps, rows=1:6, type=1)
    ci= getFormattedSummary(fits, exp=T, robust=tps, rows=1:6, type=7)
    est = paste0(est, " ", ci)
    p=  getFormattedSummary(fits, exp=T, robust=tps, rows=1:6, type=10)
    
    #generalized Wald test for whether the set of markers has any correlation (rejecting the complete null)
    stat=coef(fit)[var.ind] %*% solve(vcov(fit)[var.ind,var.ind]) %*% coef(fit)[var.ind] 
    p.gwald=pchisq(stat, length(var.ind), lower.tail = FALSE)
    
    tab=cbind(est, p)
    colnames(tab)=c("HR", "P value")
    tab=rbind(tab, "Generalized Wald Test for Itxn"=c("", formatDouble(p.gwald,3, remove.leading0 = F)))
    tab
    
    mytex(tab, file.name="CoR_region_itxn", align="c", include.colnames = T, save2input.only=T, 
          input.foldername=save.results.to)
    
  }
  
  
  ###################################################################################################
  # interaction_models
  
  if (!is.null(config$interaction)) {
    if(verbose) cat("Interaction models Cox models\n")    
    itxn.pvals=c()      
    for (ab in config$interaction) {
      tmp=trim(strsplit(ab, " *\\* *")[[1]])
      aold=tmp[1]
      bold=tmp[2]            
      a=paste0("Day",tpeak,aold)
      b=paste0("Day",tpeak,bold)
      
      # fit the interaction model and save regression results to a table
      f= update(form.0, as.formula(paste0("~.+", a," + ",b," + ",a,":",b)))
      if (tps) {
        fit=svycoxph(f, design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat)) 
      } else {
        fits[[a]]=coxph(f, dat) 
      }
      
      fits=list(fit)
      est=getFormattedSummary(fits, exp=T, robust=tps, type=1)
      ci= getFormattedSummary(fits, exp=T, robust=tps, type=7)
      est = paste0(est, " ", ci)
      p=  getFormattedSummary(fits, exp=T, robust=tps, type=10)
      # generalized Wald test for whether the set of markers has any correlation (rejecting the complete null)
      var.ind=length(coef(fit))-2:0
      stat=coef(fit)[var.ind] %*% solve(vcov(fit)[var.ind,var.ind]) %*% coef(fit)[var.ind] 
      p.gwald=pchisq(stat, length(var.ind), lower.tail = FALSE)
      # put together the table
      tab=cbind(est, p)
      colnames(tab)=c("HR", "P value")
      tab=rbind(tab, "Generalized Wald Test for Markers"=c("", formatDouble(p.gwald,3, remove.leading0 = F))); tab
      mytex(tab, file.name=paste0("CoR_itxn_",aold,"_",bold), align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to)
      
      itxn.pvals=c(itxn.pvals, mylast(getFixedEf(fit)[,"p.value"]))
    }    
    
    names(itxn.pvals)=config$interaction
    itxn.pvals=itxn.pvals[!contain(config$interaction, "ICS4AnyEnv")] # remove the ones with ICS4AnyEnv
    itx.pvals.adj.fdr=p.adjust(itxn.pvals, method="fdr")
    itx.pvals.adj.hol=p.adjust(itxn.pvals, method="holm")
    tab=cbind(itxn.pvals, itx.pvals.adj.hol, itx.pvals.adj.fdr)
    colnames(tab)=c("interaction P value", "FWER", "FDR")
    mytex(tab, file.name="CoR_itxn_multitesting", align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to)
    
  }
  
  
  
  ###################################################################################################
  if(run.trichtom) {
    save(fits.cont.coef.ls, fits.tri.coef.ls, file=paste0(save.results.to, paste0("coxph_fits_", fname.suffix, ".Rdata")))
    # save.s.1, save.s.2
    save (tab.cont, tab.cat, tab.cont.scaled, file=paste0(save.results.to, paste0("coxph_slopes_", fname.suffix, ".Rdata")))
  } else {
    save(fits.cont.coef.ls, file=paste0(save.results.to, paste0("coxph_fits_", fname.suffix, ".Rdata")))
    # save.s.1, save.s.2
    save (tab.cont, tab.cont.scaled, file=paste0(save.results.to, paste0("coxph_slopes_", fname.suffix, ".Rdata")))
  }
  
  
}

