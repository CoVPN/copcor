# competing risk is handled transparently through form.0, e.g., within marginalized.risk.svycoxph.boot
# multiple imputation is hardcoded by TRIAL, e.g. within marginalized.risk.svycoxph.boot


cor_coxph_risk_bootstrap = function(
  form.0,
  dat,
  fname.suffix, #used in the file names to save results
  save.results.to,
  config,
  config.cor,
  tfinal.tpeak=NULL,
  
  markers,
  
  run.Sgts=T, # whether to get risk conditional on continuous S>=s
  verbose=FALSE
) {
  
cat("Inside function cor_coxph_risk_bootstrap ...............\n")
  
comp.risk=is.list(form.0) # competing risk
  
tpeak=config.cor$tpeak

numCores <- unname(ifelse(Sys.info()["sysname"] == "Windows", 1, min(20, config$num_boot_replicates, future::availableCores())))

B=config$num_boot_replicates

# get tfinal.tpeak
if (is.null(tfinal.tpeak)) tfinal.tpeak=config.cor$tfinal.tpeak
if (is.null(tfinal.tpeak)) stop("tfinal.tpeak should be passed in or in config.cor")

myprint(fname.suffix, tfinal.tpeak, B, numCores, comp.risk, run.Sgts)


###################################################################################################
cat("bootstrap vaccine arm risk, conditional on continuous S=s\n")

fname = paste0(save.results.to, "risks.all.1_", fname.suffix, ".Rdata")
myprint(fname)
  
if(!file.exists(fname)) {    
  risks.all.1=lapply(markers, function (a) {
    if(verbose) myprint(a)
    marginalized.risk.svycoxph.boot(form.0, marker.name=a, type=1, data=dat, t=tfinal.tpeak, B=B, ci.type="quantile", numCores=numCores)
  })
  names(risks.all.1)=markers
  save(risks.all.1, file=fname)
} else {
  load(fname)
}
assign("risks.all.1", risks.all.1, envir = .GlobalEnv) # make it available outside this function


write(ncol(risks.all.1[[1]]$boot), file=paste0(save.results.to, "bootstrap_replicates"))


###################################################################################################
cat("bootstrap vaccine arm, conditional on categorical S\n")

fname = paste0(save.results.to, "risks.all.3_", fname.suffix, ".Rdata")
myprint(fname)

if(!file.exists(fname)) {    
  risks.all.3=lapply(markers, function (a) {
    if(verbose) myprint(a)
    marginalized.risk.svycoxph.boot(form.0, marker.name=a%.%"cat", type=3, data=dat, t=tfinal.tpeak, B=B, ci.type="quantile", numCores=numCores)                
  })    
  names(risks.all.3)=markers
  save(risks.all.3, file=fname)
} else {
  load(fname)
}
assign("risks.all.3", risks.all.3, envir = .GlobalEnv) # make it available outside this function



###################################################################################################
if (run.Sgts) {
  cat("bootstrap vaccine arm risk, conditional on continuous S>=s\n")
  
  fname = paste0(save.results.to, "risks.all.2_", fname.suffix, ".Rdata")
  myprint(fname)
  
  
  if(!file.exists(fname)) {    
    risks.all.2=lapply(markers, function (a) {
      if(verbose) myprint(a)
      marginalized.risk.svycoxph.boot(form.0, marker.name=a, type=2, data=dat, tfinal.tpeak, B=B, ci.type="quantile", numCores=numCores)
    }) 
    names(risks.all.2)=markers
    save(risks.all.2, file=fname)
  } else {
    load(fname)
  }
  assign("risks.all.2", risks.all.2, envir = .GlobalEnv) # make it available outside this function
}


###################################################################################################
# interaction models 

if (!is.null(config$interaction)) {
  if(verbose) cat("Interaction models\n")
  
  if(!file.exists(paste0(save.results.to, "itxn.marginalized.risk.Rdata"))) {    
    cat("make itxn.marginalized.risk\n")
    
    risks.itxn=list()      
    for (ab in config$interaction) {
      tmp=trim(strsplit(ab, " *\\* *")[[1]])
      aold=tmp[1]
      bold=tmp[2]            
      a=paste0("Day",tpeak,aold)
      b=paste0("Day",tpeak,bold)
      
      # idx=2: only use vaccine arm. idx=1 uses placebo data and structural knowledge; it is commented out and moved to the end of the file
      dat.ph1=dat            
      data.ph2=dat.ph1[dat.ph1$ph2==1,]     
      
      # fit the interaction model and save regression results to a table
      f= update(form.0, as.formula(paste0("~.+", a," + ",b," + ",a,":",b)))
      fit=svycoxph(f, design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.ph1)) 
      
      # first treat a as the x axis variable, second treat b as the x axis variable
      for (inner.id in 1:2) {
        if (inner.id == 1) {
          vx=a; vthree=b
        } else {
          vx=b; vthree=a
        }        
        
        # compute risks at three values of vthree
        three.val=wtd.quantile(dat.ph1[[vthree]][dat.ph1$Trt==1], dat.ph1$wt[dat.ph1$Trt==1], c(.15, .5, .85))                    
        # compute risks at a sequence of vx values for each of the three vthree values
        ss=sort(c(
          wtd.quantile(dat.ph1[[vx]], dat.ph1$wt, c(0.025,0.05,0.95,0.975)), # will be included in the table
          seq(min(dat.ph1[[vx]], na.rm=TRUE), max(dat.ph1[[vx]], na.rm=TRUE), length=100) # equally spaced between min and max so that the curves look good
        ))    
        
        #                # special handling code
        #                if (attr(config, "config")=="hvtn705second") {
        #                    if (inner.id == 1) {
        #                        three.val=c(min=-2, wtd.quantile(dat.ph1[[vthree]][dat.ph1$Trt==1], dat.ph1$wt[dat.ph1$Trt==1], c(.5, .9)))
        #                    }
        #                }
        
        
        # estimate marginalized risks, return a matrix
        prob.ls=sapply (three.val, function(val) {
          marginalized.risk.cont.2(fit, marker.name  =vx, data=data.ph2, weights=data.ph2$wt, t=tfinal.tpeak, ss=ss, 
                                   marker.name.2=vthree, s.2=val)
        })
        
        #### bootstrap
        # store the current rng state
        save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
        if (inherits(save.seed,"try-error")) {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }         
        
        seeds=1:B; names(seeds)=seeds
        out=mclapply(seeds, mc.cores = numCores, FUN=function(seed) {   
          seed=seed+560
          if (verbose>=2) myprint(seed)
          
          #                        if (idx==1) {
          #                            # bootstrap vaccine and placebo arm separately
          #                            dat.b = rbind(bootstrap.case.control.samples(subset(dat.ph1, Trt==1), seed, delta.name="EventIndPrimary", strata.name="tps.stratum", ph2.name="ph2"),
          #                                          subset(dat.ph1, Trt==0)[sample.int(nrow(subset(dat.ph1, Trt==0)), r=TRUE),])         
          dat.b = bootstrap.case.control.samples(dat.ph1, seed, delta.name="EventIndPrimary", strata.name="tps.stratum", ph2.name="ph2")
          
          dat.b.ph2=dat.b[dat.b$ph2==1,]
          with(dat.b, table(Wstratum, ph2))     
          
          # inline design object b/c it may also throw an error
          fit.b=try(svycoxph(f, design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.b)))
          
          if (inherits(fit.b,"try-error")) {
            probs=sapply (three.val, function(val) {
              marginalized.risk.cont.2(fit.b, marker.name  =vx, data=dat.b.ph2, weights=dat.b.ph2$wt, t=tfinal.tpeak, ss=ss, 
                                       marker.name.2=vthree, s.2=val)
            })
          } else {
            matrix(NA, length(ss), length(three.val))
          }
          
        })
        
        # restore rng state 
        assign(".Random.seed", save.seed, .GlobalEnv)    
        
        # organize bootstrap results into a list of 3, each element of which is a matrix of n.s by n.seeds
        res.ls=lapply (1:length(three.val), function(i) {
          res=sapply(out, function (x) x[,i])
          res[,!is.na(res[1,])] # remove NA's
        })
        if (verbose>=2) str(res.ls)
        # put lb and ub into matrices
        lb.ls=sapply(res.ls, function (res) t(apply(res, 1, function(x) quantile(x, c(.025)))) )
        ub.ls=sapply(res.ls, function (res) t(apply(res, 1, function(x) quantile(x, c(.975)))) )
        
        risks.itxn[[paste0(vx,vthree)]]=list(marker=ss, prob=prob.ls, boot=res.ls, lb=lb.ls, ub=ub.ls, marker.2=three.val)
      } # end inner.id
      
    }
    save(risks.itxn, file=paste0(save.results.to, "itxn.marginalized.risk_",fname.suffix,".Rdata"))
    
  } else {
    load(paste0(save.results.to, "itxn.marginalized.risk_",fname.suffix,".Rdata"))
  }
  
  assign("risks.itxn", risks.itxn, envir = .GlobalEnv) # make it available outside this function
  
}



}


