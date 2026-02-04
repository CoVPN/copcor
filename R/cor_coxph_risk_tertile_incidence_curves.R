# makes tertile incidence curves figure and logclog figure

# competing risk is handled through form.0, and multiple imputation is hardcoded by TRIAL

# janssen_partA_VL 
#   outcome is multiple imputed; depends on variant, which is defined globally
#   marker is multiple imputed
#   ph2 and wt depends on marker

# vat08_combined  
#   outcome is multiple imputed


cor_coxph_risk_tertile_incidence_curves = function(
  form.0,
  dat,
  fname.suffix,
  save.results.to,
  config,
  config.cor,
  tfinal.tpeak,

  markers,
  markers.names.short,
  markers.names.long,
  marker.cutpoints,
  assay_metadata,
  
  dat.plac = NULL,
  for.title="",
  
  trt.label = "Vaccine",
  cmp.label = "Placebo",
  verbose=FALSE
) {
  


if(verbose) print(paste0("Running cor_coxph_risk_tertile_incidence_curves: ", fname.suffix))

# load overall.p.tri
tmp=paste0(save.results.to, paste0("coxph_overall_p_tri_", fname.suffix, ".Rdata"))
if (file.exists(tmp)) load(tmp)  
    
#### define mi and comp.risk
  
mi = TRIAL %in% c("janssen_partA_VL", "vat08_combined")
  
comp.risk = is.list(form.0)

set.mi.data = function(dat, config.cor, imp, marker.name=NULL) {
  if (TRIAL %in% c("janssen_partA_VL")) {
    dat$EventIndOfInterest = ifelse(dat$EventIndPrimary==1 & dat[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
    dat$EventIndCompeting  = ifelse(dat$EventIndPrimary==1 & dat[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
    if(!is.null(marker.name)) dat[[marker.name]] = dat[[substr(marker.name, 1, nchar(marker.name)-3)%.%"_"%.%imp%.%"cat"]]
    
  } else  if (TRIAL == c("vat08_combined")) {
    dat$EventIndOfInterest  = dat[[config.cor$EventIndPrimary  %.% imp]]
    dat$EventTimeOfInterest = dat[[config.cor$EventTimePrimary %.% imp]]
    if (endsWith(COR,"st1.nAb.batch0and1")) {
      if(!is.null(marker.name)) dat[[marker.name]] = dat[[substr(marker.name, 1, nchar(marker.name)-3)%.%"_"%.%imp%.%"cat"]]
    }
  }
  dat
}
  
if(is.null(tfinal.tpeak)) tfinal.tpeak = config.cor$tfinal.tpeak
if (is.null(tfinal.tpeak)) stop("missing tfinal.tpeak")
  
has.plac=!is.null(dat.plac)
  
myprint(comp.risk, has.plac)
  
# form.s is the null model, i.e. rhs is ~1
form.s=as.formula(deparse((if(comp.risk) form.0[[1]] else form.0)[[2]])%.%"~1")

# # compute prevalence
# # note that these do not have CI
# prev.vacc = if (mi) {
#   mean(sapply(1:10, function(imp) {
#     dat = set.mi.data(dat, config.cor, imp)
#     get.marginalized.risk.no.marker(form.0, dat, tfinal.tpeak)
#   }))
# } else {
#   get.marginalized.risk.no.marker(form.0, dat, tfinal.tpeak)
# }


.mfrow <- c(1, 1)

# origin of followup days, may be different from tpeak
tpeak1 = config.cor$torigin
if (is.null(tpeak1)) tpeak1 = tpeak

tpeak = config.cor$tpeak
tpeaklag = config.cor$tpeaklag

study_name = config$study_name

assays=assay_metadata$assay
llox_labels=assay_metadata$llox_label; names(llox_labels)=assays
lloqs=assay_metadata$lloq; names(lloqs)=assays
uloqs=assay_metadata$uloq; names(uloqs)=assays
lods=assay_metadata$lod; names(lods)=assays
lloxs=ifelse(llox_labels=="lloq", lloqs, lods)
lloxs=ifelse(llox_labels=="pos", assay_metadata$pos.cutoff, lloxs)


if (TRIAL=="janssen_partA_VL") {
  # This is needed because we don't do multiple imputation when computing hist
  for(a in markers) {
    if (!is.null(dat[[a%.%"_1"]])) dat[[a]] = dat[[a%.%"_1"]]
  }
}


###################################################################################################
# some parameters

{
  s2="85%"; s1="15%" # these two reference quantiles are used in the next two blocks of code
  RRud=RReu=2
  bias.factor=bias.factor(RRud, RReu)
  
  # load ylims.cor if exists
  # by changing which file to load, we can enforce the same ylim for different sets of tables/figures
  tmp=paste0(here::here(), paste0(save.results.to, "/ylims.cor_", fname.suffix, ".Rdata"))
  if (file.exists(tmp)) {
    load(tmp)
    create.ylims.cor=F
  } else {
    ylims.cor=list()
    ylims.cor[[1]]=list(NA,NA)
    ylims.cor[[2]]=list(NA,NA)
    create.ylims.cor=T
  }
  
  report.ve.levels=c(.65,.9,.95)
  digits.risk=4
  
  # # yy is used to chose the boundary for showing >=s risk curves 
  # if (TRIAL=="janssen_partA_VL") {
  #   dat$yy = ifelse(dat$EventIndPrimary==1 & dat[["seq1.variant.hotdeck1"]]==variant, 1, 0)    
  # } else {
  #   dat$yy = ifelse(dat$EventIndPrimary==1, 1, 0)
  # }
}



################################################################################
cat("plot marginalized risk curves over time for trichotomized markers\n")

# no bootstrap

data.ph2 <- dat[dat$ph2==1,]
risks.all.ter=list()

# compute risks.all.ter
for (a in markers) {        
  
  marker.name=a%.%"cat"    
  ss=unique(dat[[marker.name]]); ss=sort(ss[!is.na(ss)])
  if(length(ss)==3) {
    names(ss)=c("low","med","high")
  } else {
    names(ss)=c("low","high")
  }
  
  # make risks.all.ter
  if (comp.risk) {
    
    f1=lapply(form.0, function(x) update(x, as.formula(paste0("~.+",marker.name))))
    
    if (TRIAL=="janssen_partA_VL") { 
      
      # ancestral markers have different weights and ph2 indicators
      if (a %in% c("Day29bindSpike","Day29pseudoneutid50")) {
        data.ph2$ph2 = data.ph2$ph2.D29
        data.ph2$wt = data.ph2$wt.D29
      } else {
        data.ph2$ph2 = data.ph2$ph2.D29variant
        data.ph2$wt = data.ph2$wt.D29variant
      }
      
      # multiple imputation
      out=lapply(1:10, function(imp) {
        # impute outcome
        data.ph2 = set.mi.data(data.ph2, config.cor, imp, marker.name)

        newdata=data.ph2
        out=lapply(ss, function(s) {
          newdata[[marker.name]]=s
          risks = pcr2(f1, data.ph2, tfinal.tpeak, weights=data.ph2$wt, newdata=newdata)
          cbind(t=attr(risks,"time"), 
                cumulative=apply(attr(risks,"cumulative"), 1, weighted.mean, weights=data.ph2$wt))
        })
        cbind(out[[1]], out[[2]][,"cumulative"], if(length(out)==3) out[[3]][,"cumulative"])
      })
      
      all.t=lapply(out, function(x) x[,"t"])
      common.t = Reduce(intersect, all.t)
      risks = sapply(out, simplify="array", function(x) x[x[,"t"] %in% common.t,-1,drop=F])
      risks.all.ter[[a]]=list(time=common.t, risk=apply(risks, 1:2, mean, na.rm=T))
      
    } else {
      
      out=lapply(ss, function(s) {
        newdata = data.ph2
        newdata[[marker.name]]=s
        risks = pcr2(f1, data.ph2, tfinal.tpeak, weights=data.ph2$wt, newdata=newdata)
        cbind(t=attr(risks,"time"), 
              cumulative=apply(attr(risks,"cumulative"), 1, weighted.mean, weights=data.ph2$wt))
        
      })
      risks.all.ter[[a]] = list(time=out[[1]][,"t"], risk=cbind(out[[1]][,"cumulative"], out[[2]][,"cumulative"], if(length(out)==3) out[[3]][,"cumulative"]))

    }
    

  } else { # not competing risk
    
    f1=update(form.0, as.formula(paste0("~.+",marker.name)))        
    
    if (TRIAL == c("vat08_combined")) {
      
      # multiple imputation
      out=lapply(1:10, function(imp) {
        dat = set.mi.data(dat, config.cor, imp, marker.name)
        fit.risk=try(svycoxph(f1, design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat)))
        # fit.risk=try(svycoxph(f1, design=svydesign(id=~1, strata=~Wstratum, weights=~wt, data=dat[dat$ph2==1,])))
        if (inherits(fit.risk, "try-error")) {
          NULL
        } else {
          marginalized.risk(fit.risk, marker.name, dat[dat$ph2==1,], weights=dat[dat$ph2==1,"wt"], categorical.s=T, t.end=tfinal.tpeak)
        }
      })
      
      all.t=lapply(out, function(x) x$time)
      common.t = Reduce(intersect, all.t)
      risks = sapply(out, simplify="array", function(x) x$risk[x$time %in% common.t,])
      risks.all.ter[[a]]=list(time=common.t, risk=apply(risks, 1:2, mean, na.rm=T))
      
    } else {
      
      if (all(dat$wt==1)) {
        fit.risk=try(coxph(f1, dat))         # all ph1 are ph2; syvcoxph throws an error, thus use coxph
      } else {
        fit.risk=try(svycoxph(f1, design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat)))
      }
      
      if(inherits(fit.risk,"try-error")) {
        risks.all.ter[[a]]=NA
      } else {
        risks.all.ter[[a]]=marginalized.risk(fit.risk, marker.name, dat[dat$ph2==1,], 
                                             weights=dat[dat$ph2==1,"wt"], categorical.s=T, t.end=tfinal.tpeak)
      }
      
    }
    
  }
} # end: compute risks.all.ter


# compute time.0 and risk.0 from placebo
if(has.plac) {
  
  if (TRIAL=="janssen_partA_VL") {
    out=lapply(1:10, function(imp) {
      dat.plac = set.mi.data(dat.plac, config.cor, imp)
      risks = pcr2(form.0, dat.plac, tfinal.tpeak)
      cbind(t=attr(risks,"time"), cumulative=apply(attr(risks,"cumulative"), 1, mean))
    })
    all.t=lapply(out, function(x) x[,"t"])
    common.t = Reduce(intersect, all.t)
    risks = sapply(out, simplify="array", function(x) x[x[,"t"] %in% common.t,-1])
    time.0=common.t
    risk.0=apply(risks, 1, mean)
    
    
  } else if (TRIAL=="vat08_combined") {
    
    out=lapply(1:10, function(imp) {
      dat.plac = set.mi.data(dat.plac, config.cor, imp)
      
      # CODE: without model=T, there are errors: object EventTimeOfInterest not found
      fit.0=coxph(form.s, dat.plac, model=T) 
      risk.0= 1 - exp(-predict(fit.0, type="expected"))
      time.0= dat.plac[[config.cor$EventTimePrimary %.% imp]]
      # risk.0 for 7 and 7+ are different
      keep=dat.plac[[config.cor$EventIndPrimary %.% imp]]==1 & time.0<=tfinal.tpeak 

      # # use ph2 data to make inference due to issue with prev_inf
      # svycoxph fails because it cannot handle formula with rhs ~1
      # fit.0=coxph(form.s, dat.plac[dat.plac$ph2==1,], model=T, weights=dat.plac[dat.plac$ph2==1,"wt"])
      # time.0= dat.plac[dat.plac$ph2==1,config.cor$EventTimePrimary %.% imp]
      # keep=dat.plac[dat.plac$ph2==1,config.cor$EventTimePrimary %.% imp] & time.0<=tfinal.tpeak
      
      risk.0 = risk.0[keep]
      time.0 = time.0[keep]
      unique_indices <- which(!duplicated(time.0))
      list(time=time.0[unique_indices], risk=risk.0[unique_indices])
    })
    
    all.t=lapply(out, function(x) x$time)
    common.t = Reduce(intersect, all.t)
    risks = sapply(out, simplify="array", function(x) x$risk[x$time %in% common.t])
    time.0=common.t; risk.0=apply(risks, 1, mean, na.rm=T)
    
      
  } else {
    fit.0=coxph(form.s, dat.plac) 
    risk.0= 1 - exp(-predict(fit.0, type="expected"))
    time.0= dat.plac[[config.cor$EventTimePrimary]]
    # risk.0 for 7 and 7+ are different
    keep=dat.plac[[config.cor$EventIndPrimary]]==1 & time.0<=tfinal.tpeak
    risk.0 = risk.0[keep]
    time.0 = time.0[keep]
    
  }
  
}


{
#### set ylim

lwd=2
ylim=c(0,
       max(
         max(sapply(markers, function(a) {
           tmp = risks.all.ter[[a]]
           if (is.null(tmp) || (length(tmp) == 1 && is.na(tmp))) {
             # tmp is NA
             NA
           } else {
             # tmp is a list
             max(tmp$risk[tmp$time<=tfinal.tpeak,])
           }
         }
       ), na.rm=T),
         if(has.plac) risk.0)
       )

# special cases

# same ylim for bAb and nAb
if (COR=="D57azd1222_stage2_delta_nAb" | COR=="D57azd1222_stage2_delta_bAb") {
  ylim=c(0,.09)
} else if (COR=="D57azd1222_stage2_severe_nAb" | COR=="D57azd1222_stage2_severe_bAb") {
  ylim=c(0,.009)
}


if (config$is_ows_trial) {
  x.time<-seq(0,tfinal.tpeak,by=30)
  if(tfinal.tpeak-mylast(x.time)>15) x.time=c(x.time, tfinal.tpeak) else x.time[length(x.time)]=tfinal.tpeak
} else {
  x.time<-floor(seq(0,tfinal.tpeak,length=8))
}
#
if(.mfrow[1]==1)  height=7.5/2*1.5 else height=7.5/2*.mfrow[1]*1.3

assay_units = sapply(assay_metadata$assay_label_short, function(x) {
  out = gsub("\\(|\\)", "", regmatches(x, gregexpr("\\((.*?)\\)", x))[[1]])
  if (length(out)==0) out=""
  out
})
if (length(assay_units)==0) {
  assay_units=rep("",nrow(assay_metadata))
}
names(assay_units)=assay_metadata$assay

}


if (file.exists(paste0(save.results.to, "svycoxph_cat_overall_pvalues_",fname.suffix,".Rdata"))) {
  load(paste0(save.results.to, "svycoxph_cat_overall_pvalues_",fname.suffix,".Rdata"))
}

# make plot for one marker at a time till the end of tertile incidence curves
for (a in markers) {        

  mypdf(oma=c(1,0,0,0), onefile=F, file=paste0(save.results.to, a, "_marginalized_risks_cat_", fname.suffix), mfrow=.mfrow, mar=c(12,4,5,2))
  par(las=1, cex.axis=0.9, cex.lab=1)# axis label 
  
  marker.name=a%.%"cat"    
  if(verbose>=2) myprint(a)
  
  has.3.levels = length(levels(dat[[marker.name]]))==3
  
  out=risks.all.ter[[a]]
  # cutpoints
  q.a=marker.cutpoints[[a]]
  
  # make cumulative incidence plot
  if(length(out)==1) empty.plot() else {
    
    # the first point is (tpeaklag, 0)
    mymatplot(c(tpeaklag, out$time[out$time<=tfinal.tpeak]), 
              rbind(0,out$risk[out$time<=tfinal.tpeak,]), 
              lty=c(1, if(has.3.levels) 5, 2), 
              col=c("darkgreen", if(has.3.levels) "green3", "green"), 
              type="l", lwd=lwd, make.legend=F, ylab=paste0("Probability of ",config.cor$txt.endpoint), ylim=ylim, xlab="", las=1, 
              xlim=c(tpeaklag,tfinal.tpeak), at=x.time, xaxt="n")
    
    title(main=markers.names.long[a], cex.main=1.1, line=1.5)
    title(main=for.title, cex.main=.9, line=.4)
    title(xlab="Days Since Day "%.%tpeak1%.%" Visit", line=2)
    
    # add tpeaklag on the x axis 
    axis(side=1, at=tpeaklag, label=tpeaklag)
    
    # if (has.3.levels) {
    #   mtext(bquote(cutpoints: list(.(formatDouble(10^q.a[1]/10^floor(q.a[1]),1)) %*% 10^ .(floor(q.a[1])), 
    #                                .(formatDouble(10^q.a[2]/10^floor(q.a[2]),1)) %*% 10^ .(floor(q.a[2])))), line= 12.2, cex=.8, side=1)
    # } else { 
    #   mtext(bquote(cutpoints: list(.(formatDouble(10^q.a[1]/10^floor(q.a[1]),1)) %*% 10^ .(floor(q.a[1])))), line= 12.2, cex=.8, side=1)
    # }
    
    # fold change should not have units
    if (startsWith(a,"Delta")) assay_unit="" else assay_unit = assay_units[marker.name.to.assay(a)]
    # add a space if assay_unit is not empty
    if (assay_unit!="") assay_unit = paste0(" ", assay_units[marker.name.to.assay(a)])
    
    if(has.3.levels) {
      legend=c(paste0(trt.label, " low (<",   prettyprint(10^q.a[1]), assay_unit, ")"), 
               paste0(trt.label, " medium (", prettyprint(10^q.a[1]), " to <",
                                          prettyprint(10^q.a[2]), assay_unit,")"),
               paste0(trt.label, " high (>=", prettyprint(10^q.a[2]), assay_unit,")"),
               if(has.plac) cmp.label)
    } else {
      legend=c(paste0(trt.label, " low (<",   prettyprint(10^q.a[1]), assay_unit, ")"), 
               paste0(trt.label, " high (>=", prettyprint(10^q.a[1]), assay_unit, ")"), 
               if(has.plac) cmp.label)
      
    }
    mylegend(x=1, legend=legend, lty=c(1, if(has.3.levels) 5, 2, if(has.plac) 1), 
             col=c("darkgreen", if(has.3.levels) "green3","green",if(has.plac) "gray"), lwd=2, cex=.8)
    if(has.plac) mylines(time.0, risk.0, col="gray", lwd=2, type="l")
    
    # add a legend to show overall p value
    if (exists("overall.p.tri")) mylegend(x=4, legend=paste0("Overall P value: ", formatDouble(overall.p.tri[a],3)), lty=1, col="white", cex=0.8)
    
    
  }
  
  # save source data for images per some journals' requirements
  if (is.null(out) || (length(out) == 1 && is.na(out))) {
    # out is NA
    # write an empty file for caption
    write("", 
          file=paste0(save.results.to, a, "_tertile_incidences_", fname.suffix, ".txt"))
    
  } else {
    img.dat=cbind(out$time[out$time<=tfinal.tpeak], out$risk[out$time<=tfinal.tpeak,])
    rownames(img.dat)=img.dat[,1]
    
    # write a file for caption
    write(paste0(concatList(formatDouble(
      c(img.dat[nrow(img.dat),-1], if(has.plac) max(risk.0, na.rm=T))
      , 6), ", "), "%"), 
          file=paste0(save.results.to, a, "_tertile_incidences_", fname.suffix, ".txt"))
    
    if(has.plac) {
      tmp=cbind(time.0, risk.0)
      tmp=tmp[order (tmp[,1]),]
      tmp=unique(tmp)    
      rownames(tmp)=tmp[,1]
      
      # comment out next line b/c img.dat made from predict using curve type has more rows in tmp
      # next line generates an error about duplicate row names
      # img.dat=cbinduneven(list(img.dat, tmp)) # combine 
    }
    # sort
    img.dat=img.dat[order(img.dat[,1]),]
    mywrite.csv(img.dat, paste0(save.results.to, a, "_marginalized_risks_cat_", fname.suffix))
  
  }
  
  # add data ribbon
  
  if (TRIAL=="janssen_partA_VL") {
    # ancestral markers have different weights and ph2 indicators
    if (a %in% c("Day29bindSpike","Day29pseudoneutid50")) {
      dat$ph2 = dat$ph2.D29
      dat$wt = dat$wt.D29
    } else {
      dat$ph2 = dat$ph2.D29variant
      dat$wt = dat$wt.D29variant
    }
  }
  
  # the use of cbinduneven helps to get around these exceptions if they do occur
  # stopifnot(all(tmp$time[1:length(x.time)]==x.time))
  # stopifnot(tmp$time[1:length(x.time)+length(x.time)]==x.time)
  # stopifnot(tmp$time[1:length(x.time)+length(x.time)*2]==x.time)
  
  # out is an array. if not multiple imputation, the mylast dimension has length 1
  out=sapply(1:ifelse(mi, 10, 1), simplify="array", function(imp) {
    
    if (mi) dat = set.mi.data(dat, config.cor, imp, marker.name)
    
    f1=update(form.s, as.formula(paste0("~.+",marker.name)))
    km <- survfit(f1, dat[dat$ph2==1,], weights=dat[dat$ph2==1,"wt"])
    tmp=summary(km, times=x.time)            
    
    if (has.3.levels) {
      L.idx=which(tmp$time==0)[1]:(which(tmp$time==0)[2]-1)
      n.risk.L <- round(tmp$n.risk[L.idx])
      cum.L <- round(cumsum(tmp$n.event[L.idx]))
      tmp.L = cbind(n.risk.L, cum.L)
      rownames(tmp.L)=tmp$time[L.idx]
      
      M.idx=which(tmp$time==0)[2]:(which(tmp$time==0)[3]-1)
      n.risk.M <- round(tmp$n.risk[M.idx])
      cum.M <- round(cumsum(tmp$n.event[M.idx]))
      tmp.M = cbind(n.risk.M, cum.M)
      rownames(tmp.M)=tmp$time[M.idx]
      
      H.idx=which(tmp$time==0)[3]:length(tmp$time==0)
      n.risk.H <- round(tmp$n.risk[H.idx])
      cum.H <- round(cumsum(tmp$n.event[H.idx]))
      tmp.H = cbind(n.risk.H, cum.H)
      rownames(tmp.H)=tmp$time[H.idx]
      
    } else {
      
      L.idx=which(tmp$time==0)[1]:(which(tmp$time==0)[2]-1)
      n.risk.L <- round(tmp$n.risk[L.idx])
      cum.L <- round(cumsum(tmp$n.event[L.idx]))
      tmp.L = cbind(n.risk.L, cum.L)
      rownames(tmp.L)=tmp$time[L.idx]
      
      H.idx=which(tmp$time==0)[2]:length(tmp$time==0)
      n.risk.H <- round(tmp$n.risk[H.idx])
      cum.H <- round(cumsum(tmp$n.event[H.idx]))
      tmp.H = cbind(n.risk.H, cum.H)
      rownames(tmp.H)=tmp$time[H.idx]
    }
    
    
    # add placebo
    if (has.plac) {
      
      if(mi) dat.plac = set.mi.data(dat.plac, config.cor, imp)
      
      survfit.P=summary(survfit(form.s, dat.plac), times=x.time)            
      n.risk.P <- round(survfit.P$n.risk)
      cum.P <- round(cumsum(survfit.P$n.event))  
      tmp.P = cbind(n.risk.P, cum.P)
      rownames(tmp.P)=survfit.P$time
      if(has.3.levels) {
        data.ribbon = cbinduneven(list(tmp.L, tmp.M, tmp.H, tmp.P))
      } else {
        data.ribbon = cbinduneven(list(tmp.L, tmp.H, tmp.P))
      }
    } else {
      if(has.3.levels) {
        data.ribbon = cbinduneven(list(tmp.L, tmp.M, tmp.H))
      } else {
        data.ribbon = cbinduneven(list(tmp.L, tmp.H))
      }
    }
    
    as.matrix(data.ribbon)
    
  })
    
  if (is.list(out)) stop("cor_coxph_risk_plotting.R: data ribon issue - getting different length output for imputed data")
  data.ribbon = apply(out, 1:2, mean)
  data.ribbon=as.data.frame(data.ribbon)
    
  # data.ribbon=data.ribbon[as.numeric(rownames(data.ribbon))>tpeaklag, ] # truncate x.time to after tpeakkag, thus time 0 is not shown
  
  x.time.1=as.numeric(rownames(data.ribbon))
  # since the population at risk starts at tpeaklag, the counts really start at tpeaklag
  if (x.time.1[1]==0) x.time.1[1] = tpeaklag 
  
  cex.text <- 0.7
  x.label= tpeaklag - (tfinal.tpeak-tpeaklag)/8
  
  line.start=2.2
  
  mtext("No. at risk",side=1,outer=FALSE,line=0.5+line.start,at=tpeaklag-2,adj=0,cex=cex.text)
  mtext(paste0("Low:"),side=1,outer=F,line=1.4+line.start,at=x.label,adj=0,cex=cex.text);  mtext(data.ribbon$n.risk.L,side=1,outer=FALSE,line=1.4+line.start,at=x.time.1,cex=cex.text)
  if (has.3.levels) {
    mtext(paste0("Med:"),side=1,outer=F,line=2.3+line.start,at=x.label,adj=0,cex=cex.text);  mtext(data.ribbon$n.risk.M,side=1,outer=FALSE,line=2.3+line.start,at=x.time.1,cex=cex.text)
  }
  mtext(paste0("High:"),side=1,outer=F,line=3.2+line.start,at=x.label,adj=0,cex=cex.text); mtext(data.ribbon$n.risk.H,side=1,outer=FALSE,line=3.2+line.start,at=x.time.1,cex=cex.text)
  if (has.plac) {
    mtext(paste0(cmp.label, ":"),side=1,outer=F,line=4.2+line.start,at=x.label,adj=0,cex=cex.text); mtext(data.ribbon$n.risk.P,side=1,outer=FALSE,line=4.2+line.start,at=x.time.1,cex=cex.text)
  }
  
  mtext(paste0("Cumulative No. of ",config.cor$txt.endpoint," Endpoints"),side=1,outer=FALSE,line=5.4+line.start,at=tpeaklag-2,adj=0,cex=cex.text)
  mtext(paste0("Low:"),side=1,outer=FALSE,line=6.3+line.start,at=x.label,adj=0,cex=cex.text);  mtext(data.ribbon$cum.L,side=1,outer=FALSE,line=6.3+line.start, at=x.time.1,cex=cex.text)
  if (has.3.levels) {
    mtext(paste0("Med:"),side=1,outer=FALSE,line=7.2+line.start,at=x.label,adj=0,cex=cex.text);  mtext(data.ribbon$cum.M,side=1,outer=FALSE,line=7.2+line.start ,at=x.time.1,cex=cex.text)
  }
  mtext(paste0("High:"),side=1,outer=FALSE,line=8.1+line.start,at=x.label,adj=0,cex=cex.text);mtext(data.ribbon$cum.H,side=1,outer=FALSE,line=8.1+line.start,at=x.time.1,cex=cex.text)
  if (has.plac) {
    mtext(paste0(cmp.label, ":"),side=1,outer=FALSE,line=9.1+line.start,at=x.label,adj=0,cex=cex.text);mtext(data.ribbon$cum.P,side=1,outer=FALSE,line=9.1+line.start,at=x.time.1,cex=cex.text)
  }
  
  dev.off()    
}
#mtext(toTitleCase(study_name), side = 1, line = 2, outer = T, at = NA, adj = NA, padj = NA, cex = .8, col = NA, font = NA)
#
#cumsum(summary(survfit(form.s, subset(dat, ph2==1)), times=x.time)$n.event)
#table(subset(dat, yy==1)[["Day"%.%tpeak%.%"pseudoneutid80cat"]])





###################################################################################################
cat("plot trichotomized markers, log(-log) marginalized survival curves\n")
# for goodness of fit check on PH assumptions

for (a in markers) {        
  mypdf(onefile=F, file=paste0(save.results.to, a, "_marginalized_risks_cat_logclog_",fname.suffix), mfrow=.mfrow)
  marker.name=a%.%"cat"    
  
  out=risks.all.ter[[a]]
  # cutpoints
  q.a=marker.cutpoints[[a]]
  
  if(length(out)==1) empty.plot() else {
    mymatplot(out$time[out$time<=tfinal.tpeak], log(-log(out$risk[out$time<=tfinal.tpeak,])), 
              lty=1:3, col=c("darkgreen","green3","green"), type="l", lwd=lwd, make.legend=F, 
              ylab=paste0("log(-log( Probability of ",config.cor$txt.endpoint," by Day "%.%tfinal.tpeak, " ))"), xlab="", 
              las=1, xlim=c(0,tfinal.tpeak), at=x.time, xaxt="n")
    title(xlab="Days Since Day "%.%tpeak1%.%" Visit", line=2)
    title(main=markers.names.long[a], cex.main=.9, line=2)
    title(main=for.title, line=.4, cex.main=.9)
    mtext(bquote(cutpoints: list(.(formatDouble(10^q.a[1]/10^floor(q.a[1]),1)) %*% 10^ .(floor(q.a[1])), 
                                 .(formatDouble(10^q.a[2]/10^floor(q.a[2]),1)) %*% 10^ .(floor(q.a[2])))), line= 2.8, cex=.8, side=1)   
    legend=trt.label%.%c(" low"," medium"," high")
    mylegend(x=3, legend=legend, lty=c(1:3), col=c("green3","green","darkgreen"), lwd=2)
  }
  
  dev.off()    
}



}
