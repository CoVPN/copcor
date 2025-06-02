# _2arms means putting 3 curves for both arms together in the same plot
cor_coxph_risk_tertile_incidence_curves_2arms = function(
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
  verbose=FALSE,
  
  fname.suffix.2
) {
  

# load overall.p.tri from both arms
overall.p.tri.ls = list()
tmp=paste0(save.results.to, paste0("coxph_overall_p_tri_", fname.suffix, ".Rdata"))
if (file.exists(tmp)) load(tmp)  
overall.p.tri.ls[[1]] = overall.p.tri
# load overall.p.tri
tmp=paste0(save.results.to, paste0("coxph_overall_p_tri_", fname.suffix.2, ".Rdata"))
if (file.exists(tmp)) load(tmp)  
overall.p.tri.ls[[2]] = overall.p.tri

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

# form.s is the null model, i.e. rhs is ~1
form.s=as.formula(deparse((if(comp.risk) form.0[[1]] else form.0)[[2]])%.%"~1")


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

compute.risks.all.ter = function (dat) {
  data.ph2 <- dat[dat$ph2==1,]
  risks.all.ter=list()
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
  }
  risks.all.ter
}

risks.all.ter.arm1=compute.risks.all.ter(dat)
risks.all.ter.arm2=compute.risks.all.ter(dat.plac)



{
#### set ylim

lwd=2
ylim=c(0, max(
         sapply(markers, function(a) max(risks.all.ter.arm1[[a]]$risk [risks.all.ter.arm1[[a]]$time<=tfinal.tpeak,])),
         sapply(markers, function(a) max(risks.all.ter.arm2[[a]]$risk [risks.all.ter.arm2[[a]]$time<=tfinal.tpeak,]))
       ))
       

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


# make plot for one marker at a time till the end of tertile incidence curves
for (a in markers) {        
  mypdf(oma=c(1,0,0,0), onefile=F, mfrow=.mfrow, mar=c(3,4,3,2), width=6.7, height=4.7, file=paste0(save.results.to, a, "_marginalized_risks_cat_", fname.suffix, "_", fname.suffix.2))
  par(las=1, cex.axis=0.9, cex.lab=1)# axis label 
  
  marker.name=a%.%"cat"    
  if(verbose>=2) myprint(a)
  
  has.3.levels = length(levels(dat[[marker.name]]))==3
  
  col1 = c("darkgreen", if(has.3.levels) "green3", "green")
  col2 = c("gray40", if(has.3.levels) "gray60", "gray80")
  
  # use the same colors for low med high
  if (startsWith(TRIAL, "nextgen")) {
    col1 = c("darkgreen", if(has.3.levels) "darkgreen", "darkgreen")
    col2 = c("gray40", if(has.3.levels) "gray40", "gray40")
  }
  
  for (i in 1:2) { # 1: arm1; 2: arm2
    add = i==2
    
    out=get("risks.all.ter.arm"%.%i)[[a]]
    
    # make cumulative incidence plot
    if(length(out)==1) {
      empty.plot() 
      break
      
    } else {
      # the first point is (tpeaklag, 0)
      mymatplot(c(tpeaklag, out$time[out$time<=tfinal.tpeak]), 
                rbind(0,out$risk[out$time<=tfinal.tpeak,]), 
                lty=c(1, if(has.3.levels) 5, 2), 
                col=get("col"%.%i), 
                type="l", lwd=lwd, make.legend=F, ylab=paste0("Probability of ",config.cor$txt.endpoint), ylim=ylim, xlab="", las=1, 
                xlim=c(tpeaklag,tfinal.tpeak), at=x.time, xaxt="n", add=add)
      
    }
    
  } # end for i in 1:2

  q.a=marker.cutpoints[[a]]
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
    legend=c(paste0(trt.label, " low  (<",   prettyprint(10^q.a[1]), assay_unit, ")"), 
             paste0(trt.label, " med (", prettyprint(10^q.a[1]), " to <", prettyprint(10^q.a[2]), assay_unit,")"),
             paste0(trt.label, " high (>=", prettyprint(10^q.a[2]), assay_unit,")"),
             paste0(cmp.label, " low  (<",   prettyprint(10^q.a[1]), assay_unit, ")"), 
             paste0(cmp.label, " med (", prettyprint(10^q.a[1]), " to <", prettyprint(10^q.a[2]), assay_unit,")"),
             paste0(cmp.label, " high (>=", prettyprint(10^q.a[2]), assay_unit,")")
    )
  } else {
    legend=c(paste0(trt.label, " low  (<",   prettyprint(10^q.a[1]), assay_unit, ")"), 
             paste0(trt.label, " high (>=", prettyprint(10^q.a[1]), assay_unit, ")"),
             paste0(cmp.label, " low  (<",   prettyprint(10^q.a[1]), assay_unit, ")"), 
             paste0(cmp.label, " high (>=", prettyprint(10^q.a[1]), assay_unit, ")")
    )
    
  }
  mylegend(x=1, legend=legend, lty=c(1, if(has.3.levels) 5, 2), 
           col=c(col1, col2), lwd=2, cex=.8)
  
  # # add a legend to show overall p value
  # if (exists("overall.p.tri")) mylegend(x=4, legend=paste0("Overall P value: ", formatDouble(overall.p.tri[a],3)), lty=1, col="white", cex=0.8)'
  
  dev.off()    
}
}