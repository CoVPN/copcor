# competing risk is handled transparently through form.0, e.g.
# multiple imputation is hardcoded by TRIAL

# janssen_partA_VL handling is hard coded. 
#   outcome is multiple imputed; depends on variant, which is defined globally
#   marker is multiple imputed
#   ph2 and wt depends on marker

cor_coxph_risk_plotting = function(
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
  res.plac.cont = NULL,
  prev.plac=NULL,
  overall.ve=NULL,
  
  show.ve.curves=T,
  plot.geq = F, # whether to plot risk vs S>=s
  plot.w.plac = T, # whether to plot a version with plac line
  for.title="",
  
  verbose=FALSE
) {
  

  
if(verbose) print("Running cor_coxph_risk_plotting")
  
comp.risk = is.list(form.0)
  
if(is.null(tfinal.tpeak)) tfinal.tpeak = config.cor$tfinal.tpeak
if (is.null(tfinal.tpeak)) stop("missing tfinal.tpeak")
  
has.plac=!is.null(dat.plac)
  
eq.geq.ub=ifelse(plot.geq, 2, 1)
wo.w.plac.ub=ifelse(plot.w.plac, 2, 1)

myprint(comp.risk, has.plac, plot.geq, plot.w.plac, for.title)
  
# make form.s from form.0
form.s=as.formula(deparse((if(comp.risk) form.0[[1]] else form.0)[[2]])%.%"~1")

# compute prevalence
# note that these do not have CI
prev.vacc = if (TRIAL %in% c("janssen_partA_VL")) {
  mean(sapply(1:10, function(imp) {
    dat$EventIndOfInterest = ifelse(dat$EventIndPrimary==1 & dat[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
    dat$EventIndCompeting  = ifelse(dat$EventIndPrimary==1 & dat[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
    get.marginalized.risk.no.marker(form.0, dat, tfinal.tpeak)
  }))
  
} else  if (TRIAL == c("vat08_combined")) {
  mean(sapply(1:10, function(imp) {
    dat$EventIndOfInterest  = dat[[config.cor$EventIndPrimary  %.% imp]]
    dat$EventTimeOfInterest = dat[[config.cor$EventTimePrimary %.% imp]]
    get.marginalized.risk.no.marker(form.0, dat, tfinal.tpeak)
  }))
  
} else {
  get.marginalized.risk.no.marker(form.0, dat, tfinal.tpeak)
}

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


###################################################################################################
cat("make trichotomized markers, marginalized risk and controlled risk table\n")

res=sapply (markers, function(a) {        
  risks=get("risks.all.3")[[a]]
  has.3.levels = length(risks$prob)==3
  if (has.3.levels) {
    with(risks, c(prob[3]/prob[1], quantile(boot[3,]/boot[1,], c(.025,.975), na.rm=T)))
  } else {
    with(risks, c(prob[2]/prob[1], quantile(boot[2,]/boot[1,], c(.025,.975), na.rm=T)))
  }
})
#    
tab=sapply (markers, function(a) {
  paste0(
    markers.names.short[a], "&",
    # marginal RR and ci
    formatDouble(res[1,a],2,remove.leading0=F), "&", formatDouble(res[2,a],2,remove.leading0=F), "--", ifelse(res[3,a]>1000,"Inf",formatDouble(res[3,a],2,remove.leading0=F))
    , "&" ,
    # causal RR and ci
    formatDouble(res[1,a]*bias.factor,2,remove.leading0=F), "&", formatDouble(res[2,a]*bias.factor,2,remove.leading0=F), "--", ifelse(res[3,a]*bias.factor>1000,"Inf",formatDouble(res[3,a]*bias.factor,2,remove.leading0=F))
    , "&" ,
    # E-value and ub
    formatDouble(E.value(res[1,a]),1), "&", formatDouble(E.value(res[3,a]),1)
  )
})
write(concatList(tab, "\\\\"), file=paste0(save.results.to, "marginalized_risks_cat_", fname.suffix,".tex"))




###################################################################################################
cat("plot marginalized risk curves for continuous markers\n")

for (eq.geq in 1:eq.geq.ub) {  # 1 conditional on s,   2 is conditional on S>=s
  for (wo.w.plac in 1:wo.w.plac.ub) { # 1 with placebo lines, 2 without placebo lines. Implementation-wise, the main difference is in ylim
    # eq.geq=1; wo.w.plac=2; a=markers[1]
    
    risks.all=get("risks.all."%.%eq.geq)
    
    if (!create.ylims.cor) {
      ylim=ylims.cor[[eq.geq]][[wo.w.plac]]
      
    } else {
      if(verbose>=2) print("no ylims.cor found")        
      if (eq.geq==2 & wo.w.plac==1) {
        # later values in prob may be wildly large due to lack of samples
        ylim=range(sapply(risks.all, function(x) x$prob[1]), if(wo.w.plac==2) prev.plac, prev.vacc, 0)
        # may need to add some white space at the top to write placebo overall risk
      } else {
        ylim = range(lapply(markers, function(a) {
          risks=risks.all[[a]]
          # in some datasets, not all markers have risks
          if (!is.null(risks)) {
            shown=risks$marker>=wtd.quantile(dat[[a]], dat$wt, 2.5/100) & risks$marker<=wtd.quantile(dat[[a]], dat$wt, 1-2.5/100)
            risks$prob[shown]
          } else {
            NULL
          }
        }))
        ylim=range(ylim, prev.vacc, 0, if(wo.w.plac==2) prev.plac)
      }
      ylims.cor[[eq.geq]][[wo.w.plac]]=ylim
    }
    # the following is commented out because it depends on COR
    # make the ylim look comparable to ID50 in this trial
    # if (attr(config,"config")=="azd1222_bAb" & eq.geq==1 & wo.w.plac==2 & COR=="D57") ylim=c(0,0.05)    
    if(verbose>=2) myprint(ylim)
    lwd=2
    
    for (a in markers) {        
      mypdf(oma=c(0,0,0,0), onefile=F, file=paste0(save.results.to, a, "_marginalized_risks", ifelse(eq.geq==1,"_eq","_geq"), ifelse(wo.w.plac==2,"","_woplacebo"), "_"%.%fname.suffix), mfrow=.mfrow)
      
      par(las=1, cex.axis=0.9, cex.lab=1)# axis label orientation
      risks=risks.all[[a]]
      assay=marker.name.to.assay(a)
      is.delta=startsWith(a,"Delta")
      
      ncases=sapply(risks$marker, function(s) sum(dat$yy[dat[[a]]>=s], na.rm=T))
      
      if (!is.delta) xlim=get.range.cor(dat, assay, tpeak) else xlim=range(dat[[a]], na.rm=T)
      shown=risks$marker>=ifelse(study_name=="COVE",log10(10),wtd.quantile(dat[[a]], dat$wt, 2.5/100)) & 
        risks$marker<=wtd.quantile(dat[[a]], dat$wt, 1-2.5/100)
      plot(risks$marker[shown], risks$prob[shown], 
           xlab=markers.names.short[a]%.%ifelse(eq.geq==1," (=s)"," (>=s)"), 
           xlim=xlim, ylab=paste0("Probability* of ",config.cor$txt.endpoint," by ", tfinal.tpeak, " days post Day ", tpeak1, " Visit"), lwd=lwd, ylim=ylim, 
           type="n", main=paste0(markers.names.long[a]), xaxt="n")
      title(main=for.title, line=.6, cex.main=.9)
      draw.x.axis.cor(xlim, llox=lloxs[assay], if(is.delta) "delta" else llox_labels[assay])
      
      # prevalence lines
      if (has.plac & wo.w.plac==2) abline(h=prev.plac[1], col="gray", lty=c(2,3,3), lwd=lwd)
      
      # risks
      if (eq.geq==1) {
        abline(h=prev.vacc, col="darkgray", lty=c(1,3,3), lwd=lwd)
        lines(risks$marker[shown], risks$prob[shown], lwd=lwd)
        lines(risks$marker[shown], risks$lb[shown],   lwd=lwd, lty=3)
        lines(risks$marker[shown], risks$ub[shown],   lwd=lwd, lty=3)    
        img.dat=cbind(risks$marker[shown], risks$prob[shown], risks$lb[shown], risks$ub[shown])
      } else {
        abline(h=prev.vacc[1], col="gray", lty=c(1), lwd=lwd)
        lines(risks$marker[ncases>=5], risks$prob[ncases>=5], lwd=lwd)
        lines(risks$marker[ncases>=5], risks$lb[ncases>=5],   lwd=lwd, lty=3)
        lines(risks$marker[ncases>=5], risks$ub[ncases>=5],   lwd=lwd, lty=3)    
        img.dat=cbind(risks$marker[ncases>=5], risks$prob[ncases>=5], risks$lb[ncases>=5], risks$ub[ncases>=5])
      }
      
      # save to satisfy some journal requirements
      mywrite.csv(img.dat, file=paste0(save.results.to, a, "_risk_curves",ifelse(eq.geq==1,"_eq","_geq"),"_"%.%fname.suffix))    
      
      # text overall risks
      if (wo.w.plac==2) {
        text(x=par("usr")[2]-diff(par("usr")[1:2])/3.5, y=prev.plac[1]+(prev.plac[1]-prev.plac[2])/2, "placebo overall "%.%formatDouble(prev.plac[1],3,remove.leading0=F))        
        text(x=par("usr")[2]-diff(par("usr")[1:2])/3.5, y=prev.vacc[1]+(prev.plac[1]-prev.plac[2])/2, "vaccine overall "%.%formatDouble(prev.vacc[1],3,remove.leading0=F))
      } else {
        if (has.plac) text(x=par("usr")[2]-diff(par("usr")[1:2])/3.5, y=par("usr")[4]-diff(par("usr")[3:4])/20,     "placebo overall "%.%formatDouble(prev.plac[1],3,remove.leading0=F))
        text(x=par("usr")[2]-diff(par("usr")[1:2])/4.5, y=prev.vacc[1]+diff(par("usr")[3:4])/20, "overall "%.%formatDouble(prev.vacc[1],3,remove.leading0=F))
      }
      
      # add histogram
      par(new=TRUE) 
      col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
      col <- rgb(col[1], col[2], col[3], alpha=255*0.4, maxColorValue=255)
      tmp.x=dat[[a]][dat$ph2==1]
      tmp.w=dat$wt[dat$ph2==1]
      tmp=get.marker.histogram(tmp.x, tmp.w, attr(config,"config"))
      if (any(is.nan(tmp$density))) tmp=hist(tmp.x, plot=F)
      # plot
      plot(tmp,col=col,axes=F,labels=F,main="",xlab="",ylab="",border=0,freq=F, xlim=xlim, ylim=c(0,max(tmp$density*1.25)))
      #axis(side=4, at=axTicks(side=4)[1:5])
      #mtext("Density", side=4, las=0, line=2, cex=1, at=.3)  
      #mylegend(x=6, fill=col, border=col, legend="Vaccine Group", bty="n", cex=0.7)      
      #mtext(toTitleCase(study_name), side = 1, line = 0, outer = T, at = NA, adj = NA, padj = NA, cex = NA, col = NA, font = NA)
      
      dev.off()    
    } # end assays
    
  }
}
save(ylims.cor, file=paste0(save.results.to, "ylims.cor_"%.%fname.suffix%.%".Rdata"))


# show the results at select assay values
risks.all=get("risks.all.1") 
for (a in markers) {
  risks=risks.all[[a]]
  table.order=which(names(risks$marker) %in% c(" 2.5%", " 5.0%", "10.0%", "50.0%", "90.0%", "95.0%", "97.5%")); table.order=c(setdiff(1:length(risks$marker), table.order), table.order)
  tmp=10**risks$marker[table.order]
  tmp=ifelse(tmp<100, signif(tmp,3), round(tmp))
  out=with(risks, cbind("s"=tmp, "Estimate"=paste0(formatDouble(prob[table.order],digits.risk), " (", formatDouble(lb[table.order],digits.risk), ",", formatDouble(ub[table.order],digits.risk), ")")))
  while (nrow(out)%%4!=0) out=rbind(out, c("s"="", "Estimate"=""))
  tab=cbind(out[1:(nrow(out)/4), ], out[1:(nrow(out)/4)+(nrow(out)/4), ], out[1:(nrow(out)/4)+(nrow(out)/4*2), ], out[1:(nrow(out)/4)+(nrow(out)/4*3), ])
  mytex(tab, file.name=paste0(a, "_marginalized_risks_eq", "_"%.%fname.suffix), align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to, include.rownames = F,
        longtable=T, caption.placement = "top", label=paste0("tab marginalized_risks_eq ", fname.suffix), caption=paste0("Marginalized cumulative risk by Day ",tfinal.tpeak," as functions of Day ",
                                                                                                                         tpeak, " ", markers.names.short[a], " (=s) among baseline negative vaccine recipients with 95\\% bootstrap point-wise confidence intervals (",
                                                                                                                         ncol(risks.all[[1]]$boot)," replicates). ",
                                                                                                                         "Last seven values correspond to 2.5\\%, 5.0\\%, 10.0\\%, 50.0\\%, 90.0\\%, 95.0\\%, 97.5\\%, respectively.")
        #, col.headers=paste0("\\hline\n", concatList(paste0("\\multicolumn{2}{c}{", labels.axis[1,], "}"), "&"), "\\\\\n")
  )
}




###################################################################################################
if(show.ve.curves) {

cat("continuous markers, controlled VE curves\n")

# 1 conditional on s
# 2 conditional on S>=s
# 3 same as 1 except that no sens curve is shown
# 4 same as 3 except that y axis on -log(1-) scale
# 5 same as 1 except that y axis on -log(1-) scale
curve.set=1:5
if(!plot.geq) curve.set=setdiff(curve.set, 2)
for (eq.geq in curve.set) {  
  # eq.geq=4; a=markers[1]
  
  names(markers)=markers # required so that outs has names
  outs=lapply (markers, function(a) {        
    if (verbose>=2) myprint(a)
    is.delta=startsWith(a,"Delta")
    assay=marker.name.to.assay(a)
    
    tmp.1=ifelse(eq.geq==1,"_eq",ifelse(eq.geq==2,"_geq","_eq_manus")); if(eq.geq %in% c(4,5)) tmp.1=eq.geq; 
    
    mypdf(onefile=F, file=paste0(save.results.to, a, "_controlled_ve_curves",tmp.1,"_"%.%fname.suffix), mfrow=.mfrow, oma=c(0,0,0,0))
    
    lwd=2.5
    par(las=1, cex.axis=0.9, cex.lab=1)# axis label orientation
    
    # load risks
    risks=get("risks.all."%.%ifelse(eq.geq==2,2,1))[[a]]   
    table.order=which(names(risks$marker) %in% c(" 2.5%", " 5.0%", "50.0%", "90.0%", "95.0%", "97.5%")); 
    table.order=c(setdiff(1:length(risks$marker), table.order), table.order)
    
    #xlim=quantile(dat[["Day"%.%tpeak%.%a]],if(eq.geq==1) c(.025,.975) else c(0,.95),na.rm=T)
    if (!is.delta) xlim=get.range.cor(dat, assay, tpeak) else xlim=range(dat[[a]], na.rm=T)            
    
    # compute Bias as a vector, which is a function of s
    # choose a reference marker value based on matching the overall risk
    which=which.min(abs(risks$prob-prev.vacc[1]))
    s.ref=risks$marker[which]
    Bias=controlled.risk.bias.factor(ss=risks$marker, s.cent=s.ref, s1=risks$marker[s1], s2=risks$marker[s2], RRud) 
    if (is.nan(Bias[1])) Bias=rep(1,length(Bias))
    
    if (eq.geq==2) {
      if (study_name %in% c("COVE", "MockCOVE")) {
        ylim=c(0.8,1)
      } else if (study_name %in% c("ENSEMBLE", "MockENSEMBLE")) {
        ylim=c(0.5,1)
      } else if (study_name %in% c("HVTN705")) {
        ylim=c(-1,1)
      } else {
        ylim=c(0,1)
      }
    } else if (eq.geq %in% c(4,5)) {
      ylim=-log(1-config$ve_ylim_log)
    } else {
      ylim=config$ve_ylim
    }
    
    ncases=sapply(risks$marker, function(s) sum(dat$yy[dat[[a]]>=s], na.rm=T))        
    .subset=if(eq.geq!=2) {
      risks$marker>=ifelse(study_name=="COVE",log10(10),wtd.quantile(dat[[a]], dat$wt, 2.5/100)) & 
        risks$marker<=wtd.quantile(dat[[a]], dat$wt, 1-2.5/100)
    } else ncases>=5
    
    
    # CVE with sensitivity analysis
    est = 1 - risks$prob*Bias/res.plac.cont["est"]
    boot = 1 - t( t(risks$boot*Bias)/res.plac.cont[2:(1+ncol(risks$boot))] ) # res.plac.cont may have more bootstrap replicates than risks$boot
    ci.band=apply(boot, 1, function (x) quantile(x, c(.025,.975), na.rm=T))
    # for table
    tmp=10**risks$marker[table.order];     tmp=ifelse(tmp<100, signif(tmp,3), round(tmp))
    ret=cbind("s"=tmp, "Estimate"=paste0(formatDouble(est[table.order],digits.risk), " (", formatDouble(ci.band[1,table.order],digits.risk), ",", formatDouble(ci.band[2,table.order],digits.risk), ")"))
    
    # draw CVE curve with sensitivity analysis
    y= t(rbind(est, ci.band))[.subset,]
    if(eq.geq%in% c(4,5)) y=-log(1-y)
    mymatplot(risks$marker[.subset], y, type="l", lty=c(1,2,2), 
              col=ifelse(eq.geq%in% c(1,5),"red","white"), # white is no plot
              lwd=lwd, make.legend=F, 
              ylab=paste0("Controlled VE against ",config.cor$txt.endpoint," by ", tfinal.tpeak, " days post Day ", tpeak1, " Visit"), 
              main=paste0(markers.names.long[a]),
              xlab=markers.names.short[a]%.%ifelse(eq.geq!=2," (=s)"," (>=s)"), 
              ylim=ylim, xlim=xlim, yaxt="n", xaxt="n", draw.x.axis=F)
    # y axis labels
    if (!eq.geq%in% c(4,5)) {
      yat=seq(-1,1,by=.1)
      axis(side=2,at=yat,labels=(yat*100)%.%"%")            
    } else {
      yat=c(seq(-2,0,by=.5),seq(0,.90,by=.1),.95)
      axis(side=2,at=-log(1-yat),labels=(yat*100)%.%"%")            
    }
    # x axis
    draw.x.axis.cor(xlim, lloxs[assay], if(is.delta) "delta" else llox_labels[assay])
    title(main=for.title, line=.6, cex.main=.9)
    
    img.dat=cbind(risks$marker[.subset], t(rbind(est, ci.band))[.subset,])
    
    # add overall CVE horizontal line
    abline(h=if(eq.geq%in% c(4,5)) -log(1-overall.ve) else overall.ve, col="gray", lwd=2, lty=c(1,3,3))
    #text(x=par("usr")[1], y=overall.ve[1]+(overall.ve[1]-overall.ve[2])/2,     "overall VE "%.%round(overall.ve[1]*100)%.%"%", adj=0)
    
    
    # add CVE curve
    est = 1 - risks$prob/res.plac.cont["est"]; boot = 1 - t( t(risks$boot)/res.plac.cont[2:(1+ncol(risks$boot))] )
    #est = 1 - (risks$prob+0.00227)/res.plac.cont["est"]; boot = 1 - t( t(risks$boot+0.00227)/res.plac.cont[2:(1+ncol(risks$boot))] )
    ci.band=apply(boot, 1, function (x) quantile(x, c(.025,.975), na.rm=T))  
    y= t(rbind(est, ci.band))[.subset,]
    if(eq.geq%in% c(4,5)) y=-log(1-y)
    mymatplot(risks$marker[.subset], y, type="l", lty=c(1,2,2), col=if(!eq.geq%in% c(1,5)) "black" else "pink", lwd=lwd, make.legend=F, add=T)
    #            if (config$is_ows_trial) {
    #                # find marker values under specific VE
    #                tmpind=sapply(report.ve.levels, function (x) ifelse (x>min(est)-0.01 & x<max(est)+0.01, which.min(abs(est-x)), NA))
    #                tmp=10**risks$marker[tmpind]; tmp=c(round(tmp[1],1), round(tmp[-1]))
    #                ret=rbind(ret, cbind("s"=tmp, "Estimate"=paste0(formatDouble(est[tmpind],digits.risk), " (", formatDouble(ci.band[1,tmpind],digits.risk), ",", formatDouble(ci.band[2,tmpind],digits.risk), ")")))            
    #            }
    
    img.dat=cbind(img.dat, y)
    mywrite.csv(img.dat, file=paste0(save.results.to, a, "_controlled_ve_curves",tmp.1,"_"%.%fname.suffix))    
    
    
    # legend
    tmp=formatDouble(overall.ve*100,1)%.%"%"        
    legend.x=9; if(eq.geq %in% c(1,3) & config$low_efficacy) legend.x=1; if(eq.geq%in% c(4,5)) legend.x=1
    mylegend(x=legend.x,legend=c(
      paste0("Overall VE ",tmp[1]," (",tmp[2],", ",tmp[3],")"), 
      "Controlled VE",
      if(eq.geq%in% c(1,5)) "Controlled VE Sens. Analysis"), 
      col=c("gray", if(eq.geq%in% c(1,5)) "pink" else "black", if(eq.geq%in% c(1,5)) "red"), 
      lty=1, lwd=2, cex=.8)
    
    
    #            # add segments if needed
    #            newx=log10(c(54,247,563))
    #            fit.tmp=try(svycoxph(update(form.0, as.formula(paste0("~.+",a))), design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat)))        
    #            out=marginalized.risk(fit.tmp, a, dat.ph2, t=tfinal.tpeak, ss=newx, weights=dat.ph2$wt, categorical.s=F)
    #            out=-log(out/res.plac.cont["est"])
    #            segments(newx, -1, newx, out, col=c("darkgreen","darkorchid3","deepskyblue3"), lwd=2)
    #            segments(rep(-2,3),out, newx,out, col=c("darkgreen","darkorchid3","deepskyblue3"), lwd=2)
    
    
    # add histogram
    par(new=TRUE) 
    col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
    col <- rgb(col[1], col[2], col[3], alpha=255*0.4, maxColorValue=255)            
    tmp.x=dat[[a]][dat$ph2==1]
    tmp.w=dat$wt[dat$ph2==1]
    tmp=get.marker.histogram(tmp.x, tmp.w, attr(config,"config"))
    if (is.nan(tmp$density[1])) tmp=hist(tmp.x, plot=F)
    if(eq.geq%in% c(4,5)) tmp$density=tmp$density*3
    plot(tmp,col=col,axes=F,labels=F,main="",xlab="",ylab="",border=0,freq=F,xlim=xlim, ylim=c(0,max(tmp$density*1.25))) 
    
    
    dev.off()    
    
    ret        
  })
  
  
  if(eq.geq==1) {
    # show the results at select assay values
    for (a in markers) { 
      out=outs[[a]]
      while (nrow(out)%%4!=0) out=rbind(out, c("s"="", "Estimate"=""))
      tab=cbind(out[1:(nrow(out)/4), ], out[1:(nrow(out)/4)+(nrow(out)/4), ], out[1:(nrow(out)/4)+(nrow(out)/4*2), ], out[1:(nrow(out)/4)+(nrow(out)/4*3), ])        
      mytex(tab, file.name=paste0(a, "_controlled_ve_sens_eq", "_"%.%fname.suffix), align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to, include.rownames = F,
            longtable=T, caption.placement = "top", label=paste0("tab controlled_ve_sens_eq ", fname.suffix), caption=paste0("Controlled VE with sensitivity analysis as functions of Day ",
                                                                                                                    tpeak," ", markers.names.short[a], " (=s) among baseline negative vaccine recipients with 95\\% bootstrap point-wise confidence intervals (",
                                                                                                                    ncol(risks.all[[1]]$boot)," replicates). ",
                                                                                                                    "Last six values correspond to 2.5\\%, 5.0\\%, 50.0\\%, 90.0\\%, 95.0\\%, 97.5\\%, respectively.")
            #, col.headers=paste0("\\hline\n", concatList(paste0("\\multicolumn{2}{c}{", labels.axis[1,], "}"), "&"), "\\\\\n")
      )
    }
  }
  
} # end for eq.geq


# show tables of controlled ve without sensitivity at select assay values
digits.risk=4
risks.all=get("risks.all.1")
for(a in markers) {        
  risks=risks.all[[a]]
  table.order=which(names(risks$marker) %in% c(" 2.5%", " 5.0%", "50.0%", "90.0%", "95.0%", "97.5%"))
  table.order=c(setdiff(1:length(risks$marker), table.order), table.order)
  
  est = 1 - risks$prob/res.plac.cont["est"]
  boot = 1 - t( t(risks$boot)/res.plac.cont[2:(1+ncol(risks$boot))] )                         
  ci.band=apply(boot, 1, function (x) quantile(x, c(.025,.975), na.rm=T))        
  
  tmp=10**risks$marker[table.order];     tmp=ifelse(tmp<100, signif(tmp,3), round(tmp))
  ret = cbind("s"=tmp, "Estimate"=paste0(formatDouble(est[table.order],digits.risk), " (", formatDouble(ci.band[1,table.order],digits.risk), ",", formatDouble(ci.band[2,table.order],digits.risk), ")"))
  
  #    if (config$is_ows_trial) {
  #        # find marker values under specific VE
  #        tmpind=sapply(report.ve.levels, function (x) ifelse (x>min(est)-0.01 & x<max(est)+0.01, which.min(abs(est-x)), NA))
  #        tmp=10**risks$marker[tmpind]; tmp=c(round(tmp[1],1), round(tmp[-1]))
  #        out=rbind(ret, cbind("s"=tmp, "Estimate"=paste0(formatDouble(est[tmpind],digits.risk), " (", formatDouble(ci.band[1,tmpind],digits.risk), ",", formatDouble(ci.band[2,tmpind],digits.risk), ")")))
  #    } 
  
  while (nrow(ret)%%4!=0) ret=rbind(ret, c("s"="", "Estimate"=""))
  tab=cbind(ret[1:(nrow(ret)/4), ], ret[1:(nrow(ret)/4)+(nrow(ret)/4), ], ret[1:(nrow(ret)/4)+(nrow(ret)/4*2), ], ret[1:(nrow(ret)/4)+(nrow(ret)/4*3), ])
  
  mytex(tab, file.name=paste0(a, "_controlled_ve_eq", "_"%.%fname.suffix), align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to, include.rownames = F,
        longtable=T, caption.placement = "top", label=paste0("tab controlled_ve_eq ", fname.suffix), 
        caption=paste0("Controlled VE as functions of Day ",
                       tpeak," ", markers.names.short[a], " (=s) among baseline negative vaccine recipients with 95\\% bootstrap point-wise confidence intervals (",
                       ncol(risks.all[[1]]$boot)," replicates). ", 
                       if(has.plac) paste0("Overall cumulative incidence from ", tpeaklag, " to ",tfinal.tpeak," days post Day ",tpeak1," was ",
                                           formatDouble(prev.vacc[1], 3, remove.leading0=F)," in vaccine recipients compared to ",
                                           formatDouble(prev.plac[1], 3, remove.leading0=F)," in placebo recipients, with cumulative vaccine efficacy ",
                                           formatDouble(overall.ve[1]*100,1),"\\% (95\\% CI ",formatDouble(overall.ve[2]*100,1)," to ",formatDouble(overall.ve[3]*100,1),"\\%). "),
                       "Last six values correspond to 2.5\\%, 5.0\\%, 50.0\\%, 90.0\\%, 95.0\\%, 97.5\\%, respectively.")
        #, col.headers=paste0("\\hline\n", concatList(paste0("\\multicolumn{2}{c}{", labels.axis[1,], "}"), "&"), "\\\\\n")
  )
  
}


}




###################################################################################################
# interaction models and risk curves

if (!is.null(config$interaction)) {
  risks.itxn = get("risks.itxn")
  
  for (ab in config$interaction) {
    tmp = trim(strsplit(ab, " *\\* *")[[1]])
    aold = tmp[1]
    bold = tmp[2]
    a = paste0("Day", tpeak, aold)
    b = paste0("Day", tpeak, bold)
    
    for (inner.id in 1:2) {
      vx = ifelse(inner.id == 1, a, b)
      vthree = ifelse(inner.id == 1, b, a)
      risks = risks.itxn[[paste0(vx, vthree)]]
      
      mypdf(
        oma = c(0, 0, 0, 0),
        onefile = F,
        file = paste0(
          save.results.to,
          "itxn_marginalized_risks_",
          ifelse(inner.id == 1, aold, bold),
          "_",
          ifelse(inner.id == 1, bold, aold)
        ),
        mfrow = .mfrow
      )
      
      par(las = 1,
          cex.axis = 0.9,
          cex.lab = 1)# axis label orientation
      lwd = 2
      
      shown = risks$marker >= wtd.quantile(dat[[vx]], dat$wt, 2.5 /
                                             100) &
        risks$marker <= wtd.quantile(dat[[vx]], dat$wt, 1 -
                                       2.5 / 100)
      
      # hard code ylim to make the plot look better
      ylim = c(0, 0.11)
      #ylim=range(risks$lb[shown,], risks$ub[shown,], 0) # [shown] so that there is not too much empty space
      xlim = get.range.cor(dat, marker.name.to.assay(vx), tpeak)
      if (verbose >= 2)
        myprint(xlim, ylim)
      
      # set up an empty plot
      plot(
        risks$marker[shown],
        risks$prob[shown, 1],
        xlab = paste0(markers.names.short[vx], " (=s)"),
        ylab = paste0(
          "Probability* of ",
          config.cor$txt.endpoint,
          " by Day ",
          tfinal.tpeak
        ),
        lwd = lwd,
        xlim = xlim,
        ylim = ylim,
        type = "n",
        main = "",
        xaxt = "n"
      )
      title(main = for.title,
            line = .6,
            cex = .9)
      draw.x.axis.cor(xlim, lloxs[vx], llox_labels[vx])
      
      # draw risk lines and confidence bands
      for (i in 1:length(risks$marker.2)) {
        lines(
          risks$marker[shown],
          risks$prob[shown, i],
          lwd = lwd,
          col = i,
          lty = ifelse(i == 2, 2, 1)
        )# use dashed line for the middle so that overlaps can be seen
        lines(
          risks$marker[shown],
          risks$lb[shown, i],
          lwd = lwd,
          col = i,
          lty = 3
        )
        lines(
          risks$marker[shown],
          risks$ub[shown, i],
          lwd = lwd,
          col = i,
          lty = 3
        )
      }
      
      # legend for the three lines
      legend.txt = c("(15th percentile)", "(median)", "(85th percentile)")
      #                # special handling code
      #                if (attr(config, "config")=="hvtn705second") {
      #                    if (inner.id==1) legend.txt=c("(min)","(median)","(90th percentile)")
      #                }
      mylegend(
        x = 3,
        legend = paste(signif(10 ** risks$marker.2, 3), legend.txt),
        col = 1:3,
        lty = c(1, 2, 1),
        title = markers.names.short[vthree],
        lwd = lwd
      )
      
      # placebo prevelance lines
      abline(
        h = prev.plac[1],
        col = "gray",
        lty = c(1, 3, 3),
        lwd = lwd
      )
      text(
        x = par("usr")[2] - diff(par("usr")[1:2]) / 5,
        y = prev.plac[1] + diff(par("usr")[3:4]) / 30,
        "placebo arm " %.% formatDouble(prev.plac[1], 3, remove.leading0 = F)
      )
      #abline(h=risks.itxn.1$prob[1,1], col="gray", lty=c(1), lwd=lwd)
      #text(x=par("usr")[2]-diff(par("usr")[1:2])/5, y=risks.itxn.1$prob[1,1]+diff(par("usr")[3:4])/30, "placebo arm "%.%formatDouble(risks.itxn.1$prob[1,1],3,remove.leading0=F))
      
      # add histogram
      par(new = TRUE)
      col <-
        c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
      col <-
        rgb(col[1],
            col[2],
            col[3],
            alpha = 255 * 0.4,
            maxColorValue = 255)
      tmp.x = dat[[vx]][dat$ph2 == 1]
      tmp.w = dat$wt[dat$ph2 == 1]
      tmp = get.marker.histogram(tmp.x, tmp.w, attr(config, "config"))
      if (is.nan(tmp$density)[1])
        tmp = hist(tmp.x, plot = F)
      plot(
        tmp,
        col = col,
        axes = F,
        labels = F,
        main = "",
        xlab = "",
        ylab = "",
        border = 0,
        freq = F,
        xlim = xlim,
        ylim = c(0, max(tmp$density * 1.25))
      )
      
      dev.off()
    }
  }
}


}
