# Output:
# table to report HR, CI and P value of the last nCoef coefficients. One model per row

# nCoef=1, used for model with one continuous effect of interest
# nCoef=2, used for model with two continuous effects of interest
# nCoef=3, used for interaction between two continuous models
# no multitesting b/c we usually do multitesting across continuous and discrete variables
# all.markers, which defines the formula and can be flexible, e.g., Bmarkercat * scale(Day15marker,scale=F)

cor_coxph_coef_n = function(
  form.0,
  design_or_dat, # either design when tps is T, or data frame otherwise
  fname.suffix, #used in the file names to save results
  save.results.to,
  config,
  config.cor,
  all.markers,
  all.markers.names.short, 
  
  nCoef,
  col.headers, # length n coef, e.g., c("TrtA Moderna~Pfizer", "center(D15)", "TrtA:center(D15)")
  
  verbose=FALSE
) {
  
  if(verbose) {
    cat("Running cor_coxph_coef_n \n")
    myprint(col.headers)
  }

  tps=inherits(design_or_dat, "survey.design")
  if (tps) dat = design_or_dat$phase1$full$variables

  fits=list()
  for (a in all.markers) {
    f= update(form.0, as.formula(paste0("~.+", a)))
    if (tps) {
      fits[[a]]=svycoxph(f, design=design_or_dat) 
    } else {
      fits[[a]]=try(coxph(f, design_or_dat) )
      if (inherits(fits[[a]], "try-error")) {
        fits[[a]]=NULL
        all.markers.names.short = all.markers.names.short[names(all.markers.names.short)!=a]
      }
    }
  }
  
  
  rows=length(coef(fits[[1]]))
  rows=(rows-(nCoef-1)):rows
  # robust=tps b/c not an option to use robust=T for coxph, but it is a required argument for getFormattedSummary
  est=getFormattedSummary(fits, exp=T, robust=tps, rows=rows, type=1)
  ci= getFormattedSummary(fits, exp=T, robust=tps, rows=rows, type=7)
  p=  getFormattedSummary(fits, exp=T, robust=tps, rows=rows, type=10)
  
  
  natrisk=nrow(dat)
  nevents=sum(dat$yy==1)

  tab.1=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), 
              est[1,], ci[1,], p[1,])
  if (nCoef>1) {
    for (i in 2:nCoef) tab.1=cbind(tab.1, 
                                   est[i,], ci[i,], p[i,])
  }
  rownames(tab.1)=all.markers.names.short
  tab.1
  
  header=paste0("\\hline\n 
       \\multicolumn{1}{l}{", '', "} & \\multicolumn{1}{c}{}", concatList(paste0("& \\multicolumn{3}{c}{",col.headers,"} ")), "   \\\\ 
       \\multicolumn{1}{l}{", '', "} & \\multicolumn{1}{c}{No. cases /}   ", concatList(rep("& \\multicolumn{2}{c}{HR per 10-fold incr.} & \\multicolumn{1}{c}{P-value}",nCoef)), "  \\\\ 
       \\multicolumn{1}{l}{Immunologic Marker}  & \\multicolumn{1}{c}{No. at-risk**} ", concatList(rep("& \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{}", nCoef)), "  \\\\ 
       \\hline\n 
  ")
  
  # keep svycoxph in the name so that we can reuse rmd file
  mytex(tab.1, file.name="CoR_univariable_svycoxph_pretty_"%.%fname.suffix, align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
        col.headers=header,
        longtable=T, 
        label=paste0("tab:CoR_univariable_svycoxph_pretty"), 
        caption.placement = "top", 
        caption=paste0("Inference for Day ", config.cor$tpeak, " antibody marker covariate-adjusted correlates of risk of ", config.cor$txt.endpoint, " in the ", escape(fname.suffix), " group: Hazard ratios per 10-fold increment in the marker*")
  )
  
  
  pvals.cont = sapply(fits, function(x) {
    tmp=getFixedEf(x)
    p.val.col=which(startsWith(tolower(colnames(tmp)),"p"))
    tmp[nrow(tmp),p.val.col]
  })
  
  pvals.cont
}

