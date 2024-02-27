
# bootstrap case cohort samples
# data is assumed to contain only ph1 ptids
get.bootstrap.data.cor = function(data, ptids.by.stratum, seed) {
  set.seed(seed)    
  
  # For each sampling stratum, bootstrap samples in subcohort and not in subchort separately
  tmp=lapply(ptids.by.stratum, function(x) c(sample(x$subcohort, replace=TRUE), sample(x$nonsubcohort, replace=TRUE)))
  
  dat.b=data[match(unlist(tmp), data$Ptid),]
  
  # compute weights
  tmp=with(dat.b, table(Wstratum, ph2))
  weights=rowSums(tmp)/tmp[,2]
  dat.b$wt=weights[""%.%dat.b$Wstratum]
  # we assume data only contains ph1 ptids, thus weights is defined for every bootstrapped ptids
  
  dat.b
}


# bootstrap for COVE boost
# Within each quadrant (2 Trt * 2 naive status):
#   1. resample the cohort and count the number of cases and controls: n1 and n0
#         if n1 < 32 or n0 < 32, redo
#   2. resample 32 cases and 32 controls from ph2 samples 
#   3. resample n1-32 cases and n0-32 controls from non-ph2 samples
# 4. Collapse strata if needed and compute inverse probability sampling weights
# Thus, the number of cases may vary across bootstrap replicates, but the ph2 sample size remains constant

# dat.ph1 is ph1 data and need to have, in addition to markers and covariates columns:
#   Ptid, Trt, naive, EventIndPrimary, ph2, demo.stratum, CalendarBD1Interval
# return a dataframe with wt column

# bootstrap.cove.boost.2 is faster than bootstrap.cove.boost and results are similar
bootstrap.cove.boost=function(dat.ph1, seed) {
  
  set.seed(seed)
  
  # perform bootstrap within each quadrant (2 Trt * 2 naive status)
  dat.b=NULL
  for (idat in 1:4) {
    if (idat==1) {dat.tmp = subset(dat.ph1, dat.ph1$Trt==1 & dat.ph1$naive==1)}
    if (idat==2) {dat.tmp = subset(dat.ph1, dat.ph1$Trt==0 & dat.ph1$naive==1)}
    if (idat==3) {dat.tmp = subset(dat.ph1, dat.ph1$Trt==1 & dat.ph1$naive==0)}
    if (idat==4) {dat.tmp = subset(dat.ph1, dat.ph1$Trt==0 & dat.ph1$naive==0)}
    if (nrow(dat.tmp)==0) next
    
    dat.tmp.nph2=subset(dat.tmp, !dat.tmp$ph2)
    dat.tmp.ph2=subset(dat.tmp, dat.tmp$ph2)
    
    # n1.ph2 and n0.ph2 are expected to be 32 in COVE Boost
    # we make it data-dependent here to be more flexible
    n1.ph2 = sum(dat.tmp.ph2$EventIndPrimary)
    n0.ph2 = sum(1-dat.tmp.ph2$EventIndPrimary)
    
    # 1. 
    dat.2=dat.tmp[sample.int(nrow(dat.tmp), replace=TRUE),]
    n1 = nrow(subset(dat.2, dat.2$EventIndPrimary==1))
    n0 = nrow(subset(dat.2, dat.2$EventIndPrimary==0))
    
    while(n1<n1.ph2 | n0<n0.ph2) {   
      dat.2=dat.tmp[sample.int(nrow(dat.tmp), replace=TRUE),]
      n1 = nrow(subset(dat.2, dat.2$EventIndPrimary==1))
      n0 = nrow(subset(dat.2, dat.2$EventIndPrimary==0))
    }
    
    # 2.
    dat.ph2.cases=subset(dat.tmp.ph2, dat.tmp.ph2$EventIndPrimary==1)
    dat.ph2.cases.b=dat.ph2.cases[sample.int(nrow(dat.ph2.cases), size=n1.ph2, replace=TRUE),]
    
    dat.ph2.ctrls=subset(dat.tmp.ph2, dat.tmp.ph2$EventIndPrimary==0)
    dat.ph2.ctrls.b=dat.ph2.ctrls[sample.int(nrow(dat.ph2.ctrls), size=n0.ph2, replace=TRUE),]
    
    # 3.
    dat.nph2.cases=subset(dat.tmp.nph2, dat.tmp.nph2$EventIndPrimary==1)
    dat.nph2.cases.b=dat.nph2.cases[sample.int(nrow(dat.nph2.cases), size=n1-n1.ph2, replace=TRUE),]
    
    dat.nph2.ctrls=subset(dat.tmp.nph2, dat.tmp.nph2$EventIndPrimary==0)
    dat.nph2.ctrls.b=dat.nph2.ctrls[sample.int(nrow(dat.nph2.ctrls), size=n0-n0.ph2, replace=TRUE),]
    
    dat.b=rbind(dat.b, dat.ph2.cases.b, dat.ph2.ctrls.b, dat.nph2.cases.b, dat.nph2.ctrls.b)
  }
  
  
  # 4. 
  n.demo = length(table(dat.b$demo.stratum))
  # adjust Wstratum
  ret = cove.boost.collapse.strata (dat.b, n.demo)
  
  # compute inverse probability sampling weights
  tmp = with(ret, ph1)
  wts_table <- with(ret[tmp,], table(Wstratum, ph2))
  wts_norm <- rowSums(wts_table) / wts_table[, 2]
  ret[["wt"]] = ifelse(ret$ph1, wts_norm[ret$Wstratum %.% ""], NA)
  
  assertthat::assert_that(
    all(!is.na(subset(ret, ret$tmp & !is.na(ret$Wstratum))[["wt"]])),
    msg = "missing wt.BD for D analyses ph1 subjects")
  
  return (ret)
}




# a second version, simpler, faster, results are close to bootstrap.cove.boost
bootstrap.cove.boost.2=function(dat.ph1, seed) {
  
  set.seed(seed)
  
  # perform bootstrap within each quadrant (2 Trt * 2 naive status)
  dat.b=NULL
  for (idat in 1:4) {
    if (idat==1) {dat.tmp = subset(dat.ph1, dat.ph1$Trt==1 & dat.ph1$naive==1)}
    if (idat==2) {dat.tmp = subset(dat.ph1, dat.ph1$Trt==0 & dat.ph1$naive==1)}
    if (idat==3) {dat.tmp = subset(dat.ph1, dat.ph1$Trt==1 & dat.ph1$naive==0)}
    if (idat==4) {dat.tmp = subset(dat.ph1, dat.ph1$Trt==0 & dat.ph1$naive==0)}
    if (nrow(dat.tmp)==0) next
    
    dat.b=dat.tmp[sample.int(nrow(dat.tmp), replace=TRUE),]
  }
  
  # 4. adjust Wstratum
  n.demo = length(table(dat.b$demo.stratum))
  ret = cove.boost.collapse.strata (dat.b, n.demo)
  
  # compute inverse probability sampling weights
  tmp = with(ret, ph1)
  wts_table <- with(ret[tmp,], table(Wstratum, ph2))
  wts_norm <- rowSums(wts_table) / wts_table[, 2]
  ret[["wt"]] = ifelse(ret$ph1, wts_norm[ret$Wstratum %.% ""], NA)
  
  assertthat::assert_that(
    all(!is.na(subset(ret, ret$tmp & !is.na(ret$Wstratum))[["wt"]])),
    msg = "missing wt.BD for D analyses ph1 subjects")
  
  return (ret)
}




# bootstrap from case control studies is done by resampling cases, ph2 controls, and non-ph2 controls separately. 

# when all cases are sampled, bootstrap from case control studies is done by resampling cases, ph2 controls, and non-ph2 controls separately. 
# e.g. hvtn705

# Across bootstrap replicates, the number of cases does not stay constant, neither do the numbers of ph2 controls by demographics strata. 
# Specifically,
# 1) sample with replacement to get dat.b. From this dataset, take the cases and count ph2 and non-ph2 controls by strata
# 2) sample with replacement ph2 and non-ph2 controls by strata
bootstrap.case.control.samples=function(dat.ph1, seed, delta.name="EventIndPrimary", strata.name="tps.stratum", ph2.name="ph2", min.cell.size=1) {
  #dat.ph1=dat.tmp; delta.name="EventIndPrimary"; strata.name="tps.stratum"; ph2.name="ph2"; min.cell.size=0
  
  set.seed(seed)
  
  dat.tmp=data.frame(ptid=1:nrow(dat.ph1), delta=dat.ph1[,delta.name], strata=dat.ph1[,strata.name], ph2=dat.ph1[,ph2.name])
  
  nn.ph1=with(dat.tmp, table(strata, delta))
  strat=rownames(nn.ph1); names(strat)=strat
  # ctrl.ptids is a list of lists
  ctrl.ptids = with(subset(dat.tmp, dat.tmp$delta==0), lapply(strat, function (i) list(ph2=ptid[strata==i & ph2], nonph2=ptid[strata==i & !ph2])))
  
  # 1. resample dat.ph1 to get dat.b, but only take the cases 
  dat.b=dat.tmp[sample.int(nrow(dat.tmp), replace=TRUE),]
  
  # re-do resampling if the bootstrap dataset has too few samples in a cell in nn.ctrl.b
  while(TRUE) {   
    nn.ctrl.b=with(subset(dat.b, !dat.b$delta), table(strata, ph2))
    if (min(nn.ctrl.b)<min.cell.size | ncol(nn.ctrl.b)<2) dat.b=dat.tmp[sample.int(nrow(dat.tmp), replace=TRUE),] else break
  }
  
  # take the case ptids
  case.ptids.b = dat.b$ptid[dat.b$delta==1]
  
  # 2. resample controls in dat.ph1 (numbers determined by dat.b) stratified by strata and ph2/nonph2
  # ph2 and non-ph2 controls by strata
  nn.ctrl.b=with(subset(dat.b, !dat.b$delta), table(strata, ph2))
  # sample the control ptids
  ctrl.ptids.by.stratum.b=lapply(strat, function (i) {
    c(sample(ctrl.ptids[[i]]$ph2, nn.ctrl.b[i,2], replace=TRUE),
      sample(ctrl.ptids[[i]]$nonph2, nn.ctrl.b[i,1], replace=TRUE))
  })
  ctrl.ptids.b=do.call(c, ctrl.ptids.by.stratum.b)    
  
  # return data frame
  dat.ph1[c(case.ptids.b, ctrl.ptids.b), ]
}

## testing
#dat.b=bootstrap.case.control.samples(dat.vac.seroneg)
#with(dat.vac.seroneg, table(ph2, tps.stratum, EventIndPrimary))
#with(dat.b, table(ph2, tps.stratum, EventIndPrimary))
#> with(dat.vac.seroneg, table(ph2, tps.stratum, EventIndPrimary))
#, , EventIndPrimary = 0
#
#       tps.stratum
#ph2       33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48
#  FALSE 1483  915  759  439 1677 1138  894  591 3018 1973 1559 1051 1111  693  511  329
#  TRUE    57   53   55   57   56   57   57   56   58   55   55   57   57   56   56   56
#
#, , EventIndPrimary = 1
#
#       tps.stratum
#ph2       33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48
#  FALSE    1    0    0    1    0    1    0    0    2    1    2    1    0    0    0    1
#  TRUE     3    7    7   10    8   11    2   13   17   23   15   23    5    6    4    6
#
#> with(dat.b, table(ph2, tps.stratum, EventIndPrimary))
#, , EventIndPrimary = 0
#
#       tps.stratum
#ph2       33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48
#  FALSE 1487  911  750  462 1675 1181  884  570 3058 2023 1499 1034 1094  694  487  329
#  TRUE    47   57   65   62   50   53   50   64   55   61   65   53   64   53   54   60
#
#, , EventIndPrimary = 1
#
#       tps.stratum
#ph2       33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48
#  FALSE    0    0    0    0    0    2    0    0    1    1    3    3    0    0    0    2
#  TRUE     2    6    8    5    9   13    0   11   20   26   10   20    4    3    4    5


# for bootstrap use
get.ptids.by.stratum.for.bootstrap = function(data) {
  strat=sort(unique(data$tps.stratum))
  ptids.by.stratum=lapply(strat, function (i) 
    list(subcohort=subset(data, data$tps.stratum==i & data$SubcohortInd==1)[["Ptid"]], 
         nonsubcohort=subset(data, data$tps.stratum==i & data$SubcohortInd==0)[['Ptid']])
  )    
  # add a pseudo-stratum for subjects with NA in tps.stratum (not part of Subcohort). 
  # we need this group because it contains some cases with missing tps.stratum
  # if data is ph2 only, then this group is only cases because ph2 = subcohort + cases
  tmp=list(subcohort=subset(data, is.na(data$tps.stratum))[['Ptid']],               nonsubcohort=NULL)
  ptids.by.stratum=append(ptids.by.stratum, list(tmp))    
  ptids.by.stratum
}




bootstrap.cohort=function(dat, seed) {
  
  set.seed(seed)
  
  dat.b=dat[sample.int(nrow(dat), replace=TRUE),]
  
  dat.b

}