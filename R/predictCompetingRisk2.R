predictCompetingRisk2=function(formula.list, data, t0, newdata=data, ...){
    
    # weights might be in the optional arguments
    extra.args <- list(...)
    
    fits=lapply (formula.list, function(formula) {
        # fit a cause-specific model to get hazard
        do.call("coxph",c(list(formula=formula, data=data, model=TRUE),extra.args))
        #direct calling coxph cannot handle weights
        #coxph(formula=formula, data=data, model=TRUE)
    })
    
    bhazs=lapply(fits, function(fit) {
        ret=survival::basehaz(fit, centered=F)
        names(ret)[1]="cumhaz"
        ret
    }) # stype=2 and ctype=2 when calling basehaz
    
    Fs=lapply (1:length(formula.list), function(i) {
        X=model.matrix(formula.list[[i]], newdata)[,-1,drop=FALSE]# -1 in formula does not work
        if (ncol(X)>0) {
            F=drop(exp(X %*% coef(fits[[i]])))
            if (all(is.na(F))) {
              # this could happen if there are no cases of this type
              # set all values to a value other than NA
              F=rep(1,nrow(X)) 
            }
        } else {
            F=rep(1,nrow(X))
        }
        F
    })
		
    if (all(Fs[[1]]==1)) {
			stop("return because there are no cases of interest")
	  }

    # assume the first is the cause of interest
    bhaz.1=bhazs[[1]]

    # compute hazard from cumhazard
    bhaz.1=cbind(bhaz.1, hazard=c(bhaz.1$cumhaz[1], diff(bhaz.1$cumhaz)))

    # time points
    tt=bhaz.1$time[bhaz.1$time<=t0]
        # idx is used to subset to where hazard is not 0
        idx=which(bhaz.1$hazard[1:length(tt)]!=0)        

    # get hazard from cause-specific model
    h=outer(bhaz.1$hazard[1:length(tt)][idx], Fs[[1]])# dim: n_times x n_subj
    #print(h[,1])
    
    # get survival prob from each model    
    mat=lapply (1:length(formula.list), function (i) {    
        outer(c(0,bhazs[[i]][1:(length(tt)-1),1])[idx], Fs[[i]])
    })    
    S.1.mat=exp(-do.call("+", mat))# dim: n_times x n_subj
    
    
    # cumulative incidence 
    cum.ind = S.1.mat * h # dim: n_subj
  	out = colSums(cum.ind)
  	attr(out, "cumulative") <- apply(cum.ind, 2, cumsum)
  	attr(out, "time") <- tt[idx]

	out    
}
pcr2=predictCompetingRisk2
