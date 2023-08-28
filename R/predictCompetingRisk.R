predictCompetingRisk=function(formula, formula.all, data, t0, newdata=data, stype=2, ctype=2, ...){
    
    # weights might be in the optional arguments
    extra.args <- list(...)
    
    # fit a cause 1-specific model to get hazard
    fit.1=do.call("coxph",c(list(formula=formula, data=data, model=TRUE),extra.args))
    #direct calling coxph cannot handle weights
    #fit.1 <- coxph(formula=formula, data=data, model=TRUE)
    
    # fit a all-cause model to get survival prob
    fit.a=do.call("coxph",c(list(formula=formula.all, data=data, model=TRUE),extra.args))
    # cannot directly call coxph since that cannot handle weights
    #fit.a <- coxph(formula=formula.all, data=data, model=TRUE)
    
    X.1=model.matrix(formula, newdata)[,-1,drop=FALSE]# -1 in formula does not work
    X.a=model.matrix(formula.all, newdata)[,-1,drop=FALSE]# -1 in formula does not work
    
    if (ncol(X.1)>0) {
        F.1=drop(exp(X.1 %*% coef(fit.1)))
        F.a=drop(exp(X.a %*% coef(fit.a)))
    } else {
        F.1=rep(1,nrow(X.1))
        F.a=rep(1,nrow(X.a))
    }
    
    #### get hazard from cause-specific model
    
#    bhaz.1=basehaz(fit.1, centered=F)# stype=2 and ctype=2 when calling basehaz
#    names(bhaz.1)[1]="cumhaz"
    
    # mimic basehaz, this allows more flexibility since stype and ctype can be changed
    # somehow stype also affects cumhaz
    sfit=survival::survfit(fit.1, se.fit = FALSE, stype=stype, ctype=ctype)
    # centered=FALSE requires the following
    zcoef <- ifelse(is.na(coef(fit.1)), 0, coef(fit.1))
    offset <- sum(fit.1$means * zcoef)
    bhaz.1 <- data.frame(time = sfit$time, cumhaz = sfit$cumhaz * exp(-offset), surv=sfit$surv ^ exp(-offset))
    #plot(sfit$surv, exp(-chaz)); abline(0,1)
    strata <- sfit$strata
    if (!is.null(strata)) bhaz.1$strata <- factor(rep(names(strata), strata), levels = names(strata))
    
    # compute hazard from cumhazard
    bhaz.1=cbind(bhaz.1, hazard=c(bhaz.1$cumhaz[1], diff(bhaz.1$cumhaz)))

    # time points
    tt=bhaz.1$time[bhaz.1$time<=t0]
    # idx is used to subset to where hazard is not 0
    idx=which(bhaz.1$hazard[1:length(tt)]!=0)        

    # hazard from cause-specific model
    h=outer(bhaz.1$hazard[1:length(tt)][idx], F.1)# dim: n_times x n_subj
    #print(h[,1])

    
    #### get survival prob from all-cauase model
    
#    bhaz.a=basehaz(fit.a, centered=F)# stype=2 and ctype=2 when calling basehaz
#    names(bhaz.a)[1]="cumhaz"
    
    # mimic basehaz, this allows more flexibility since stype and ctype can be changed
    sfit=survival::survfit(fit.a, se.fit = FALSE, stype=stype, ctype=ctype)
    zcoef <- ifelse(is.na(coef(fit.a)), 0, coef(fit.a))
    offset <- sum(fit.a$means * zcoef)
    bhaz.a <- data.frame(time = sfit$time, cumhaz = sfit$cumhaz * exp(-offset), surv=sfit$surv ^ exp(-offset))
    #plot(sfit$surv, exp(-chaz)); abline(0,1)
    strata <- sfit$strata
    if (!is.null(strata)) bhaz.a$strata <- factor(rep(names(strata), strata), levels = names(strata))
    
    # survival prob from all-cause model
    if(stype==1) {
        # prod limit
        S.1.mat=outer(c(1,bhaz.a$surv[1:(length(tt)-1)])[idx], F.a, "^")# dim: n_times x n_subj
    } else {
        # exp(cum hazard)
        S.1.mat=exp(-outer(c(0,bhaz.a$cumhaz[1:(length(tt)-1)])[idx], F.a))# dim: n_times x n_subj
    }
    #print(head(cbind(t(S.1.mat), t(S.1.mat.2))))
    #plot(S.1.mat, S.1.mat.2)
    #print(bhaz.a$cumhaz[1:3]*F.a[1])
    
    
    #### cumulative incidence 
    colSums(S.1.mat * h)# dim: n_subj
    
}
pcr=predictCompetingRisk
