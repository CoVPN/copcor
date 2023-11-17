##This file includes the code used to run the simulations in the paper.

## functions for upper and lower limits of MI-Wilson interval
#Upper limit function for the MI-Wilson interval for multiple imputation
upper.limit.new = function(z, n, r.m, mean_q) {
  ul = ((((2*mean_q) + ((z^2)/n) + (((z^2) * r.m)/n))/
           (2*(1 + ((z^2)/n) + (((z^2)*r.m)/n)))) +
          sqrt(((((2*mean_q) + ((z^2)/n) + (((z^2)*r.m)/n))^2)/
                  (4*(1 + ((z^2)/n) + (((z^2)*r.m)/n))^2)) - 
                 ((mean_q^2)/(1 + ((z^2)/n) + (((z^2)*r.m)/n)))))
  return(ul)
}

#Lower limit function for the MI-Wilson interval for multiple imputation
lower.limit.new = function(z, n, r.m, mean_q) {
  ll = ((((2*mean_q) + ((z^2)/n) + (((z^2) * r.m)/n))/
           (2*(1 + ((z^2)/n) + (((z^2)*r.m)/n)))) -
          sqrt(((((2*mean_q) + ((z^2)/n) + (((z^2)*r.m)/n))^2)/
                  (4*(1 + ((z^2)/n) + (((z^2)*r.m)/n))^2)) - 
                 ((mean_q^2)/(1 + ((z^2)/n) + (((z^2)*r.m)/n)))))
  return(ll)
}

## functions for limits for MI-Plug intervals
#Center for the plug-in Wilson method for multiple imputation
plugin.center = function(z, n, mean_q){
  center = 
    (mean_q+((1/(2*n))* (z^2)))/(1+((1/n)*(z^2)))
  
  return(center)
}

#Spread for the plug-in Wilson method for multiple imputation
plugin.spread = function(z, n, mean_q){
  spread =
    ((z*sqrt(((1/n)*mean_q*(1-mean_q))+
               ((1/(4*(n^2)))* (z^2))))/
       (1+((1/n)*(z^2))))
  
  return(spread)
}

## functions for upper and lower limits of MI-Li intervals
#Function for upper limit of Li, Mehrotra and Barnard interval
li.et.al.ul = function(mean_q, Tm, z){
  ul = (mean_q + ((z^2)/(2*((mean_q*(1-mean_q))/Tm))) + 
          (z*sqrt(((mean_q*(1-mean_q))/((mean_q*(1-mean_q))/Tm)) + 
                    ((z^2)/(4*((mean_q*(1-mean_q))/Tm)^2)))))/
    (1 + ((z^2)/((mean_q*(1-mean_q))/Tm)))
  return(ul)
}

#Function for lower limit of Li, Mehrotra and Barnard interval
li.et.al.ll = function(mean_q, Tm, z){
  ll = (mean_q + ((z^2)/(2*((mean_q*(1-mean_q))/Tm))) -
          (z*sqrt(((mean_q*(1-mean_q))/((mean_q*(1-mean_q))/Tm)) + 
                    ((z^2)/(4*((mean_q*(1-mean_q))/Tm)^2)))))/
    (1 + ((z^2)/((mean_q*(1-mean_q))/Tm)))
  return(ll)
}

#function for complete data Wilson interval limits
#Upper limit function for the standard Wilson interval
wilson.upper.limit = function(p.hat, z, n){
  ul = (((p.hat + ((1/(2*n))*(z^2)))/(1 + ((1/n)*(z^2)))) + 
          (z*sqrt(((1/n)*p.hat*(1-p.hat)) + ((1/(4*(n^2)))*(z^2))))/(1+((1/n)*(z^2))))
  return(ul)
    
}

#Lower limit function for the standard Wilson interval
wilson.lower.limit = function(p.hat, z, n){
  ll = (((p.hat + ((1/(2*n))*(z^2)))/(1 + ((1/n)*(z^2)))) -
          (z*sqrt(((1/n)*p.hat*(1-p.hat)) + ((1/(4*(n^2)))*(z^2))))/(1+((1/n)*(z^2))))
  return(ll)
  
}


## functions for an interval that is no longer used in the paper -- kept in the script to document entire process used in coding
#Upper limit function for obsolete Wilson interval for multiple imputation
upper.limit = function(z, m, n, r, mean_q){
  ul = (((2*mean_q +((z^2)/n) + (((z^2)*r)/n) + (((z^2)*r)/(m*n)))/
           (2*(1+ ((z^2)/n) + (((z^2)*r)/n) + (((z^2)*r)/(m*n))))) +
          sqrt(((((-2*mean_q)- ((z^2)/n)- (((z^2)*r)/n)- (((z^2)*r)/(m*n)))^2)/
                  (4*(1+((z^2)/n) + (((z^2)*r)/n)+ (((z^2)*r)/(m*n)))^2))-
                 ((mean_q^2)/(1+ ((z^2)/n) + (((z^2)*r)/n)+(((z^2)*r)/(m*n))))))
  return(ul)
}

#Lower limit function for obsolete Wilson interval for multiple imputation
lower.limit = function(z, m, n, r, mean_q){
  ll = (((2*mean_q +((z^2)/n) + (((z^2)*r)/n) + (((z^2)*r)/(m*n)))/
           (2*(1+ ((z^2)/n) + (((z^2)*r)/n) + (((z^2)*r)/(m*n))))) -
          sqrt(((((-2*mean_q)- ((z^2)/n)- (((z^2)*r)/n)- (((z^2)*r)/(m*n)))^2)/
                  (4*(1+((z^2)/n) + (((z^2)*r)/n)+ (((z^2)*r)/(m*n)))^2))-
                 ((mean_q^2)/(1+ ((z^2)/n) + (((z^2)*r)/n)+(((z^2)*r)/(m*n))))))
  return(ll)
}


####  Code to evaluate different confidence intervals with multiple  imputation using repeated sampling

mult.imp = function(n, p, a, b, it, m, mp){
  #n is the number of binomial trials
  #p is the true proportion of successes
  #a is the first beta prior parameter
  #b is the second beta prior parameter
  #it is the number of iterations 
  #m is the number of imputed data sets
  #mp is the proportion of missingness for
  #MCAR data
  
  #Initializing vectors
  qm.v = rep(NA, it)
  um.v = rep(NA, it)
  bm.v = rep(NA, it)
  Tm.v = rep(NA, it)
  v.v = rep(NA, it)
  ll.wz.v = rep(NA, it)
  ul.wz.v = rep(NA,it)
  ll.wt.v = rep(NA,it)
  ul.wt.v = rep(NA, it)
  ll.wald.v = rep(NA,it)
  ul.wald.v = rep(NA, it)
  ll.wn.v = rep(NA, it)
  ul.wn.v = rep(NA, it)
  ul.wm.v = rep(NA, it)
  ll.wm.v = rep(NA, it)
  ul.li.v = rep(NA, it)
  ll.li.v = rep(NA, it)
  p.in.int.wz = rep(NA, it)
  p.in.int.wt = rep(NA, it)
  p.in.int.wald = rep(NA, it)
  p.in.int.wn = rep(NA, it)
  p.in.int.wm = rep(NA, it)
  p.in.int.li = rep(NA, it)
  nonzero.width.wt = rep(NA, it)
  nonzero.width.wz = rep(NA, it)
  nonzero.width.wn = rep(NA, it)
  nonzero.width.wald = rep(NA, it)
  nonzero.width.wm = rep(NA, it)
  nonzero.width.li = rep(NA, it)
  int.length.wt = rep(NA,it)
  int.length.wz = rep(NA, it)
  int.length.wald = rep(NA, it)
  int.length.wn = rep(NA, it)
  int.length.wm = rep(NA, it)
  int.length.li = rep(NA, it)
  int.length.wald.c = rep(NA, it)
  int.length.wn.c = rep(NA, it)
  int.length.wt.c = rep(NA, it)
  int.length.wz.c = rep(NA, it)
  int.length.wm.c = rep(NA, it)
  int.length.li.c = rep(NA, it)
  prop.miss.v = rep(NA, it)
  
  qm.v.mcar = rep(NA, it)
  um.v.mcar = rep(NA, it)
  bm.v.mcar = rep(NA, it)
  Tm.v.mcar = rep(NA, it)
  v.v.mcar = rep(NA, it)
  ll.wz.v.mcar = rep(NA, it)
  ul.wz.v.mcar = rep(NA,it)
  ll.wt.v.mcar = rep(NA,it)
  ul.wt.v.mcar = rep(NA, it)
  ll.wald.v.mcar = rep(NA,it)
  ul.wald.v.mcar = rep(NA, it)
  ll.wn.v.mcar = rep(NA, it)
  ul.wn.v.mcar = rep(NA, it)
  ul.wm.v.mcar = rep(NA, it)
  ll.wm.v.mcar = rep(NA, it)
  ul.li.v.mcar = rep(NA, it)
  ll.li.v.mcar = rep(NA, it)
  p.in.int.wz.mcar = rep(NA, it)
  p.in.int.wt.mcar = rep(NA, it)
  p.in.int.wald.mcar = rep(NA, it)
  p.in.int.wn.mcar = rep(NA, it)
  p.in.int.wm.mcar = rep(NA, it)
  p.in.int.li.mcar = rep(NA, it)
  nonzero.width.wt.mcar = rep(NA, it)
  nonzero.width.wz.mcar = rep(NA, it)
  nonzero.width.wn.mcar = rep(NA, it)
  nonzero.width.wald.mcar = rep(NA, it)
  nonzero.width.wm.mcar = rep(NA, it)
  nonzero.width.li.mcar = rep(NA, it)
  int.length.wt.mcar = rep(NA,it)
  int.length.wz.mcar = rep(NA, it)
  int.length.wald.mcar = rep(NA, it)
  int.length.wn.mcar = rep(NA, it)
  int.length.wm.mcar = rep(NA, it)
  int.length.li.mcar = rep(NA, it)
  int.length.wald.c.mcar = rep(NA, it)
  int.length.wn.c.mcar = rep(NA, it)
  int.length.wt.c.mcar = rep(NA, it)
  int.length.wz.c.mcar = rep(NA, it)
  int.length.wm.c.mcar = rep(NA, it)
  int.length.li.c.mcar = rep(NA, it)
  
  int.length.wt.nzw = rep(NA,it)
  int.length.wz.nzw = rep(NA, it)
  int.length.wald.nzw = rep(NA, it)
  int.length.wn.nzw = rep(NA, it)
  int.length.wm.nzw = rep(NA, it)
  int.length.li.nzw = rep(NA, it)
  int.length.wald.c.nzw = rep(NA, it)
  int.length.wn.c.nzw = rep(NA, it)
  int.length.wt.c.nzw = rep(NA, it)
  int.length.wz.c.nzw = rep(NA, it)
  int.length.wm.c.nzw = rep(NA, it)
  int.length.li.c.nzw = rep(NA, it)
  int.length.wt.mcar.nzw = rep(NA,it)
  int.length.wz.mcar.nzw = rep(NA, it)
  int.length.wald.mcar.nzw = rep(NA, it)
  int.length.wn.mcar.nzw = rep(NA, it)
  int.length.wm.mcar.nzw = rep(NA, it)
  int.length.li.mcar.nzw = rep(NA, it)
  int.length.wald.c.mcar.nzw = rep(NA, it)
  int.length.wn.c.mcar.nzw = rep(NA, it)
  int.length.wt.c.mcar.nzw = rep(NA, it)
  int.length.wz.c.mcar.nzw = rep(NA, it)
  int.length.wm.c.mcar.nzw = rep(NA, it)
  int.length.li.c.mcar.nzw = rep(NA, it)
  pz.bm.um = rep(NA, it)
  pz.bm.um.mcar = rep(NA, it)
  wilson.ll = rep(NA, it)
  wilson.ul = rep(NA, it)
  nonzero.width.wilson = rep(NA, it)
  p.in.int.wilson = rep(NA, it)
  int.length.wilson = rep(NA, it)
  int.length.wilson.c = rep(NA, it)
  int.length.wilson.nzw = rep(NA,it)
  int.length.wilson.c.nzw = rep(NA, it)
  s.wald.ll = rep(NA, it)
  s.wald.ul = rep(NA, it)
  nonzero.width.s.wald = rep(NA, it)
  p.in.int.s.wald = rep(NA, it)
  int.length.s.wald = rep(NA, it)
  int.length.s.wald.c = rep(NA, it)
  int.length.s.wald.nzw = rep(NA, it)
  int.length.s.wald.c.nzw = rep(NA, it)
  
  for (i in 1:it){
    
    #Initialize x.d, generate y.d
    x.d = rep(NA, n)
    y.d = rbinom(n, 1, p)
    mcar.data = y.d
    no.missing.data = y.d
    #Generate x.d based on y.d
    for (u in 1:length(y.d)){
      if (y.d[u] ==1){
        x.d[u] = rbinom(1, 1, 0.6)
      }else {
        x.d[u] = rbinom(1,1, 0.6)
      }
    }
    
    #Create vector of indices indicating
    #missingness based on x.d
    miss = rep(NA, n)
    for (k in 1:length(x.d)){
      if (x.d[k] == 1){
        miss[k] = rbinom(1,1, 0.06) 
      }else{                       
        miss[k] = rbinom(1,1, 0.16)
      }
    }
    
    
    #Seperately compiling all observed y.d's 
    #corresponding to x.d = 1 and x.d = 0 (where
    #missingness index equals 0)
    y_x1 = c()
    y_x0 = c()
    for (j2 in 1:length(x.d)){
      if (miss[j2] == 0 & x.d[j2] ==1){
        y_x1 = c(y_x1, y.d[j2])
      }else if (miss[j2] == 0 & x.d[j2] == 0){
        y_x0 = c(y_x0, y.d[j2])
      }
    }
    
    #Determining r, the proportion of missingness
    count = 0 
    for (xi in 1:length(miss)){
      if (miss[xi] ==1){
        count = count + 1
      }
    }
    prop_miss = count/n
    
    prop.miss.v[i] = prop_miss
    
    #Initializing vectors for multiple imputation
    p_hat = rep(NA, m)
    var_p = rep(NA, m)
    #Below: creating m multiply imputed data sets,
    #each data set with missing y.d imputed based on 
    #whether x.d = 1 or x.d = 0 for that index,
    #with different imputation parameters for the 
    #two conditions
    an_x1 = a + sum(y_x1)
    bn_x1 = b + length(y_x1) - sum(y_x1)
    an_x0 = a + sum(y_x0)
    bn_x0 = b + length(y_x0) - sum(y_x0)
    for (j in 1:m) {
      theta_x1 = rbeta(1, an_x1, bn_x1)
      theta_x0 = rbeta(1, an_x0, bn_x0)
      for (k2 in 1:length(y.d)){
        if (miss[k2] ==1 & x.d[k2] == 1){
          y.d[k2] = rbinom(1, 1, theta_x1)
        }else if (miss[k2] == 1 & x.d[k2] == 0){
          y.d[k2] = rbinom(1,1, theta_x0)
        }
      }
      p_hat[j] = mean(y.d)
      var_p[j] = p_hat[j]*(1-p_hat[j])/n
    }
    #Rubin's rules for MAR data
    qm= sum(p_hat)/m
    um = sum(var_p)/m
    bm = var(p_hat)
    Tm = ((1 + (1/m))*bm) + um
    if (bm == 0){
      v = 1000000000
      rm = 0 #convert to rm = 0, for both bm = 0 and um = 0 and 
      #bm = 0 and um != 0
    }else{
      v = (m-1)*(1+(um/((1+ (1/m))*bm)))^2 #convert to normal if bm = 0 and um != 0
      rm = (1 + (1/m))*(bm/um)
    }
    
    #Multiple Imputation for MCAR data
    s = n*(1-mp)
    s_plus_1 = s + 1
    missing_vec = c(s_plus_1:n)
    obs = mcar.data[-missing_vec]
    an.mcar = a + sum(obs)
    bn.mcar = b + length(obs) - sum(obs)
    p.hat.mcar = rep(NA, m)
    var.p.mcar = rep(NA, m)
    for (ji in 1:m) {
      theta.mcar = rbeta(1, an.mcar, bn.mcar)
      imp.mcar = rbinom(n*mp, 1, theta.mcar)
      binom_data = c(obs, imp.mcar)
      p.hat.mcar[ji] = mean(binom_data)
      var.p.mcar[ji] = p.hat.mcar[ji]*(1-p.hat.mcar[ji])/n
    }
    #Rubin's rules for MCAR data
    qm.mcar= sum(p.hat.mcar)/m
    um.mcar = sum(var.p.mcar)/m
    bm.mcar = var(p.hat.mcar)
    Tm.mcar = ((1 + (1/m))*bm.mcar) + um.mcar
    if (bm.mcar == 0){
      v.mcar = 1000000000
      rm.mcar = 0
    }else{
      v.mcar = (m-1)*(1+(um.mcar/((1+ (1/m))*bm.mcar)))^2
      rm.mcar = (1 + (1/m))*(bm.mcar/um.mcar)
    }
    
    #Store values of bm, um, bm.mcar, and um.mcar to check if they are equal to 0.
    bm.v[i] = bm
    um.v[i] = um
    bm.v.mcar[i] = bm.mcar
    um.v.mcar[i] = um.mcar
    
    #Check if both bm and um, or both bm.mcar and um.mcar, are equal to 0 at the 
    #same time
    if (bm.v[i] == 0 & um.v[i] == 0){
      pz.bm.um[i] = TRUE
    }else{
      pz.bm.um[i] = FALSE
    }
    
    if (bm.v.mcar[i] == 0 & um.v.mcar[i] == 0){
      pz.bm.um.mcar[i] = TRUE
    }else{
      pz.bm.um.mcar[i] = FALSE
    }
    
    #Compute the lower limits and upper limits for each interval based on MAR data
    
    ll.wald.v[i] = qm - qt(0.975,df = v)*sqrt(Tm)
    ul.wald.v[i] = qm + qt(0.975,df = v)*sqrt(Tm)
    ll.wn.v[i] = lower.limit.new(z = qt(0.975, df = v), n = n, r.m = rm, mean_q = qm)
    ul.wn.v[i] = upper.limit.new(z = qt(0.975, df = v), n = n, r.m = rm , mean_q = qm)
    ll.wz.v[i] = (plugin.center(z = qnorm(0.975), n = n, mean_q = qm) -
                    plugin.spread(z = qnorm(0.975), n = n, mean_q = qm))
    ul.wz.v[i] = (plugin.center(z = qnorm(0.975), n = n, mean_q = qm) +
                    plugin.spread(z = qnorm(0.975), n = n, mean_q = qm))
    ll.wt.v[i] = (plugin.center(z = qt(0.975, df = v), n = n, mean_q = qm) -
                    plugin.spread(z = qt(0.975, df = v), n = n, mean_q = qm))
    ul.wt.v[i] = (plugin.center(z = qt(0.975, df = v), n = n, mean_q = qm) +
                    plugin.spread(z = qt(0.975, df = v), n = n, mean_q = qm))
    ll.wm.v[i] = lower.limit(z = qt(0.975, df = v), m = m, n = n, r = prop.miss.v[i],
                             mean_q = qm)
    ul.wm.v[i] = upper.limit(z = qt(0.975, df = v), m = m, n = n, r = prop.miss.v[i],
                             mean_q = qm)
    ll.li.v[i] = li.et.al.ll(mean_q = qm, Tm = Tm, z = qnorm(0.975))
    ul.li.v[i] = li.et.al.ul(mean_q =qm, Tm = Tm, z = qnorm(0.975))
    
    #Resetting NaNs in the Li et al. interval to 0 or 1 to signify 
    #zero-width intervals
    if (is.na(ll.li.v[i]) == TRUE &
        is.na(ul.li.v[i]) == TRUE){
      if (qm == 1 & Tm == 0){
        ll.li.v[i] = 1
        ul.li.v[i] = 1
      }else{
        ll.li.v[i] = 0
        ul.li.v[i] = 0
      }
    }
    
    #Compute the lower limits and upper limits of each interval based on MCAR data
    ll.wald.v.mcar[i] = qm.mcar - qt(0.975, df = v.mcar)*sqrt(Tm.mcar)
    ul.wald.v.mcar[i] = qm.mcar + qt(0.975, df = v.mcar)*sqrt(Tm.mcar)
    ll.wn.v.mcar[i] = lower.limit.new(z = qt(0.975, df = v.mcar),
                                      n = n, r.m = rm.mcar, mean_q = qm.mcar)
    ul.wn.v.mcar[i] = upper.limit.new(z = qt(0.975, df = v.mcar), 
                                      n = n, r.m = rm.mcar, mean_q = qm.mcar)
    ll.wz.v.mcar[i] = (plugin.center(z = qnorm(0.975), n = n, mean_q = qm.mcar) -
                    plugin.spread(z = qnorm(0.975), n = n, mean_q = qm.mcar))
    ul.wz.v.mcar[i] = (plugin.center(z = qnorm(0.975), n = n, mean_q = qm.mcar) +
                    plugin.spread(z = qnorm(0.975), n = n, mean_q = qm.mcar))
    ll.wt.v.mcar[i] = (plugin.center(z = qt(0.975, df = v.mcar), n = n, mean_q = qm.mcar) -
                    plugin.spread(z = qt(0.975, df = v.mcar), n = n, mean_q = qm.mcar))
    ul.wt.v.mcar[i] = (plugin.center(z = qt(0.975, df = v.mcar), n = n, mean_q = qm.mcar) +
                    plugin.spread(z = qt(0.975, df = v.mcar), n = n, mean_q = qm.mcar))
    ll.wm.v.mcar[i] = lower.limit(z = qnorm(0.975), 
                                  m = m, n = n, r = 0.3, mean_q = qm.mcar)
    ul.wm.v.mcar[i] = upper.limit(z = qnorm(0.975),
                                  m = m, n = n, r = 0.3, mean_q = qm.mcar)
    ul.li.v.mcar[i] = li.et.al.ul(mean_q = qm.mcar, Tm = Tm.mcar, z = qnorm(0.975))
    ll.li.v.mcar[i] = li.et.al.ll(mean_q = qm.mcar, Tm = Tm.mcar, z = qnorm(0.975))
    
    #Resetting NaNs in the Li et al. interval to 0 or 1 to signal 
    #zero-width intervals
    if (is.na(ul.li.v.mcar[i]) == TRUE & 
        is.na(ll.li.v.mcar[i]) == TRUE){
      if (qm.mcar == 1 & Tm.mcar == 0){
        ul.li.v.mcar[i] = 1
        ll.li.v.mcar[i] = 1
      }else{
        ul.li.v.mcar[i] = 0
        ll.li.v.mcar[i] = 0
      }
    }
    
    #Compute the interval lengths without truncation for the MAR data for each interval
    int.length.wt[i] = ul.wt.v[i] - ll.wt.v[i]
    int.length.wz[i] = ul.wz.v[i] -ll.wz.v[i]
    int.length.wald[i] = ul.wald.v[i] - ll.wald.v[i]
    int.length.wn[i] = ul.wn.v[i] - ll.wn.v[i]
    int.length.wm[i] = ul.wm.v[i] - ll.wm.v[i]
    int.length.li[i] = ul.li.v[i] - ll.li.v[i]
    
    #Compute the interval lengths without truncation for the MAR data for each interval, 
    #setting zero-width interval lengths to NA
    if (int.length.wt[i] != 0){
      int.length.wt.nzw[i] = int.length.wt[i]
    }else{
      int.length.wt.nzw[i] = NA
    }
    if (int.length.wz[i] != 0){
      int.length.wz.nzw[i] = int.length.wz[i]
    }else{
      int.length.wz.nzw[i] = NA
    }
    if (int.length.wald[i] != 0){
      int.length.wald.nzw[i] = int.length.wald[i]
    }else{
      int.length.wald.nzw[i] = NA
    }
    if (int.length.wn[i] != 0){
      int.length.wn.nzw[i] = int.length.wn[i]
    }else{
      int.length.wn.nzw[i] = NA
    }
    if (int.length.wm[i] != 0){
      int.length.wm.nzw[i] = int.length.wm[i]
    }else{
      int.length.wm.nzw[i] = NA
    }
    if (int.length.li[i] != 0){
      int.length.li.nzw[i] = int.length.li[i]
    }else{
      int.length.li.nzw[i] = NA
    }
    
    #Compute the upper and lower limits of the Wilson and Wald intervals using the
    #initial data set, without any missing data
    p_hat.2 = mean(no.missing.data)
    wilson.ll[i] = wilson.lower.limit(z = qnorm(0.975), n = n, p.hat = p_hat.2)
    wilson.ul[i] = wilson.upper.limit(z = qnorm(0.975), n = n, p.hat = p_hat.2)
    s.wald.ul[i] = p_hat.2 + (qnorm(0.975)*sqrt((1/n)*p_hat.2*(1-p_hat.2)))
    s.wald.ll[i] = p_hat.2 - (qnorm(0.975)*sqrt((1/n)*p_hat.2*(1-p_hat.2)))
    
    #Store the number of times the width of the intervals computed using no missing 
    #data are zero width
    nonzero.width.wilson[i] = (wilson.ul[i] - wilson.ll[i] > 0)
    nonzero.width.s.wald[i] = (s.wald.ul[i] - s.wald.ll[i] > 0)
    
    #Determine whether or not p is in the intervals created without missing data
    p.in.int.s.wald[i] = (s.wald.ul[i] >= p & p >= s.wald.ll[i])
    p.in.int.wilson[i] = (wilson.ul[i] >= p & p >= wilson.ll[i])
    
    #Compute the interval lengths without truncation for no missing data, 
    #including the zero-width intervals
    int.length.s.wald[i] = s.wald.ul[i] - s.wald.ll[i]
    int.length.wilson[i] = wilson.ul[i] - wilson.ll[i]
    
    #Compute the truncated interval lengths for no missing data, including zero-width
    #intervals
    if (wilson.ul[i] >1 & wilson.ll[i] >= 0){
      int.length.wilson.c[i] = 1 - wilson.ll[i]
    }else if (wilson.ll[i] < 0 & wilson.ul[i] <=1){
      int.length.wilson.c[i] = wilson.ul[i]
    }else if (wilson.ll[i] < 0 & wilson.ul[i] > 1){
      int.length.wilson.c[i] = 1
    }else{
      int.length.wilson.c[i] = wilson.ul[i] - wilson.ll[i]
    }
    
    if (s.wald.ul[i] >1 & s.wald.ll[i] >= 0){
      int.length.s.wald.c[i] = 1 - s.wald.ll[i]
    }else if (s.wald.ll[i] < 0 & s.wald.ul[i] <=1){
      int.length.s.wald.c[i] = s.wald.ul[i]
    }else if (s.wald.ll[i] < 0 & s.wald.ul[i] > 1){
      int.length.s.wald.c[i] = 1
    }else{
      int.length.s.wald.c[i] = s.wald.ul[i] - s.wald.ll[i]
    }
    
    #Compute the interval lengths without truncation for the MCAR data for each interval
    int.length.wt.mcar[i] = ul.wt.v.mcar[i] - ll.wt.v.mcar[i]
    int.length.wz.mcar[i] = ul.wz.v.mcar[i] -ll.wz.v.mcar[i]
    int.length.wald.mcar[i] = ul.wald.v.mcar[i] - ll.wald.v.mcar[i]
    int.length.wn.mcar[i] = ul.wn.v.mcar[i] - ll.wn.v.mcar[i]
    int.length.wm.mcar[i] = ul.wm.v.mcar[i] - ll.wm.v.mcar[i]
    int.length.li.mcar[i] = ul.li.v.mcar[i] - ll.li.v.mcar[i]
    
    #Compute the interval lengths without truncation for the MCAR data for each interval,
    #setting zero-width interval lengths to NA
    if (int.length.wt.mcar[i] != 0){
      int.length.wt.mcar.nzw[i] = int.length.wt.mcar[i]
    }else{
      int.length.wt.mcar.nzw[i] = NA
    }
    if (int.length.wz.mcar[i] != 0){
      int.length.wz.mcar.nzw[i] = int.length.wz.mcar[i]
    }else{
      int.length.wz.mcar.nzw[i] = NA
    }
    if (int.length.wald.mcar[i] != 0){
      int.length.wald.mcar.nzw[i] = int.length.wald.mcar[i]
    }else{
      int.length.wald.mcar.nzw[i] = NA
    }
    if (int.length.wn.mcar[i] != 0){
      int.length.wn.mcar.nzw[i] = int.length.wn.mcar[i]
    }else{
      int.length.wn.mcar.nzw[i] = NA
    }
    if (int.length.wm.mcar[i] != 0){
      int.length.wm.mcar.nzw[i] = int.length.wm.mcar[i]
    }else{
      int.length.wm.mcar.nzw[i] = NA
    }
    if (int.length.li.mcar[i] != 0){
      int.length.li.mcar.nzw[i] = int.length.li.mcar[i]
    }else{
      int.length.li.mcar.nzw[i] = NA
    }
    
    #Create a logical vector checking whether p is in the interval for each type 
    #of interval, using MAR data
    p.in.int.wald[i] = (ll.wald.v[i] <= p & ul.wald.v[i] >= p)
    p.in.int.wm[i] = (ll.wm.v[i] <= p & ul.wm.v[i] >= p)
    p.in.int.wn[i] = (ll.wn.v[i] <= p & ul.wn.v[i] >= p)
    p.in.int.wt[i] = (ll.wt.v[i] <= p & ul.wt.v[i] >= p)
    p.in.int.wz[i] = (ll.wz.v[i] <= p & ul.wz.v[i] >= p)
    p.in.int.li[i] = (ll.li.v[i] <= p & ul.li.v[i] >= p)
    
    #Create a logical vector checking whether p is in the interval for each type 
    #of interval, using MCAR data
    p.in.int.wald.mcar[i] = (ll.wald.v.mcar[i] <= p & ul.wald.v.mcar[i] >= p)
    p.in.int.wm.mcar[i] = (ll.wm.v.mcar[i] <= p & ul.wm.v.mcar[i] >= p)
    p.in.int.wn.mcar[i] = (ll.wn.v.mcar[i] <= p & ul.wn.v.mcar[i] >= p)
    p.in.int.wt.mcar[i] = (ll.wt.v.mcar[i] <= p & ul.wt.v.mcar[i] >= p)
    p.in.int.wz.mcar[i] = (ll.wz.v.mcar[i] <= p & ul.wz.v.mcar[i] >= p)
    p.in.int.li.mcar[i] = (ll.li.v.mcar[i] <= p & ul.li.v.mcar[i] >= p)
    
    #Create a logical vector checking whether each type of interval is zero-width 
    #for this iteration, using MAR data
    nonzero.width.wald[i] = (ul.wald.v[i] - ll.wald.v[i] > 0)
    nonzero.width.wt[i] = (ul.wt.v[i] - ll.wt.v[i] > 0)
    nonzero.width.wz[i] = (ul.wz.v[i] - ll.wz.v[i] > 0)
    nonzero.width.wm[i] = (ul.wm.v[i] - ll.wm.v[i] > 0)
    nonzero.width.wn[i] = (ul.wn.v[i] - ll.wn.v[i] > 0)
    nonzero.width.li[i] = (ul.li.v[i] - ll.li.v[i] > 0)
    
    #Create a logical vector checking whether each type of interval is zero-width 
    #for this iteration, using MCAR data
    nonzero.width.wald.mcar[i] = (ul.wald.v.mcar[i] - ll.wald.v.mcar[i] > 0)
    nonzero.width.wt.mcar[i] = (ul.wt.v.mcar[i] - ll.wt.v.mcar[i] > 0)
    nonzero.width.wz.mcar[i] = (ul.wz.v.mcar[i] - ll.wz.v.mcar[i] > 0)
    nonzero.width.wm.mcar[i] = (ul.wm.v.mcar[i] - ll.wm.v.mcar[i] > 0)
    nonzero.width.wn.mcar[i] = (ul.wn.v.mcar[i] - ll.wn.v.mcar[i] > 0)
    nonzero.width.li.mcar[i] = (ul.li.v.mcar[i] - ll.li.v.mcar[i] > 0)
    
    #Truncate the interval lengths by setting the upper and lower interval limits to 1 and 0 
    #respectively, for MAR data, including zero-width intervals
    if (ul.wt.v[i] >1 & ll.wt.v[i] >= 0){
      int.length.wt.c[i] = 1 - ll.wt.v[i]
    }else if (ll.wt.v[i] < 0 & ul.wt.v[i] <=1){
      int.length.wt.c[i] = ul.wt.v[i]
    }else if (ll.wt.v[i] < 0 & ul.wt.v[i] > 1){
      int.length.wt.c[i] = 1
    }else{
      int.length.wt.c[i] = ul.wt.v[i] - ll.wt.v[i]
    }
    
    if (ul.wz.v[i] >1 & ll.wz.v[i] >= 0){
      int.length.wz.c[i] = 1 - ll.wz.v[i]
    }else if (ll.wz.v[i] < 0 & ul.wz.v[i] <=1){
      int.length.wz.c[i] = ul.wz.v[i]
    }else if (ll.wz.v[i] < 0 & ul.wz.v[i] > 1){
      int.length.wz.c[i] = 1
    }else{
      int.length.wz.c[i] = ul.wz.v[i] - ll.wz.v[i]
    }
    
    if (ul.wald.v[i] > 1 & ll.wald.v[i] >= 0){
      int.length.wald.c[i] = 1 - ll.wald.v[i]
    }else if (ll.wald.v[i] < 0 & ul.wald.v[i] <= 1){
      int.length.wald.c[i] = ul.wald.v[i]
    }else if (ll.wald.v[i] < 0 & ul.wald.v[i] > 1){
      int.length.wald.c[i] = 1
    }else{
      int.length.wald.c[i] = ul.wald.v[i] - ll.wald.v[i]
    }
    
    if (ul.wn.v[i] > 1 & ll.wn.v[i] >= 0){
      int.length.wn.c[i] = 1 - ll.wn.v[i]
    }else if (ll.wn.v[i] < 0 & ul.wn.v[i] <= 1){
      int.length.wn.c[i] = ul.wn.v[i]
    }else if (ll.wn.v[i] < 0 & ul.wn.v[i] > 1){
      int.length.wn.c[i] = 1
    }else{
      int.length.wn.c[i] = ul.wn.v[i] - ll.wn.v[i]
    }
    
    if (ul.wm.v[i] > 1 & ll.wm.v[i] >= 0){
      int.length.wm.c[i] = 1 - ll.wm.v[i]
    }else if (ll.wm.v[i] < 0 & ul.wm.v[i] <= 1){
      int.length.wm.c[i] = ul.wm.v[i]
    }else if (ll.wm.v[i] < 0 & ul.wm.v[i] > 1){
      int.length.wm.c[i] = 1
    }else{
      int.length.wm.c[i] = ul.wm.v[i] - ll.wm.v[i]
    }
    
    if (ul.li.v[i] > 1 & ll.li.v[i] >= 0){
      int.length.li.c[i] = 1 - ll.li.v.[i]
    }else if (ll.li.v[i] < 0 & ul.li.v[i] <= 1){
      int.length.li.c[i] = ul.li.v[i]
    }else if (ll.li.v[i] < 0 & ul.li.v[i] > 1){
      int.length.li.c[i] = 1
    }else{
      int.length.li.c[i] = ul.li.v[i] - ll.li.v[i]
    }
    
    #Truncate the interval lengths by setting the upper and lower interval limits to 1 and 0
    #respectively, for MAR data, while also setting zero-width interval lengths to NA
    if (int.length.wt[i] != 0){
      if (ul.wt.v[i] >1 & ll.wt.v[i] >= 0){
        int.length.wt.c.nzw[i] = 1 - ll.wt.v[i]
      }else if (ll.wt.v[i] < 0 & ul.wt.v[i] <=1){
        int.length.wt.c.nzw[i] = ul.wt.v[i]
      }else if (ll.wt.v[i] < 0 & ul.wt.v[i] > 1){
        int.length.wt.c.nzw[i] = 1
      }else{
        int.length.wt.c.nzw[i] = ul.wt.v[i] - ll.wt.v[i]
      }
    }else{
      int.length.wt.c.nzw[i] = NA
    }
    
    if (int.length.wz[i] != 0){
      if (ul.wz.v[i] >1 & ll.wz.v[i] >= 0){
        int.length.wz.c.nzw[i] = 1 - ll.wz.v[i]
      }else if (ll.wz.v[i] < 0 & ul.wz.v[i] <=1){
        int.length.wz.c.nzw[i] = ul.wz.v[i]
      }else if (ll.wz.v[i] < 0 & ul.wz.v[i] > 1){
        int.length.wz.c.nzw[i] = 1
      }else{
        int.length.wz.c.nzw[i] = ul.wz.v[i] - ll.wz.v[i]
      }
    }else{
      int.length.wz.c.nzw[i] = NA
    }
    
    if (int.length.wald[i] != 0){
      if (ul.wald.v[i] > 1 & ll.wald.v[i] >= 0){
        int.length.wald.c.nzw[i] = 1 - ll.wald.v[i]
      }else if (ll.wald.v[i] < 0 & ul.wald.v[i] <= 1){
        int.length.wald.c.nzw[i] = ul.wald.v[i]
      }else if (ll.wald.v[i] < 0 & ul.wald.v[i] > 1){
        int.length.wald.c.nzw[i] = 1
      }else{
        int.length.wald.c.nzw[i] = ul.wald.v[i] - ll.wald.v[i]
      }
    }else{
      int.length.wald.c.nzw[i] = NA
    }
    
    if (int.length.wn[i] != 0){
      if (ul.wn.v[i] > 1 & ll.wn.v[i] >= 0){
        int.length.wn.c.nzw[i] = 1 - ll.wn.v[i]
      }else if (ll.wn.v[i] < 0 & ul.wn.v[i] <= 1){
        int.length.wn.c.nzw[i] = ul.wn.v[i]
      }else if (ll.wn.v[i] < 0 & ul.wn.v[i] > 1){
        int.length.wn.c.nzw[i] = 1
      }else{
        int.length.wn.c.nzw[i] = ul.wn.v[i] - ll.wn.v[i]
      }
    }else{
      int.length.wn.c.nzw[i] = NA
    }
    
    if (int.length.wm[i] != 0){
      if (ul.wm.v[i] > 1 & ll.wm.v[i] >= 0){
        int.length.wm.c.nzw[i] = 1 - ll.wm.v[i]
      }else if (ll.wm.v[i] < 0 & ul.wm.v[i] <= 1){
        int.length.wm.c.nzw[i] = ul.wm.v[i]
      }else if (ll.wm.v[i] < 0 & ul.wm.v[i] > 1){
        int.length.wm.c.nzw[i] = 1
      }else{
        int.length.wm.c.nzw[i] = ul.wm.v[i] - ll.wm.v[i]
      }
    }else{
      int.length.wm.c.nzw[i] = NA
    }
    
    if (int.length.li[i] != 0) {
      if (ul.li.v[i] > 1 & ll.li.v[i] >= 0){
        int.length.li.c.nzw[i] = 1 - ll.li.v[i]
      }else if (ll.li.v[i] < 0 & ul.li.v[i] <= 1){
        int.length.li.c.nzw[i] = ul.li.v[i]
      }else if (ll.li.v[i] < 0 & ul.li.v[i] > 1){
        int.length.li.c.nzw[i] = 1
      }else{
        int.length.li.c.nzw[i] = ul.li.v[i] - ll.li.v[i]
      }
    }else{
      int.length.li.c.nzw[i] = NA
    }
    
    #Truncate the interval lengths by setting the lower and upper limits to 0 and 1,
    #respectively, for MCAR data, including zero-width intervals
    if (ul.wt.v.mcar[i] >1 & ll.wt.v.mcar[i] >= 0){
      int.length.wt.c.mcar[i] = 1 - ll.wt.v.mcar[i]
    }else if (ll.wt.v.mcar[i] < 0 & ul.wt.v.mcar[i] <=1){
      int.length.wt.c.mcar[i] = ul.wt.v.mcar[i]
    }else if (ll.wt.v.mcar[i] < 0 & ul.wt.v.mcar[i] > 1){
      int.length.wt.c.mcar[i] = 1
    }else{
      int.length.wt.c.mcar[i] = ul.wt.v.mcar[i] - ll.wt.v.mcar[i]
    }
    
    if (ul.wz.v.mcar[i] >1 & ll.wz.v.mcar[i] >= 0){
      int.length.wz.c.mcar[i] = 1 - ll.wz.v.mcar[i]
    }else if (ll.wz.v.mcar[i] < 0 & ul.wz.v.mcar[i] <=1){
      int.length.wz.c.mcar[i] = ul.wz.v.mcar[i]
    }else if (ll.wz.v.mcar[i] < 0 & ul.wz.v.mcar[i] > 1){
      int.length.wz.c.mcar[i] = 1
    }else{
      int.length.wz.c.mcar[i] = ul.wz.v.mcar[i] - ll.wz.v.mcar[i]
    }
    
    if (ul.wald.v.mcar[i] > 1 & ll.wald.v.mcar[i] >= 0){
      int.length.wald.c.mcar[i] = 1 - ll.wald.v.mcar[i]
    }else if (ll.wald.v.mcar[i] < 0 & ul.wald.v.mcar[i] <= 1){
      int.length.wald.c.mcar[i] = ul.wald.v.mcar[i]
    }else if (ll.wald.v.mcar[i] < 0 & ul.wald.v.mcar[i] > 1){
      int.length.wald.c.mcar[i] = 1
    }else{
      int.length.wald.c.mcar[i] = ul.wald.v.mcar[i] - ll.wald.v.mcar[i]
    }
    
    if (ul.wn.v.mcar[i] > 1 & ll.wn.v.mcar[i] >= 0){
      int.length.wn.c.mcar[i] = 1 - ll.wn.v.mcar[i]
    }else if (ll.wn.v.mcar[i] < 0 & ul.wn.v.mcar[i] <= 1){
      int.length.wn.c.mcar[i] = ul.wn.v.mcar[i]
    }else if (ll.wn.v.mcar[i] < 0 & ul.wn.v.mcar[i] > 1){
      int.length.wn.c.mcar[i] = 1
    }else{
      int.length.wn.c.mcar[i] = ul.wn.v.mcar[i] - ll.wn.v.mcar[i]
    }
    
    if (ul.wm.v.mcar[i] > 1 & ll.wm.v.mcar[i] >= 0){
      int.length.wm.c.mcar[i] = 1 - ll.wm.v.mcar[i]
    }else if (ll.wm.v.mcar[i] < 0 & ul.wm.v.mcar[i] <= 1){
      int.length.wm.c.mcar[i] = ul.wm.v.mcar[i]
    }else if (ll.wm.v.mcar[i] < 0 & ul.wm.v.mcar[i] > 1){
      int.length.wm.c.mcar[i] = 1
    }else{
      int.length.wm.c.mcar[i] = ul.wm.v.mcar[i] - ll.wm.v.mcar[i]
    }
    
    if (ul.li.v.mcar[i] > 1 & ll.li.v.mcar[i] >= 0){
      int.length.li.c.mcar[i] = 1 - ll.li.v.mcar[i]
    }else if (ll.li.v.mcar[i] < 0 & ul.li.v.mcar[i] <= 1){
      int.length.li.c.mcar[i] = ul.li.v.mcar[i]
    }else if (ll.li.v.mcar[i] < 0 & ul.li.v.mcar[i] > 1){
      int.length.li.c.mcar[i] = 1
    }else{
      int.length.li.c.mcar[i] = ul.li.v.mcar[i] - ll.li.v.mcar[i]
    }
    
    #Truncate the interval lengths by setting the upper and lower interval limits to 1 and 0
    #respectively, for MCAR data, while also setting zero-width interval lengths to NA
    if (int.length.wt.mcar[i] != 0){
      if (ul.wt.v.mcar[i] >1 & ll.wt.v.mcar[i] >= 0){
        int.length.wt.c.mcar.nzw[i] = 1 - ll.wt.v.mcar[i]
      }else if (ll.wt.v.mcar[i] < 0 & ul.wt.v.mcar[i] <=1){
        int.length.wt.c.mcar.nzw[i] = ul.wt.v.mcar[i]
      }else if (ll.wt.v.mcar[i] < 0 & ul.wt.v.mcar[i] > 1){
        int.length.wt.c.mcar.nzw[i] = 1
      }else{
        int.length.wt.c.mcar.nzw[i] = ul.wt.v.mcar[i] - ll.wt.v.mcar[i]
      }
    }else{
      int.length.wt.c.mcar.nzw[i] = NA
    }
    
    if (int.length.wz.mcar[i] != 0){
      if (ul.wz.v.mcar[i] >1 & ll.wz.v.mcar[i] >= 0){
        int.length.wz.c.mcar.nzw[i] = 1 - ll.wz.v.mcar[i]
      }else if (ll.wz.v.mcar[i] < 0 & ul.wz.v.mcar[i] <=1){
        int.length.wz.c.mcar.nzw[i] = ul.wz.v.mcar[i]
      }else if (ll.wz.v.mcar[i] < 0 & ul.wz.v.mcar[i] > 1){
        int.length.wz.c.mcar.nzw[i] = 1
      }else{
        int.length.wz.c.mcar.nzw[i] = ul.wz.v.mcar[i] - ll.wz.v.mcar[i]
      }
    }else{
      int.length.wz.c.mcar.nzw[i] = NA
    }
    
    if (int.length.wald.mcar[i] != 0){
      if (ul.wald.v.mcar[i] > 1 & ll.wald.v.mcar[i] >= 0){
        int.length.wald.c.mcar.nzw[i] = 1 - ll.wald.v.mcar[i]
      }else if (ll.wald.v.mcar[i] < 0 & ul.wald.v.mcar[i] <= 1){
        int.length.wald.c.mcar.nzw[i] = ul.wald.v.mcar[i]
      }else if (ll.wald.v.mcar[i] < 0 & ul.wald.v.mcar[i] > 1){
        int.length.wald.c.mcar.nzw[i] = 1
      }else{
        int.length.wald.c.mcar.nzw[i] = ul.wald.v.mcar[i] - ll.wald.v.mcar[i]
      }
    }else{
      int.length.wald.c.mcar.nzw[i] = NA
    }
    
    if (int.length.wn.mcar[i] != 0){
      if (ul.wn.v.mcar[i] > 1 & ll.wn.v.mcar[i] >= 0){
        int.length.wn.c.mcar.nzw[i] = 1 - ll.wn.v.mcar[i]
      }else if (ll.wn.v.mcar[i] < 0 & ul.wn.v.mcar[i] <= 1){
        int.length.wn.c.mcar.nzw[i] = ul.wn.v.mcar[i]
      }else if (ll.wn.v.mcar[i] < 0 & ul.wn.v.mcar[i] > 1){
        int.length.wn.c.mcar.nzw[i] = 1
      }else{
        int.length.wn.c.mcar.nzw[i] = ul.wn.v.mcar[i] - ll.wn.v.mcar[i]
      }
    }else{
      int.length.wn.c.mcar.nzw[i] = NA
    }
    
    if (int.length.wm.mcar[i] != 0){
      if (ul.wm.v.mcar[i] > 1 & ll.wm.v.mcar[i] >= 0){
        int.length.wm.c.mcar.nzw[i] = 1 - ll.wm.v.mcar[i]
      }else if (ll.wm.v.mcar[i] < 0 & ul.wm.v.mcar[i] <= 1){
        int.length.wm.c.mcar.nzw[i] = ul.wm.v.mcar[i]
      }else if (ll.wm.v.mcar[i] < 0 & ul.wm.v.mcar[i] > 1){
        int.length.wm.c.mcar.nzw[i] = 1
      }else{
        int.length.wm.c.mcar.nzw[i] = ul.wm.v.mcar[i] - ll.wm.v.mcar[i]
      }
    }else{
      int.length.wm.c.mcar.nzw[i] = NA
    }
    
    if (int.length.li.mcar[i] != 0){
      if (ul.li.v.mcar[i] > 1 & ll.li.v.mcar[i] >= 0){
        int.length.li.c.mcar.nzw[i] = 1 - ll.li.v.mcar[i]
      }else if (ll.li.v.mcar[i] < 0 & ul.li.v.mcar[i] <= 1){
        int.length.li.c.mcar.nzw[i] = ul.li.v.mcar[i]
      }else if (ll.li.v.mcar[i] < 0 & ul.li.v.mcar[i] > 1){
        int.length.li.c.mcar.nzw[i] = 1
      }else{
        int.length.li.c.mcar.nzw[i] = ul.li.v.mcar[i] - ll.li.v.mcar[i]
      }
    }else{
      int.length.li.c.mcar.nzw[i] = NA
    }
    
    #Set zero-width intervals equal to NA for no missing data
    if (int.length.s.wald[i] != 0){
      int.length.s.wald.nzw[i] = int.length.s.wald[i]
    }else{
      int.length.s.wald.nzw[i] = NA
    }
    if (int.length.wilson[i] != 0){
      int.length.wilson.nzw[i] = int.length.wilson[i]
    }else{
      int.length.wilson.nzw[i] = NA
    }
    
    #For no missing data, calculate interval lengths by truncating
    #intervals by cutting intervals at 0 and 1, and set zero-width intervals to NA
    if (int.length.s.wald[i] != 0){
      if (s.wald.ul[i] >1 & s.wald.ll[i] >= 0){
        int.length.s.wald.c.nzw[i] = 1 - s.wald.ll[i]
      }else if (s.wald.ll[i] < 0 & s.wald.ul[i] <=1){
        int.length.s.wald.c.nzw[i] = s.wald.ul[i]
      }else if (s.wald.ll[i] < 0 & s.wald.ul[i] > 1){
        int.length.s.wald.c.nzw[i] = 1
      }else{
        int.length.s.wald.c.nzw[i] = s.wald.ul[i] - s.wald.ll[i]
      }
    }else{
      int.length.s.wald.c.nzw[i] = NA
    }
    
    if (int.length.wilson[i] != 0){
      if (wilson.ul[i] >1 & wilson.ll[i] >= 0){
        int.length.wilson.c.nzw[i] = 1 - wilson.ll[i]
      }else if (wilson.ll[i] < 0 & wilson.ul[i] <=1){
        int.length.wilson.c.nzw[i] = wilson.ul[i]
      }else if (wilson.ll[i] < 0 & wilson.ul[i] > 1){
        int.length.wilson.c.nzw[i] = 1
      }else{
        int.length.wilson.c.nzw[i] = wilson.ul[i] - wilson.ll[i]
      }
    }else{
      int.length.wilson.c.nzw[i] = NA
    }
    
    
  }
  
  #Compute the mean truncated and untruncated interval lengths for MAR data,
  #including zero-width intervals
  mean.int.length.wt = mean(int.length.wt)
  mean.int.length.wz = mean(int.length.wz)
  mean.int.length.wn = mean(int.length.wn)
  mean.int.length.wald = mean(int.length.wald)
  mean.int.length.wm = mean(int.length.wm)
  mean.int.length.li = mean(int.length.li)
  mean.int.length.wald.c = mean(int.length.wald.c)
  mean.int.length.wt.c = mean(int.length.wt.c)
  mean.int.length.wn.c = mean(int.length.wn.c)
  mean.int.length.wz.c = mean(int.length.wz.c)
  mean.int.length.wm.c = mean(int.length.wm.c)
  mean.int.length.li.c = mean(int.length.li.c)
  
  #Compute the mean truncated and untruncated interval lengths for MAR data,
  #leaving out zero-width intervals
  mean.int.length.wt.nzw = mean(int.length.wt.nzw, na.rm = TRUE)
  mean.int.length.wz.nzw = mean(int.length.wz.nzw, na.rm = TRUE)
  mean.int.length.wn.nzw = mean(int.length.wn.nzw, na.rm = TRUE)
  mean.int.length.wald.nzw = mean(int.length.wald.nzw, na.rm = TRUE)
  mean.int.length.wm.nzw = mean(int.length.wm.nzw, na.rm = TRUE)
  mean.int.length.li.nzw = mean(int.length.li.nzw, na.rm = TRUE)
  mean.int.length.wald.c.nzw = mean(int.length.wald.c.nzw, na.rm = TRUE)
  mean.int.length.wt.c.nzw = mean(int.length.wt.c.nzw, na.rm = TRUE)
  mean.int.length.wn.c.nzw = mean(int.length.wn.c.nzw, na.rm = TRUE)
  mean.int.length.wz.c.nzw = mean(int.length.wz.c.nzw, na.rm = TRUE)
  mean.int.length.wm.c.nzw = mean(int.length.wm.c.nzw, na.rm = TRUE)
  mean.int.length.li.c.nzw = mean(int.length.li.c.nzw, na.rm = TRUE)
  
  #Compute the mean truncated and untruncated interval lengths for MCAR data,
  #including zero-width intervals
  mean.int.length.wt.mcar = mean(int.length.wt.mcar)
  mean.int.length.wz.mcar = mean(int.length.wz.mcar)
  mean.int.length.wn.mcar = mean(int.length.wn.mcar)
  mean.int.length.wald.mcar = mean(int.length.wald.mcar)
  mean.int.length.wm.mcar = mean(int.length.wm.mcar)
  mean.int.length.li.mcar = mean(int.length.li.mcar)
  mean.int.length.wald.c.mcar = mean(int.length.wald.c.mcar)
  mean.int.length.wt.c.mcar = mean(int.length.wt.c.mcar)
  mean.int.length.wn.c.mcar = mean(int.length.wn.c.mcar)
  mean.int.length.wz.c.mcar = mean(int.length.wz.c.mcar)
  mean.int.length.wm.c.mcar = mean(int.length.wm.c.mcar)
  mean.int.length.li.c.mcar = mean(int.length.li.c.mcar)
  
  #Compute the mean truncated and untruncated interval lengths for MCAR data,
  #leaving out zero-width intervals
  mean.int.length.wt.mcar.nzw = mean(int.length.wt.mcar.nzw, na.rm = TRUE)
  mean.int.length.wz.mcar.nzw = mean(int.length.wz.mcar.nzw, na.rm = TRUE)
  mean.int.length.wn.mcar.nzw = mean(int.length.wn.mcar.nzw, na.rm = TRUE)
  mean.int.length.wald.mcar.nzw = mean(int.length.wald.mcar.nzw, na.rm = TRUE)
  mean.int.length.wm.mcar.nzw = mean(int.length.wm.mcar.nzw, na.rm = TRUE)
  mean.int.length.li.mcar.nzw = mean(int.length.li.mcar.nzw, na.rm = TRUE)
  mean.int.length.wald.c.mcar.nzw = mean(int.length.wald.c.mcar.nzw, na.rm = TRUE)
  mean.int.length.wt.c.mcar.nzw = mean(int.length.wt.c.mcar.nzw, na.rm = TRUE)
  mean.int.length.wn.c.mcar.nzw = mean(int.length.wn.c.mcar.nzw, na.rm = TRUE)
  mean.int.length.wz.c.mcar.nzw = mean(int.length.wz.c.mcar.nzw, na.rm = TRUE)
  mean.int.length.wm.c.mcar.nzw = mean(int.length.wm.c.mcar.nzw, na.rm = TRUE)
  mean.int.length.li.c.mcar.nzw = mean(int.length.li.c.mcar.nzw, na.rm = TRUE)
  
  
  #Compute the proportion of intervals that go above 1 for MAR data
  pr.1.wn = length(ul.wn.v[ul.wn.v > 1])/length(ul.wn.v)
  pr.1.wm = length(ul.wm.v[ul.wm.v > 1])/length(ul.wm.v)
  pr.1.wt = length(ul.wt.v[ul.wt.v > 1])/length(ul.wt.v)
  pr.1.wz = length(ul.wz.v[ul.wz.v > 1])/length(ul.wz.v)
  pr.1.wald = length(ul.wald.v[ul.wald.v > 1])/length(ul.wald.v)
  pr.1.li = length(ul.li.v[ul.li.v > 1])/length(ul.li.v)
  
  #Compute the proportion of intervals that go above 1 for MCAR data
  pr.1.wn.mcar = length(ul.wn.v.mcar[ul.wn.v.mcar > 1])/length(ul.wn.v.mcar)
  pr.1.wm.mcar = length(ul.wm.v.mcar[ul.wm.v.mcar > 1])/length(ul.wm.v.mcar)
  pr.1.wt.mcar = length(ul.wt.v.mcar[ul.wt.v.mcar > 1])/length(ul.wt.v.mcar)
  pr.1.wz.mcar = length(ul.wz.v.mcar[ul.wz.v.mcar > 1])/length(ul.wz.v.mcar)
  pr.1.wald.mcar = length(ul.wald.v.mcar[ul.wald.v.mcar > 1])/length(ul.wald.v.mcar)
  pr.1.li.mcar = length(ul.li.v.mcar[ul.li.v.mcar > 1])/length(ul.li.v.mcar)
  
  #Compute the proportion of intervals that go below 0 for MAR data
  pr.0.wn = length(ll.wn.v[ll.wn.v < 0])/length(ll.wn.v)
  pr.0.wm = length(ll.wm.v[ll.wm.v < 0])/length(ll.wm.v)
  pr.0.wt = length(ll.wt.v[ll.wt.v < 0])/length(ll.wt.v)
  pr.0.wz = length(ll.wz.v[ll.wz.v < 0])/length(ll.wz.v)
  pr.0.wald = length(ll.wald.v[ll.wald.v < 0])/length(ll.wald.v)
  pr.0.li = length(ll.li.v[ll.li.v < 0])/length(ll.li.v)
  
  #Compute the proportion of intervals that go below 0 for MCAR data
  pr.0.wn.mcar = length(ll.wn.v.mcar[ll.wn.v.mcar < 0])/length(ll.wn.v.mcar)
  pr.0.wm.mcar = length(ll.wm.v.mcar[ll.wm.v.mcar < 0])/length(ll.wm.v.mcar)
  pr.0.wt.mcar = length(ll.wt.v.mcar[ll.wt.v.mcar < 0])/length(ll.wt.v.mcar)
  pr.0.wz.mcar = length(ll.wz.v.mcar[ll.wz.v.mcar < 0])/length(ll.wz.v.mcar)
  pr.0.wald.mcar = length(ll.wald.v.mcar[ll.wald.v.mcar < 0])/length(ll.wald.v.mcar)
  pr.0.li.mcar = length(ll.li.v.mcar[ll.li.v.mcar < 0])/length(ll.li.v.mcar)
  
  #Proportion of intervals that go below 0 for no missing data
  pr.0.s.wald = length(s.wald.ll[s.wald.ll < 0])/length(s.wald.ll)
  pr.0.wilson = length(wilson.ll[wilson.ll < 0])/length(wilson.ll)
  
  #Proportion of intervals that go above 1 for no missing data
  pr.1.s.wald = length(s.wald.ul[s.wald.ul > 1])/length(s.wald.ul)
  pr.1.wilson = length(wilson.ul[wilson.ul > 1])/length(wilson.ul)
  
  #Compute the coverage rate of intervals for MAR data
  cr.wn = length(p.in.int.wn[p.in.int.wn == TRUE])/length(p.in.int.wn)
  cr.wm = length(p.in.int.wm[p.in.int.wm == TRUE])/length(p.in.int.wm)
  cr.wt = length(p.in.int.wt[p.in.int.wt == TRUE])/length(p.in.int.wt)
  cr.wz = length(p.in.int.wz[p.in.int.wz == TRUE])/length(p.in.int.wz)
  cr.wald = length(p.in.int.wald[p.in.int.wald == TRUE])/length(p.in.int.wald)
  cr.li = length(p.in.int.li[p.in.int.li == TRUE])/length(p.in.int.li)
  
  #Compute the coverage rate of intervals for MCAR data
  cr.wn.mcar = length(p.in.int.wn.mcar[p.in.int.wn.mcar == TRUE])/length(p.in.int.wn.mcar)
  cr.wm.mcar = length(p.in.int.wm.mcar[p.in.int.wm.mcar == TRUE])/length(p.in.int.wm.mcar)
  cr.wt.mcar = length(p.in.int.wt.mcar[p.in.int.wt.mcar == TRUE])/length(p.in.int.wt.mcar)
  cr.wz.mcar = length(p.in.int.wz.mcar[p.in.int.wz.mcar == TRUE])/length(p.in.int.wz.mcar)
  cr.wald.mcar = (length(p.in.int.wald.mcar[p.in.int.wald.mcar == TRUE])/
                    length(p.in.int.wald.mcar))
  cr.li.mcar = length(p.in.int.li.mcar[p.in.int.li.mcar == TRUE])/length(p.in.int.li.mcar)
  
  #Calculate the proportion of zero-width intervals for MAR data
  pzw.wn = (length(nonzero.width.wn[nonzero.width.wn == FALSE])/
              length(nonzero.width.wn))
  pzw.wm = (length(nonzero.width.wm[nonzero.width.wm == FALSE])/
              length(nonzero.width.wm))
  pzw.wt = (length(nonzero.width.wt[nonzero.width.wt == FALSE])/
              length(nonzero.width.wt))
  pzw.wz = (length(nonzero.width.wz[nonzero.width.wz == FALSE])/
              length(nonzero.width.wz))
  pzw.wald = (length(nonzero.width.wald[nonzero.width.wald == FALSE])/
              length(nonzero.width.wald))
  pzw.li = (length(nonzero.width.li[nonzero.width.li == FALSE])/
              length(nonzero.width.li))
  
  #Compute the proportion of zero-width intervals for MCAR data
  pzw.wn.mcar = (length(nonzero.width.wn.mcar[nonzero.width.wn.mcar == FALSE])/
              length(nonzero.width.wn.mcar))
  pzw.wm.mcar = (length(nonzero.width.wm.mcar[nonzero.width.wm.mcar == FALSE])/
              length(nonzero.width.wm.mcar))
  pzw.wt.mcar = (length(nonzero.width.wt.mcar[nonzero.width.wt.mcar == FALSE])/
              length(nonzero.width.wt.mcar))
  pzw.wz.mcar = (length(nonzero.width.wz.mcar[nonzero.width.wz.mcar == FALSE])/
              length(nonzero.width.wz.mcar))
  pzw.wald.mcar = (length(nonzero.width.wald.mcar[nonzero.width.wald.mcar == FALSE])/
                length(nonzero.width.wald.mcar))
  pzw.li.mcar = (length(nonzero.width.li.mcar[nonzero.width.li.mcar == FALSE])/
                   length(nonzero.width.li.mcar))
  
  #Calculate the proportion of zero-width intervals for no missing data
  pzw.s.wald = (length(nonzero.width.s.wald[nonzero.width.s.wald == FALSE])/
                  length(nonzero.width.s.wald))
  pzw.wilson = (length(nonzero.width.wilson[nonzero.width.wilson == FALSE])/
                  length(nonzero.width.wilson))
  
  #Calculate the coverage rate of intervals for no missing data
  cr.wilson = (length(p.in.int.wilson[p.in.int.wilson == TRUE])/
                 length(p.in.int.wilson))
  cr.s.wald = (length(p.in.int.s.wald[p.in.int.s.wald == TRUE])/
                 length(p.in.int.s.wald))
  
  #Compute the mean interval lengths for no missing data, without truncation
  #and including zero-width intervals
  mean.int.length.s.wald = mean(int.length.s.wald)
  mean.int.length.wilson = mean(int.length.wilson)
  
  #Compute the mean interval lengths for no missing data with truncation,
  #including zero-width intervals
  mean.int.length.s.wald.c = mean(int.length.s.wald.c)
  mean.int.length.wilson.c = mean(int.length.wilson.c)
  
  #Compute the mean interval lengths for no missing data without zero-width intervals
  mean.int.length.s.wald.nzw = mean(int.length.s.wald.nzw, na.rm = TRUE)
  mean.int.length.wilson.nzw = mean(int.length.wilson.nzw, na.rm = TRUE)
  
  #Compute the mean truncated length of intervals without missing data,
  #leaving out zero-width intervals
  mean.int.length.s.wald.c.nzw = mean(int.length.s.wald.c.nzw, na.rm = TRUE)
  mean.int.length.wilson.c.nzw = mean(int.length.wilson.c.nzw, na.rm = TRUE)
  
  #Calculate the mean proportion of missingness
  mean.prop.miss = mean(prop.miss.v)
  
  #Calculate the proportion of times that bm, um, bm.mcar, and um.mcar are equal 
  #to 0
  pz.bm = length(bm.v[bm.v == 0])/length(bm.v)
  pz.um = length(um.v[um.v == 0])/length(um.v) 
  pz.bm.mcar = length(bm.v.mcar[bm.v.mcar == 0])/length(bm.v.mcar)
  pz.um.mcar = length(um.v.mcar[um.v.mcar == 0])/length(um.v.mcar)
  
  #Calculate the proportion of times that bm and um are both equal to 0, and the
  #proportion of times that bm.mcar and um.mcar are equal to 0
  pz.bm.and.um = length(pz.bm.um[pz.bm.um == TRUE])/length(pz.bm.um)
  pz.bm.and.um.mcar = length(pz.bm.um.mcar[pz.bm.um.mcar == TRUE])/length(pz.bm.um.mcar)

  
  
  return(data.frame(mean.int.length.wt = mean.int.length.wt,
                    mean.int.length.wz = mean.int.length.wz,
                    mean.int.length.wt.c = mean.int.length.wt.c,
                    mean.int.length.wz.c = mean.int.length.wz.c, 
                    mean.int.length.wald = mean.int.length.wald,
                    mean.int.length.wn = mean.int.length.wn,
                    mean.int.length.wald.c = mean.int.length.wald.c,
                    mean.int.length.wn.c = mean.int.length.wn.c,
                    mean.int.length.wt.mcar = mean.int.length.wt.mcar,
                    mean.int.length.wz.mcar = mean.int.length.wz.mcar,
                    mean.int.length.wt.c.mcar = mean.int.length.wt.c.mcar,
                    mean.int.length.wz.c.mcar = mean.int.length.wz.c.mcar, 
                    mean.int.length.wald.mcar = mean.int.length.wald.mcar,
                    mean.int.length.wn.mcar = mean.int.length.wn.mcar,
                    mean.int.length.wald.c.mcar = mean.int.length.wald.c.mcar,
                    mean.int.length.wn.c.mcar = mean.int.length.wn.c.mcar,
                    mean.int.length.wm = mean.int.length.wm,
                    mean.int.length.wm.c = mean.int.length.wm.c,
                    mean.int.length.wm.mcar = mean.int.length.wm.mcar,
                    mean.int.length.wm.c.mcar = mean.int.length.wm.c.mcar,
                    pr.1.wn = pr.1.wn,
                    pr.1.wm = pr.1.wm,
                    pr.1.wt = pr.1.wt,
                    pr.1.wz = pr.1.wz,
                    pr.1.wald = pr.1.wald,
                    pr.1.wn.mcar = pr.1.wn.mcar,
                    pr.1.wm.mcar = pr.1.wm.mcar,
                    pr.1.wt.mcar = pr.1.wt.mcar,
                    pr.1.wz.mcar = pr.1.wz.mcar,
                    pr.1.wald.mcar = pr.1.wald.mcar,
                    pr.0.wn = pr.0.wn,
                    pr.0.wm = pr.0.wm,
                    pr.0.wt = pr.0.wt,
                    pr.0.wz = pr.0.wz,
                    pr.0.wald = pr.0.wald,
                    pr.0.wn.mcar = pr.0.wn.mcar,
                    pr.0.wm.mcar = pr.0.wm.mcar,
                    pr.0.wt.mcar = pr.0.wt.mcar,
                    pr.0.wz.mcar = pr.0.wz.mcar,
                    pr.0.wald.mcar = pr.0.wald.mcar,
                    cr.wn = cr.wn,
                    cr.wm = cr.wm,
                    cr.wt = cr.wt,
                    cr.wz = cr.wz,
                    cr.wald = cr.wald,
                    cr.wn.mcar = cr.wn.mcar,
                    cr.wm.mcar = cr.wm.mcar,
                    cr.wt.mcar = cr.wt.mcar,
                    cr.wz.mcar = cr.wz.mcar,
                    cr.wald.mcar = cr.wald.mcar,
                    pzw.wn = pzw.wn,
                    pzw.wm = pzw.wm,
                    pzw.wt = pzw.wt,
                    pzw.wz = pzw.wz,
                    pzw.wald = pzw.wald,
                    pzw.wn.mcar = pzw.wn.mcar,
                    pzw.wm.mcar = pzw.wm.mcar,
                    pzw.wt.mcar = pzw.wt.mcar,
                    pzw.wz.mcar = pzw.wt.mcar,
                    pzw.wald.mcar = pzw.wald.mcar,
                    mean.prop.miss = mean.prop.miss,
                    mean.int.length.wt.nzw = mean.int.length.wt.nzw,
                    mean.int.length.wz.nzw = mean.int.length.wz.nzw,
                    mean.int.length.wn.nzw = mean.int.length.wn.nzw,
                    mean.int.length.wald.nzw = mean.int.length.wald.nzw,
                    mean.int.length.wm.nzw = mean.int.length.wm.nzw,
                    mean.int.length.wald.c.nzw = mean.int.length.wald.c.nzw,
                    mean.int.length.wt.c.nzw = mean.int.length.wt.c.nzw,
                    mean.int.length.wn.c.nzw = mean.int.length.wn.c.nzw,
                    mean.int.length.wz.c.nzw = mean.int.length.wz.c.nzw,
                    mean.int.length.wm.c.nzw = mean.int.length.wm.c.nzw,
                    mean.int.length.wt.mcar.nzw = mean.int.length.wt.mcar.nzw,
                    mean.int.length.wz.mcar.nzw = mean.int.length.wz.mcar.nzw,
                    mean.int.length.wn.mcar.nzw = mean.int.length.wn.mcar.nzw,
                    mean.int.length.wald.mcar.nzw = mean.int.length.wald.mcar.nzw,
                    mean.int.length.wm.mcar.nzw = mean.int.length.wm.mcar.nzw,
                    mean.int.length.wald.c.mcar.nzw = mean.int.length.wald.c.mcar.nzw,
                    mean.int.length.wt.c.mcar.nzw = mean.int.length.wt.c.mcar.nzw,
                    mean.int.length.wn.c.mcar.nzw = mean.int.length.wn.c.mcar.nzw,
                    mean.int.length.wz.c.mcar.nzw = mean.int.length.wz.c.mcar.nzw,
                    mean.int.length.wm.c.mcar.nzw = mean.int.length.wm.c.mcar.nzw,
                    pz.bm = pz.bm,
                    pz.um = pz.um,
                    pz.bm.mcar = pz.bm.mcar,
                    pz.um.mcar = pz.um.mcar,
                    pz.bm.and.um = pz.bm.and.um,
                    pz.bm.and.um.mcar = pz.bm.and.um.mcar,
                    pzw.s.wald = pzw.s.wald,
                    pzw.wilson = pzw.wilson,
                    cr.wilson = cr.wilson,
                    cr.s.wald = cr.s.wald,
                    mean.int.length.s.wald = mean.int.length.s.wald,
                    mean.int.length.wilson = mean.int.length.wilson,
                    mean.int.length.s.wald.c = mean.int.length.s.wald.c,
                    mean.int.length.wilson.c= mean.int.length.wilson.c,
                    mean.int.length.s.wald.nzw = mean.int.length.s.wald.nzw,
                    mean.int.length.wilson.nzw = mean.int.length.wilson.nzw,
                    mean.int.length.s.wald.c.nzw = mean.int.length.s.wald.c.nzw,
                    mean.int.length.wilson.c.nzw = mean.int.length.wilson.c.nzw,
                    pr.0.s.wald = pr.0.s.wald,
                    pr.0.wilson = pr.0.wilson,
                    pr.1.s.wald = pr.1.s.wald,
                    pr.1.wilson = pr.1.wilson,
                    mean.int.length.li = mean.int.length.li,
                    mean.int.length.li.c = mean.int.length.li.c,
                    mean.int.length.li.mcar = mean.int.length.li.mcar,
                    mean.int.length.li.c.mcar = mean.int.length.li.c.mcar,
                    pr.1.li = pr.1.li,
                    pr.1.li.mcar = pr.1.li.mcar,
                    pr.0.li = pr.0.li,
                    pr.0.li.mcar = pr.0.li.mcar,
                    cr.li = cr.li,
                    cr.li.mcar = cr.li.mcar,
                    pzw.li = pzw.li,
                    pzw.li.mcar = pzw.li.mcar,
                    mean.int.length.li.nzw = mean.int.length.li.nzw,
                    mean.int.length.li.c.nzw = mean.int.length.li.c.nzw,
                    mean.int.length.li.mcar.nzw = mean.int.length.li.mcar.nzw,
                    mean.int.length.li.c.mcar.nzw = mean.int.length.li.c.mcar.nzw))
}
#Applying function 100 times
x = rep(0.05, 100) #change the input to the desired value of p
y = sapply(x, mult.imp, n = 500, a = 1, b = 1, it = 1000, m = 5, mp = 0.3)

wt.mean.int.length = c()
wz.mean.int.length = c()
wt.mean.int.length.c = c()
wz.mean.int.length.c = c()
wald.mean.int.length = c()
wn.mean.int.length = c()
wald.mean.int.length.c = c()
wn.mean.int.length.c = c()
wt.mean.int.length.mcar = c()
wz.mean.int.length.mcar = c()
wt.mean.int.length.c.mcar = c()
wz.mean.int.length.c.mcar = c()
wald.mean.int.length.mcar = c()
wn.mean.int.length.mcar = c()
wald.mean.int.length.c.mcar = c()
wn.mean.int.length.c.mcar = c()
wm.mean.int.length = c()
wm.mean.int.length.c = c()
wm.mean.int.length.mcar = c()
wm.mean.int.length.c.mcar = c()
wn.prop.above.1 = c()
wm.prop.above.1 = c()
wt.prop.above.1 = c()
wz.prop.above.1 = c()
wald.prop.above.1 = c()
wn.prop.above.1.mcar = c()
wm.prop.above.1.mcar = c()
wt.prop.above.1.mcar = c()
wz.prop.above.1.mcar = c()
wald.prop.above.1.mcar = c()
wn.prop.below.0 = c()
wm.prop.below.0 = c()
wt.prop.below.0 = c()
wz.prop.below.0 = c()
wald.prop.below.0 = c()
wn.prop.below.0.mcar = c()
wm.prop.below.0.mcar = c()
wt.prop.below.0.mcar = c()
wz.prop.below.0.mcar = c()
wald.prop.below.0.mcar = c()
wn.cov.rate = c()
wm.cov.rate = c()
wt.cov.rate = c()
wz.cov.rate = c()
wald.cov.rate = c()
wn.cov.rate.mcar = c()
wm.cov.rate.mcar = c()
wt.cov.rate.mcar = c()
wz.cov.rate.mcar = c()
wald.cov.rate.mcar = c()
wn.prop.zero.width = c()
wm.prop.zero.width = c()
wt.prop.zero.width = c()
wz.prop.zero.width = c()
wald.prop.zero.width = c()
wn.prop.zero.width.mcar = c()
wm.prop.zero.width.mcar = c()
wt.prop.zero.width.mcar = c()
wz.prop.zero.width.mcar = c()
wald.prop.zero.width.mcar = c()
mean.prop.miss = c()
wt.mean.int.length.nzw = c()
wz.mean.int.length.nzw = c()
wn.mean.int.length.nzw = c()
wald.mean.int.length.nzw = c()
wm.mean.int.length.nzw = c()
wald.mean.int.length.c.nzw = c()
wt.mean.int.length.c.nzw = c()
wn.mean.int.length.c.nzw = c()
wz.mean.int.length.c.nzw = c()
wm.mean.int.length.c.nzw = c()
wt.mean.int.length.mcar.nzw = c()
wz.mean.int.length.mcar.nzw = c()
wn.mean.int.length.mcar.nzw = c()
wald.mean.int.length.mcar.nzw = c()
wm.mean.int.length.mcar.nzw = c()
wald.mean.int.length.c.mcar.nzw = c()
wt.mean.int.length.c.mcar.nzw = c()
wn.mean.int.length.c.mcar.nzw = c()
wz.mean.int.length.c.mcar.nzw = c()
wm.mean.int.length.c.mcar.nzw = c()
prop.zero.bm = c()
prop.zero.um = c()
prop.zero.bm.mcar = c()
prop.zero.um.mcar = c()
prop.zero.bm.um = c()
prop.zero.bm.um.mcar = c()
s.wald.prop.zero.width = c()
wilson.prop.zero.width = c()
wilson.cov.rate = c()
s.wald.cov.rate = c()
s.wald.mean.int.length = c()
wilson.mean.int.length = c()
s.wald.mean.int.length.c = c()
wilson.mean.int.length.c = c()
s.wald.mean.int.length.nzw = c()
wilson.mean.int.length.nzw = c()
s.wald.mean.int.length.c.nzw = c()
wilson.mean.int.length.c.nzw = c()
s.wald.prop.below.0 = c()
wilson.prop.below.0 = c()
s.wald.prop.above.1 = c()
wilson.prop.above.1 = c()
li.mean.int.length = c()
li.mean.int.length.c = c()
li.mean.int.length.mcar = c()
li.mean.int.length.c.mcar = c()
li.prop.above.1 = c()
li.prop.above.1.mcar = c()
li.prop.below.0 = c()
li.prop.below.0.mcar = c()
li.cov.rate = c()
li.cov.rate.mcar = c()
li.prop.zero.width = c()
li.prop.zero.width.mcar = c()
li.mean.int.length.nzw = c()
li.mean.int.length.c.nzw = c()
li.mean.int.length.mcar.nzw = c()
li.mean.int.length.c.mcar.nzw = c()

for (k3 in 1:ncol(y)) {
  wt.mean.int.length = c(wt.mean.int.length, y[[1, k3]])
  wz.mean.int.length = c(wz.mean.int.length, y[[2, k3]])
  wt.mean.int.length.c = c(wt.mean.int.length.c, y[[3, k3]])
  wz.mean.int.length.c = c(wz.mean.int.length.c, y[[4, k3]])
  wald.mean.int.length = c(wald.mean.int.length, y[[5, k3]])
  wn.mean.int.length = c(wn.mean.int.length, y[[6, k3]])
  wald.mean.int.length.c = c(wald.mean.int.length.c, y[[7, k3]])
  wn.mean.int.length.c = c(wn.mean.int.length.c, y[[8, k3]])
  wt.mean.int.length.mcar = c(wt.mean.int.length.mcar, y[[9, k3]])
  wz.mean.int.length.mcar = c(wz.mean.int.length.mcar, y[[10, k3]])
  wt.mean.int.length.c.mcar = c(wt.mean.int.length.c.mcar, y[[11, k3]])
  wz.mean.int.length.c.mcar = c(wz.mean.int.length.c.mcar, y[[12, k3]])
  wald.mean.int.length.mcar = c(wald.mean.int.length.mcar, y[[13, k3]])
  wn.mean.int.length.mcar = c(wn.mean.int.length.mcar, y[[14, k3]])
  wald.mean.int.length.c.mcar = c(wald.mean.int.length.c.mcar, y[[15, k3]])
  wn.mean.int.length.c.mcar = c(wn.mean.int.length.c.mcar, y[[16, k3]])
  wm.mean.int.length = c(wm.mean.int.length, y[[17, k3]])
  wm.mean.int.length.c = c(wm.mean.int.length.c, y[[18, k3]])
  wm.mean.int.length.mcar = c(wm.mean.int.length.mcar, y[[19, k3]])
  wm.mean.int.length.c.mcar = c(wm.mean.int.length.c.mcar, y[[20, k3]])
  wn.prop.above.1 = c(wn.prop.above.1, y[[21, k3]])
  wm.prop.above.1 = c(wm.prop.above.1, y[[22, k3]])
  wt.prop.above.1 = c(wt.prop.above.1, y[[23, k3]])
  wz.prop.above.1 = c(wz.prop.above.1, y[[24, k3]])
  wald.prop.above.1 = c(wald.prop.above.1, y[[25, k3]])
  wn.prop.above.1.mcar = c(wn.prop.above.1.mcar, y[[26, k3]])
  wm.prop.above.1.mcar = c(wm.prop.above.1.mcar, y[[27, k3]])
  wt.prop.above.1.mcar = c(wt.prop.above.1.mcar, y[[28, k3]])
  wz.prop.above.1.mcar = c(wz.prop.above.1.mcar, y[[29, k3]])
  wald.prop.above.1.mcar = c(wald.prop.above.1.mcar, y[[30, k3]])
  wn.prop.below.0 = c(wn.prop.below.0, y[[31, k3]])
  wm.prop.below.0 = c(wm.prop.below.0, y[[32, k3]])
  wt.prop.below.0 = c(wt.prop.below.0, y[[33, k3]])
  wz.prop.below.0 = c(wz.prop.below.0, y[[34, k3]])
  wald.prop.below.0 = c(wald.prop.below.0, y[[35, k3]])
  wn.prop.below.0.mcar = c(wn.prop.below.0.mcar, y[[36, k3]])
  wm.prop.below.0.mcar = c(wm.prop.below.0.mcar, y[[37, k3]])
  wt.prop.below.0.mcar = c(wt.prop.below.0.mcar, y[[38, k3]])
  wz.prop.below.0.mcar = c(wz.prop.below.0.mcar, y[[39, k3]])
  wald.prop.below.0.mcar = c(wald.prop.below.0.mcar, y[[40, k3]])
  wn.cov.rate = c(wn.cov.rate, y[[41, k3]])
  wm.cov.rate = c(wm.cov.rate, y[[42, k3]])
  wt.cov.rate = c(wt.cov.rate, y[[43, k3]])
  wz.cov.rate = c(wz.cov.rate, y[[44, k3]])
  wald.cov.rate = c(wald.cov.rate, y[[45, k3]])
  wn.cov.rate.mcar = c(wn.cov.rate.mcar, y[[46, k3]])
  wm.cov.rate.mcar = c(wm.cov.rate.mcar, y[[47, k3]])
  wt.cov.rate.mcar = c(wt.cov.rate.mcar, y[[48, k3]])
  wz.cov.rate.mcar = c(wz.cov.rate.mcar, y[[49, k3]])
  wald.cov.rate.mcar = c(wald.cov.rate.mcar, y[[50, k3]])
  wn.prop.zero.width = c(wn.prop.zero.width, y[[51, k3]])
  wm.prop.zero.width = c(wm.prop.zero.width, y[[52, k3]])
  wt.prop.zero.width = c(wt.prop.zero.width, y[[53, k3]])
  wz.prop.zero.width = c(wz.prop.zero.width, y[[54, k3]])
  wald.prop.zero.width = c(wald.prop.zero.width, y[[55, k3]])
  wn.prop.zero.width.mcar = c(wn.prop.zero.width.mcar, y[[56, k3]])
  wm.prop.zero.width.mcar = c(wm.prop.zero.width.mcar, y[[57, k3]])
  wt.prop.zero.width.mcar = c(wt.prop.zero.width.mcar, y[[58, k3]])
  wz.prop.zero.width.mcar = c(wz.prop.zero.width.mcar, y[[59, k3]])
  wald.prop.zero.width.mcar = c(wald.prop.zero.width.mcar, y[[60, k3]])
  mean.prop.miss = c(mean.prop.miss, y[[61, k3]])
  wt.mean.int.length.nzw = c(wt.mean.int.length.nzw, y[[62, k3]])
  wz.mean.int.length.nzw = c(wz.mean.int.length.nzw, y[[63, k3]])
  wn.mean.int.length.nzw = c(wn.mean.int.length.nzw, y[[64, k3]])
  wald.mean.int.length.nzw = c(wald.mean.int.length.nzw, y[[65, k3]])
  wm.mean.int.length.nzw = c(wm.mean.int.length.nzw, y[[66, k3]])
  wald.mean.int.length.c.nzw = c(wald.mean.int.length.c.nzw, y[[67, k3]])
  wt.mean.int.length.c.nzw = c(wt.mean.int.length.c.nzw, y[[68, k3]])
  wn.mean.int.length.c.nzw = c(wn.mean.int.length.c.nzw, y[[69, k3]])
  wz.mean.int.length.c.nzw = c(wz.mean.int.length.c.nzw, y[[70, k3]])
  wm.mean.int.length.c.nzw = c(wm.mean.int.length.c.nzw, y[[71, k3]])
  wt.mean.int.length.mcar.nzw = c(wt.mean.int.length.mcar.nzw, y[[72, k3]])
  wz.mean.int.length.mcar.nzw = c(wz.mean.int.length.mcar.nzw, y[[73, k3]])
  wn.mean.int.length.mcar.nzw = c(wn.mean.int.length.mcar.nzw, y[[74, k3]])
  wald.mean.int.length.mcar.nzw = c(wald.mean.int.length.mcar.nzw, y[[75, k3]])
  wm.mean.int.length.mcar.nzw = c(wm.mean.int.length.mcar.nzw, y[[76, k3]])
  wald.mean.int.length.c.mcar.nzw = c(wald.mean.int.length.c.mcar.nzw, y[[77, k3]])
  wt.mean.int.length.c.mcar.nzw = c(wt.mean.int.length.c.mcar.nzw, y[[78, k3]])
  wn.mean.int.length.c.mcar.nzw = c(wn.mean.int.length.c.mcar.nzw, y[[79, k3]])
  wz.mean.int.length.c.mcar.nzw = c(wz.mean.int.length.c.mcar.nzw, y[[80, k3]])
  wm.mean.int.length.c.mcar.nzw = c(wm.mean.int.length.c.mcar.nzw, y[[81, k3]])
  prop.zero.bm = c(prop.zero.bm, y[[82, k3]])
  prop.zero.um = c(prop.zero.um, y[[83, k3]])
  prop.zero.bm.mcar = c(prop.zero.bm.mcar, y[[84, k3]])
  prop.zero.um.mcar = c(prop.zero.um.mcar, y[[85, k3]])
  prop.zero.bm.um = c(prop.zero.bm.um, y[[86, k3]])
  prop.zero.bm.um.mcar = c(prop.zero.bm.um.mcar, y[[87, k3]])
  s.wald.prop.zero.width = c(s.wald.prop.zero.width, y[[88, k3]])
  wilson.prop.zero.width = c(wilson.prop.zero.width, y[[89, k3]])
  wilson.cov.rate = c(wilson.cov.rate, y[[90, k3]])
  s.wald.cov.rate = c(s.wald.cov.rate, y[[91, k3]])
  s.wald.mean.int.length = c(s.wald.mean.int.length, y[[92, k3]])
  wilson.mean.int.length = c(wilson.mean.int.length, y[[93, k3]])
  s.wald.mean.int.length.c = c(s.wald.mean.int.length.c, y[[94, k3]])
  wilson.mean.int.length.c = c(wilson.mean.int.length.c, y[[95, k3]])
  s.wald.mean.int.length.nzw = c(s.wald.mean.int.length.nzw, y[[96, k3]])
  wilson.mean.int.length.nzw = c(wilson.mean.int.length.nzw, y[[97, k3]])
  s.wald.mean.int.length.c.nzw = c(s.wald.mean.int.length.c.nzw, y[[98, k3]])
  wilson.mean.int.length.c.nzw = c(wilson.mean.int.length.c.nzw, y[[99, k3]])
  s.wald.prop.below.0 = c(s.wald.prop.below.0, y[[100, k3]])
  wilson.prop.below.0 = c(wilson.prop.below.0, y[[101, k3]])
  s.wald.prop.above.1 = c(s.wald.prop.above.1, y[[102, k3]])
  wilson.prop.above.1 = c(wilson.prop.above.1, y[[103, k3]])
  li.mean.int.length = c(li.mean.int.length, y[[104, k3]])
  li.mean.int.length.c = c(li.mean.int.length.c, y[[105, k3]])
  li.mean.int.length.mcar = c(li.mean.int.length.mcar, y[[106, k3]])
  li.mean.int.length.c.mcar = c(li.mean.int.length.c.mcar, y[[107, k3]])
  li.prop.above.1 = c(li.prop.above.1, y[[108, k3]])
  li.prop.above.1.mcar = c(li.prop.above.1.mcar, y[[109, k3]])
  li.prop.below.0 = c(li.prop.below.0, y[[110, k3]])
  li.prop.below.0.mcar = c(li.prop.below.0.mcar, y[[111, k3]])
  li.cov.rate = c(li.cov.rate, y[[112, k3]])
  li.cov.rate.mcar = c(li.cov.rate.mcar, y[[113, k3]])
  li.prop.zero.width = c(li.prop.zero.width, y[[114, k3]])
  li.prop.zero.width.mcar = c(li.prop.zero.width.mcar, y[[115, k3]])
  li.mean.int.length.nzw = c(li.mean.int.length.nzw, y[[116, k3]])
  li.mean.int.length.c.nzw = c(li.mean.int.length.c.nzw, y[[117, k3]])
  li.mean.int.length.mcar.nzw = c(li.mean.int.length.mcar.nzw, y[[118, k3]])
  li.mean.int.length.c.mcar.nzw = c(li.mean.int.length.c.mcar.nzw, y[[119, k3]])
}

#Calculate the mean and variance of each returned value of the mult.imp function
mean.wt.mean.int.length = mean(wt.mean.int.length)
var.wt.mean.int.length = var(wt.mean.int.length)

mean.wz.mean.int.length = mean(wz.mean.int.length)
var.wz.mean.int.length = var(wz.mean.int.length)

mean.wt.mean.int.length.c = mean(wt.mean.int.length.c)
var.wt.mean.int.length.c = var(wt.mean.int.length.c)

mean.wz.mean.int.length.c = mean(wz.mean.int.length.c)
var.wz.mean.int.length.c = var(wz.mean.int.length.c)

mean.wald.mean.int.length = mean(wald.mean.int.length)
var.wald.mean.int.length = var(wald.mean.int.length)

mean.wald.mean.int.length.c = mean(wald.mean.int.length.c)
var.wald.mean.int.length.c = var(wald.mean.int.length.c)

mean.wn.mean.int.length = mean(wn.mean.int.length)
var.wn.mean.int.length = var(wn.mean.int.length)

mean.wn.mean.int.length.c = mean(wn.mean.int.length.c)
var.wn.mean.int.length.c = var(wn.mean.int.length.c)

mean.wt.mean.int.length.mcar = mean(wt.mean.int.length.mcar)
var.wt.mean.int.length.mcar = var(wt.mean.int.length.mcar)

mean.wz.mean.int.length.mcar = mean(wz.mean.int.length.mcar)
var.wz.mean.int.length.mcar = var(wz.mean.int.length.mcar)

mean.wt.mean.int.length.c.mcar = mean(wt.mean.int.length.c.mcar)
var.wt.mean.int.length.c.mcar = var(wt.mean.int.length.c.mcar)

mean.wz.mean.int.length.c.mcar = mean(wz.mean.int.length.c.mcar)
var.wz.mean.int.length.c.mcar = var(wz.mean.int.length.c.mcar)

mean.wald.mean.int.length.mcar = mean(wald.mean.int.length.mcar)
var.wald.mean.int.length.mcar = var(wald.mean.int.length.mcar)

mean.wn.mean.int.length.mcar = mean(wn.mean.int.length.mcar)
var.wn.mean.int.length.mcar = var(wn.mean.int.length.mcar)

mean.wald.mean.int.length.c.mcar = mean(wald.mean.int.length.c.mcar)
var.wald.mean.int.length.c.mcar = var(wald.mean.int.length.c.mcar)

mean.wn.mean.int.length.c.mcar = mean(wn.mean.int.length.c.mcar)
var.wn.mean.int.length.c.mcar = var(wn.mean.int.length.c.mcar)

mean.wm.mean.int.length = mean(wm.mean.int.length)
var.wm.mean.int.length = var(wm.mean.int.length)

mean.wm.mean.int.length.c = mean(wm.mean.int.length.c)
var.wm.mean.int.length.c = var(wm.mean.int.length.c)

mean.wm.mean.int.length.mcar = mean(wm.mean.int.length.mcar)
var.wm.mean.int.length.mcar = var(wm.mean.int.length.mcar)

mean.wm.mean.int.length.c.mcar = mean(wm.mean.int.length.c.mcar)
var.wm.mean.int.length.c.mcar = var(wm.mean.int.length.c.mcar)

mean.wn.prop.above.1 = mean(wn.prop.above.1)
var.wn.prop.above.1 = var(wn.prop.above.1)

mean.wm.prop.above.1 = mean(wm.prop.above.1)
var.wm.prop.above.1 = var(wm.prop.above.1)

mean.wt.prop.above.1 = mean(wt.prop.above.1)
var.wt.prop.above.1 = var(wt.prop.above.1)

mean.wz.prop.above.1 = mean(wz.prop.above.1)
var.wz.prop.above.1 = var(wz.prop.above.1)

mean.wald.prop.above.1 = mean(wald.prop.above.1)
var.wald.prop.above.1 = var(wald.prop.above.1)

mean.wn.prop.above.1.mcar = mean(wn.prop.above.1.mcar)
var.wn.prop.above.1.mcar = var(wn.prop.above.1.mcar)

mean.wm.prop.above.1.mcar = mean(wm.prop.above.1.mcar)
var.wm.prop.above.1.mcar = var(wm.prop.above.1.mcar)

mean.wt.prop.above.1.mcar = mean(wt.prop.above.1.mcar)
var.wt.prop.above.1.mcar = var(wt.prop.above.1.mcar)

mean.wz.prop.above.1.mcar = mean(wz.prop.above.1.mcar)
var.wz.prop.above.1.mcar = var(wz.prop.above.1.mcar)

mean.wald.prop.above.1.mcar = mean(wald.prop.above.1.mcar)
var.wald.prop.above.1.mcar = var(wald.prop.above.1.mcar)

mean.wn.prop.below.0 = mean(wn.prop.below.0)
var.wn.prop.below.0 = var(wn.prop.below.0)

mean.wm.prop.below.0 = mean(wm.prop.below.0)
var.wm.prop.below.0 = var(wm.prop.below.0)

mean.wt.prop.below.0 = mean(wt.prop.below.0)
var.wt.prop.below.0 = var(wt.prop.below.0)

mean.wz.prop.below.0 = mean(wz.prop.below.0)
var.wz.prop.below.0 = var(wz.prop.below.0)

mean.wald.prop.below.0 = mean(wald.prop.below.0)
var.wald.prop.below.0 = var(wald.prop.below.0)

mean.wn.prop.below.0.mcar = mean(wn.prop.below.0.mcar)
var.wn.prop.below.0.mcar = var(wn.prop.below.0.mcar)

mean.wm.prop.below.0.mcar = mean(wm.prop.below.0.mcar)
var.wm.prop.below.0.mcar = var(wm.prop.below.0.mcar)

mean.wt.prop.below.0.mcar = mean(wt.prop.below.0.mcar)
var.wt.prop.below.0.mcar = var(wt.prop.below.0.mcar)

mean.wz.prop.below.0.mcar = mean(wz.prop.below.0.mcar)
var.wz.prop.below.0.mcar = var(wz.prop.below.0.mcar)

mean.wald.prop.below.0.mcar = mean(wald.prop.below.0.mcar)
var.wald.prop.below.0.mcar = var(wald.prop.below.0.mcar)

mean.wn.cov.rate = mean(wn.cov.rate)
var.wn.cov.rate = var(wn.cov.rate)

mean.wm.cov.rate = mean(wm.cov.rate)
var.wm.cov.rate = var(wm.cov.rate)

mean.wt.cov.rate = mean(wt.cov.rate)
var.wt.cov.rate = var(wt.cov.rate)

mean.wz.cov.rate =  mean(wz.cov.rate)
var.wz.cov.rate = var(wz.cov.rate)

mean.wald.cov.rate = mean(wald.cov.rate)
var.wald.cov.rate = var(wald.cov.rate)

mean.wn.cov.rate.mcar = mean(wn.cov.rate.mcar)
var.wn.cov.rate.mcar = var(wn.cov.rate.mcar)

mean.wm.cov.rate.mcar = mean(wm.cov.rate.mcar)
var.wm.cov.rate.mcar = var(wm.cov.rate.mcar)

mean.wt.cov.rate.mcar = mean(wt.cov.rate.mcar)
var.wt.cov.rate.mcar = var(wt.cov.rate.mcar)

mean.wz.cov.rate.mcar = mean(wz.cov.rate.mcar)
var.wz.cov.rate.mcar = var(wz.cov.rate.mcar)

mean.wald.cov.rate.mcar = mean(wald.cov.rate.mcar)
var.wald.cov.rate.mcar = var(wald.cov.rate.mcar)

mean.wn.prop.zero.width = mean(wn.prop.zero.width)
var.wn.prop.zero.width = var(wn.prop.zero.width)

mean.wm.prop.zero.width = mean(wm.prop.zero.width)
var.wm.prop.zero.width = var(wm.prop.zero.width)

mean.wt.prop.zero.width = mean(wt.prop.zero.width)
var.wt.prop.zero.width = var(wt.prop.zero.width)

mean.wz.prop.zero.width = mean(wz.prop.zero.width)
var.wz.prop.zero.width = var(wz.prop.zero.width)

mean.wald.prop.zero.width = mean(wald.prop.zero.width)
var.wald.prop.zero.width = var(wald.prop.zero.width)

mean.wn.prop.zero.width.mcar = mean(wn.prop.zero.width.mcar)
var.wn.prop.zero.width.mcar = var(wn.prop.zero.width.mcar)

mean.wm.prop.zero.width.mcar = mean(wm.prop.zero.width.mcar)
var.wm.prop.zero.width.mcar = var(wm.prop.zero.width.mcar)

mean.wt.prop.zero.width.mcar = mean(wt.prop.zero.width.mcar)
var.wt.prop.zero.width.mcar = var(wt.prop.zero.width.mcar)

mean.wz.prop.zero.width.mcar = mean(wz.prop.zero.width.mcar)
var.wz.prop.zero.width.mcar = var(wz.prop.zero.width.mcar)

mean.wald.prop.zero.width.mcar = mean(wald.prop.zero.width.mcar)
var.wald.prop.zero.width.mcar = var(wald.prop.zero.width.mcar)

mean.mean.prop.miss = mean(mean.prop.miss)
var.mean.prop.miss = var(mean.prop.miss)

mean.wt.mean.int.length.nzw = mean(wt.mean.int.length.nzw)
var.wt.mean.int.length.nzw = var(wt.mean.int.length.nzw)

mean.wz.mean.int.length.nzw = mean(wz.mean.int.length.nzw)
var.wz.mean.int.length.nzw = var(wz.mean.int.length.nzw)

mean.wn.mean.int.length.nzw = mean(wn.mean.int.length.nzw)
var.wn.mean.int.length.nzw = var(wn.mean.int.length.nzw)

mean.wald.mean.int.length.nzw = mean(wald.mean.int.length.nzw)
var.wald.mean.int.length.nzw = var(wald.mean.int.length.nzw)

mean.wm.mean.int.length.nzw = mean(wm.mean.int.length.nzw)
var.wm.mean.int.length.nzw = var(wm.mean.int.length.nzw)

mean.wald.mean.int.length.c.nzw = mean(wald.mean.int.length.c.nzw)
var.wald.mean.int.length.c.nzw = var(wald.mean.int.length.c.nzw)

mean.wt.mean.int.length.c.nzw = mean(wt.mean.int.length.c.nzw)
var.wt.mean.int.length.c.nzw = var(wt.mean.int.length.c.nzw)

mean.wn.mean.int.length.c.nzw = mean(wn.mean.int.length.c.nzw)
var.wn.mean.int.length.c.nzw = var(wn.mean.int.length.c.nzw)

mean.wz.mean.int.length.c.nzw = mean(wz.mean.int.length.c.nzw)
var.wz.mean.int.length.c.nzw = var(wz.mean.int.length.c.nzw)

mean.wm.mean.int.length.c.nzw = mean(wm.mean.int.length.c.nzw)
var.wm.mean.int.length.c.nzw = var(wm.mean.int.length.c.nzw)

mean.wt.mean.int.length.mcar.nzw = mean(wt.mean.int.length.mcar.nzw)
var.wt.mean.int.length.mcar.nzw = var(wt.mean.int.length.mcar.nzw)

mean.wz.mean.int.length.mcar.nzw = mean(wz.mean.int.length.mcar.nzw)
var.wz.mean.int.length.mcar.nzw = var(wz.mean.int.length.mcar.nzw)

mean.wn.mean.int.length.mcar.nzw = mean(wn.mean.int.length.mcar.nzw)
var.wn.mean.int.length.mcar.nzw = var(wn.mean.int.length.mcar.nzw)

mean.wald.mean.int.length.mcar.nzw = mean(wald.mean.int.length.mcar.nzw)
var.wald.mean.int.length.mcar.nzw = var(wald.mean.int.length.mcar.nzw)

mean.wm.mean.int.length.mcar.nzw = mean(wm.mean.int.length.mcar.nzw)
var.wm.mean.int.length.mcar.nzw = var(wm.mean.int.length.mcar.nzw)

mean.wald.mean.int.length.c.mcar.nzw = mean(wald.mean.int.length.c.mcar.nzw)
var.wald.mean.int.length.c.mcar.nzw = var(wald.mean.int.length.c.mcar.nzw)

mean.wt.mean.int.length.c.mcar.nzw = mean(wt.mean.int.length.c.mcar.nzw)
var.wt.mean.int.length.c.mcar.nzw = var(wt.mean.int.length.c.mcar.nzw)

mean.wn.mean.int.length.c.mcar.nzw = mean(wn.mean.int.length.c.mcar.nzw)
var.wn.mean.int.length.c.mcar.nzw = var(wn.mean.int.length.c.mcar.nzw)

mean.wz.mean.int.length.c.mcar.nzw = mean(wz.mean.int.length.c.mcar.nzw)
var.wz.mean.int.length.c.mcar.nzw = var(wz.mean.int.length.c.mcar.nzw)

mean.wm.mean.int.length.c.mcar.nzw = mean(wm.mean.int.length.c.mcar.nzw)
var.wm.mean.int.length.c.mcar.nzw = var(wm.mean.int.length.c.mcar.nzw)

mean.prop.zero.bm = mean(prop.zero.bm)
var.prop.zero.bm = var(prop.zero.bm)

mean.prop.zero.um = mean(prop.zero.um)
var.prop.zero.um = var(prop.zero.um)

mean.prop.zero.bm.mcar = mean(prop.zero.bm.mcar)
var.prop.zero.bm.mcar = var(prop.zero.bm.mcar)

mean.prop.zero.um.mcar = mean(prop.zero.um.mcar)
var.prop.zero.um.mcar = var(prop.zero.um.mcar)

mean.prop.zero.bm.um = mean(prop.zero.bm.um)
var.prop.zero.bm.um = var(prop.zero.bm.um)

mean.prop.zero.bm.um.mcar = mean(prop.zero.bm.um.mcar)
var.prop.zero.bm.um.mcar = var(prop.zero.bm.um.mcar)

mean.s.wald.prop.zero.width = mean(s.wald.prop.zero.width)
var.s.wald.prop.zero.width = var(s.wald.prop.zero.width)

mean.wilson.prop.zero.width = mean(wilson.prop.zero.width)
var.wilson.prop.zero.width = var(wilson.prop.zero.width)

mean.s.wald.cov.rate = mean(s.wald.cov.rate)
var.s.wald.cov.rate = var(s.wald.cov.rate)

mean.wilson.cov.rate = mean(wilson.cov.rate)
var.wilson.cov.rate = var(wilson.cov.rate)

mean.s.wald.mean.int.length = mean(s.wald.mean.int.length)
var.s.wald.mean.int.length = var(s.wald.mean.int.length)

mean.wilson.mean.int.length = mean(wilson.mean.int.length)
var.wilson.mean.int.length = var(wilson.mean.int.length)

mean.s.wald.mean.int.length.c = mean(s.wald.mean.int.length.c)
var.s.wald.mean.int.length.c = var(s.wald.mean.int.length.c)

mean.wilson.mean.int.length.c = mean(wilson.mean.int.length.c)
var.wilson.mean.int.length.c = var(wilson.mean.int.length.c)

mean.s.wald.mean.int.length.nzw = mean(s.wald.mean.int.length.nzw)
var.s.wald.mean.int.length.nzw = var(s.wald.mean.int.length.nzw)

mean.wilson.mean.int.length.nzw = mean(wilson.mean.int.length.nzw)
var.wilson.mean.int.length.nzw = var(wilson.mean.int.length.nzw)

mean.s.wald.mean.int.length.c.nzw = mean(s.wald.mean.int.length.c.nzw)
var.s.wald.mean.int.length.c.nzw = var(s.wald.mean.int.length.c.nzw)

mean.wilson.mean.int.length.c.nzw = mean(wilson.mean.int.length.c.nzw)
var.wilson.mean.int.length.c.nzw = var(wilson.mean.int.length.c.nzw)

mean.s.wald.prop.below.0 = mean(s.wald.prop.below.0)
var.s.wald.prop.below.0 = var(s.wald.prop.below.0)

mean.wilson.prop.below.0 = mean(wilson.prop.below.0)
var.wilson.prop.below.0 = var(wilson.prop.below.0)

mean.s.wald.prop.above.1 = mean(s.wald.prop.above.1)
var.s.wald.prop.above.1 = var(s.wald.prop.above.1)

mean.wilson.prop.above.1 = mean(wilson.prop.above.1)
var.wilson.prop.above.1 = var(wilson.prop.above.1)

mean.li.mean.int.length = mean(li.mean.int.length)
var.li.mean.int.length = var(li.mean.int.length)

mean.li.mean.int.length.c = mean(li.mean.int.length.c)
var.li.mean.int.length.c = var(li.mean.int.length.c)

mean.li.mean.int.length.mcar = mean(li.mean.int.length.mcar)
var.li.mean.int.length.mcar = var(li.mean.int.length.mcar)

mean.li.mean.int.length.c.mcar = mean(li.mean.int.length.c.mcar)
var.li.mean.int.length.c.mcar = var(li.mean.int.length.c.mcar)

mean.li.prop.above.1 = mean(li.prop.above.1)
var.li.prop.above.1 = var(li.prop.above.1)

mean.li.prop.above.1.mcar = mean(li.prop.above.1.mcar)
var.li.prop.above.1.mcar = var(li.prop.above.1.mcar)

mean.li.prop.below.0 = mean(li.prop.below.0)
var.li.prop.below.0 = var(li.prop.below.0)

mean.li.prop.below.0.mcar = mean(li.prop.below.0.mcar)
var.li.prop.below.0.mcar = var(li.prop.below.0.mcar)

mean.li.cov.rate = mean(li.cov.rate)
var.li.cov.rate = var(li.cov.rate)

mean.li.cov.rate.mcar = mean(li.cov.rate.mcar)
var.li.cov.rate.mcar = var(li.cov.rate.mcar)

mean.li.prop.zero.width = mean(li.prop.zero.width)
var.li.prop.zero.width = var(li.prop.zero.width)

mean.li.prop.zero.width.mcar = mean(li.prop.zero.width.mcar)
var.li.prop.zero.width.mcar = var(li.prop.zero.width.mcar)

mean.li.mean.int.length.nzw = mean(li.mean.int.length.nzw)
var.li.mean.int.length.nzw = var(li.mean.int.length.nzw)

mean.li.mean.int.length.c.nzw = mean(li.mean.int.length.c.nzw)
var.li.mean.int.length.c.nzw = var(li.mean.int.length.c.nzw)

mean.li.mean.int.length.mcar.nzw = mean(li.mean.int.length.mcar.nzw)
var.li.mean.int.length.mcar.nzw = var(li.mean.int.length.mcar.nzw)

mean.li.mean.int.length.c.mcar.nzw = mean(li.mean.int.length.c.mcar.nzw)
var.li.mean.int.length.c.mcar.nzw = var(li.mean.int.length.c.mcar.nzw)

#Code for outputting results to an Excel .csv File
designation = c("Average Percentage of Wald Intervals That Go Below Zero",
                "Average Percentage of Wald Intervals That Go Above One",
                "Average Percentage of Zero-Width Wald Intervals",
                "Average Coverage Rate of Wald Intervals",
                "Average Mean Wald Interval Length",
                "Average Mean Truncated Wald Interval Length",
                "Average Mean Wald Interval Length, Not Including Zero-Width Intervals",
                "Average Mean Truncated Wald Interval Length, Not Including Zero-Width Intervals",
                "Average Percentage of Wilson Intervals That Go Below Zero",
                "Average Percentage of Wilson Intervals That Go Above One",
                "Average Percentage of Zero-Width Wilson Intervals",
                "Average Coverage Rate of Wilson Intervals",
                "Average Mean Wilson Interval Length",
                "Average Mean Truncated Wilson Interval Length",
                "Average Mean Wilson Interval Length, Not Including Zero-Width Intervals",
                "Average Mean Truncated Wilson Interval Length, Not Including Zero-Width Intervals",
                "Average of the Mean Percentage of Missing Data",
                "Mean Percentage of Ubar_m That Equal Zero",
                "Mean Percentage of B_m That Equal Zero",
                "Mean Percentage of Times that Both Ubar_m and B_m Equal Zero",
                "Average Percentage of Multiple Imputation Wald Intervals That Go Below Zero",
                "Average Percentage of Multiple Imputation Wald Intervals That Go Above One",
                "Average Percentage of Zero-Width Multiple Imputation Wald Intervals",
                "Average Coverage Rate of Multiple Imputation Wald Intervals",
                "Average Mean Multiple Imputation Wald Interval Length",
                "Average Mean Truncated Multiple Imputation Wald Interval Length",
                "Average Mean Multiple Imputation Wald Interval Length, Not Including Zero-Width Intervals",
                "Average Mean Truncated Multiple Imputation Wald Interval Length, Not Including Zero-Width Intervals",
                "Average Percentage of MI-Wilson Intervals That Go Below Zero",
                "Average Percentage of MI-Wilson Intervals That Go Above One",
                "Average Percentage of Zero-Width MI-Wilson Intervals",
                "Average Coverage Rate of MI-Wilson Intervals",
                "Mean Average Length of MI-Wilson Intervals",
                "Mean Average Truncated Length of MI-Wilson Intervals",
                "Mean Average Length of MI-Wilson Intervals, Not Including Zero-Width Intervals",
                "Mean Average Truncated Length of MI-Wilson Intervals, Not Including Zero-Width Intervals",
                "Average Percentage of Li et al. Intervals That Go Below Zero",
                "Average Percentage of Li et al. Intervals That Go Above 1",
                "Average Percentage of Li et al. Intervals That Are Zero-Width",
                "Average Coverage Rate of Li et al. Intervals",
                "Average Mean Length of Li et al. Intervals",
                "Average Mean Truncated Length of Li et al. Intervals",
                "Average Mean Length of Li et al. Intervals, Not Including Zero-Width Intervals",
                "Average Mean Truncated Length of Li et al. Intervals, Not Including Zero-Width Intervals",
                "Average Percentage of Plug-In Intervals With The Z-Statistic That Go Below Zero",
                "Average Percentage of Plug-In Intervals With The Z-Statistic That Go Above One",
                "Average Percentage of Zero-Width Plug-In Intervals With The Z-Statistic",
                "Average Coverage Rate of Plug-In Intervals With The Z-Statistic",
                "Average Mean Length of Plug-In Intervals With The Z-Statistic",
                "Average Mean Truncated Length of Plug-In Intervals With The Z-Statistic",
                "Average Mean Length of Plug-In Intervals With The Z-Statistic, Not Including Zero-Width Intervals",
                "Average Mean Truncated Length of Plug-In Intervals With The Z-Statistic, Not Including Zero-Width Intervals",
                "Average Percentage of Plug-In Intervals With The T-Statistic That Go Below Zero",
                "Average Percentage of Plug-In Intervals With The T-Statistic That Go Above One",
                "Average Percentage of Zero-Width Plug-In Intervals With The T-Statistic",
                "Average Coverage Rate of Plug-In Intervals With The T-Statistic",
                "Average Mean Length of Plug-In Intervals With The T-Statistic",
                "Average Mean Truncated Length of Plug-In Intervals With The T-Statistic",
                "Average Mean Length of Plug-In Intervals With The T-Statistic, Not Including Zero-Width Intervals",
                "Average Mean Truncated Length of Plug-In Intervals With The T-Statistic, Not Including Zero-Width Intervals",
                "Mean Percentage of MCAR Ubar_m That Equal Zero",
                "Mean Percentage of MCAR B_m That Equal Zero",
                "Mean Percentage of Times that Both MCAR Ubar_m and MCAR B_m Equal Zero",
                "Average Percentage of MCAR Multiple Imputation Wald Intervals That Go Below Zero",
                "Average Percentage of MCAR Multiple Imputation Wald Intervals That Go Above One",
                "Average Percentage of MCAR Zero-Width Multiple Imputation Wald Intervals",
                "Average Coverage Rate of MCAR Multiple Imputation Wald Intervals",
                "Average Mean MCAR Multiple Imputation Wald Interval Length",
                "Average Mean Truncated MCAR Multiple Imputation Wald Interval Length",
                "Average Mean MCAR Multiple Imputation Wald Interval Length, Not Including Zero-Width Intervals",
                "Average Mean Truncated MCAR Multiple Imputation Wald Interval Length, Not Including Zero-Width Intervals",
                "Average Percentage of MCAR MI-Wilson Intervals That Go Below Zero",
                "Average Percentage of MCAR MI-Wilson Intervals That Go Above One",
                "Average Percentage of Zero-Width MCAR MI-Wilson Intervals",
                "Average Coverage Rate of MCAR MI-Wilson Intervals",
                "Mean Average Length of MCAR MI-Wilson Intervals",
                "Mean Average Truncated Length of MCAR MI-Wilson Intervals",
                "Mean Average Length of MCAR MI-Wilson Intervals, Not Including Zero-Width Intervals",
                "Mean Average Truncated Length of MCAR MI-Wilson Intervals, Not Including Zero-Width Intervals",
                "Average Percentage of MCAR Li et al. Intervals That Go Below Zero",
                "Average Percentage of MCAR Li et al. Intervals That Go Above One",
                "Average Percentage of MCAR Li et al. Intervals That Are Zero-Width",
                "Average Coverage Rate of MCAR Li et al. Intervals",
                "Average Mean Length of MCAR Li et al. Intervals",
                "Average Mean Truncated Length of MCAR Li et al. Intervals",
                "Average Mean Length of MCAR Li et al. Intervals, Not Including Zero-Width Intervals",
                "Average Mean Truncated Length of MCAR Li et al. Intervals, Not Including Zero-Width Intervals",
                "Average Percentage of MCAR Plug-In Intervals With The Z-Statistic That Go Below Zero",
                "Average Percentage of MCAR Plug-In Intervals With The Z-Statistic That Go Above One",
                "Average Percentage of Zero-Width MCAR Plug-In Intervals With The Z-Statistic",
                "Average Coverage Rate of MCAR Plug-In Intervals With The Z-Statistic",
                "Average Mean Length of MCAR Plug-In Intervals With The Z-Statistic",
                "Average Mean Truncated Length of MCAR Plug-In Intervals With The Z-Statistic",
                "Average Mean Length of MCAR Plug-In Intervals With The Z-Statistic, Not Including Zero-Width Intervals",
                "Average Mean Truncated Length of MCAR Plug-In Intervals With The Z-Statistic, Not Including Zero-Width Intervals",
                "Average Percentage of MCAR Plug-In Intervals With The T-Statistic That Go Below Zero",
                "Average Percentage of MCAR Plug-In Intervals With The T-Statistic That Go Above One",
                "Average Percentage of Zero-Width MCAR Plug-In Intervals With The T-Statistic",
                "Average Coverage Rate of MCAR Plug-In Intervals With The T-Statistic",
                "Average Mean Length of MCAR Plug-In Intervals With The T-Statistic",
                "Average Mean Truncated Length of MCAR Plug-In Intervals With The T-Statistic",
                "Average Mean Length of MCAR Plug-In Intervals With The T-Statistic, Not Including Zero-Width Intervals",
                "Average Mean Truncated Length of MCAR Plug-In Intervals With The T-Statistic, Not Including Zero-Width Intervals")
value = c(mean.s.wald.prop.below.0,
          mean.s.wald.prop.above.1,
          mean.s.wald.prop.zero.width,
          mean.s.wald.cov.rate,
          mean.s.wald.mean.int.length,
          mean.s.wald.mean.int.length.c,
          mean.s.wald.mean.int.length.nzw,
          mean.s.wald.mean.int.length.c.nzw,
          mean.wilson.prop.below.0,
          mean.wilson.prop.above.1,
          mean.wilson.prop.zero.width,
          mean.wilson.cov.rate,
          mean.wilson.mean.int.length,
          mean.wilson.mean.int.length.c,
          mean.wilson.mean.int.length.nzw,
          mean.wilson.mean.int.length.c.nzw,
          mean.mean.prop.miss,
          mean.prop.zero.um,
          mean.prop.zero.bm,
          mean.prop.zero.bm.um,
          mean.wald.prop.below.0,
          mean.wald.prop.above.1,
          mean.wald.prop.zero.width,
          mean.wald.cov.rate,
          mean.wald.mean.int.length,
          mean.wald.mean.int.length.c,
          mean.wald.mean.int.length.nzw,
          mean.wald.mean.int.length.c.nzw,
          mean.wn.prop.below.0,
          mean.wn.prop.above.1,
          mean.wn.prop.zero.width,
          mean.wn.cov.rate,
          mean.wn.mean.int.length,
          mean.wn.mean.int.length.c,
          mean.wn.mean.int.length.nzw,
          mean.wn.mean.int.length.c.nzw,
          mean.li.prop.below.0,
          mean.li.prop.above.1,
          mean.li.prop.zero.width,
          mean.li.cov.rate,
          mean.li.mean.int.length,
          mean.li.mean.int.length.c,
          mean.li.mean.int.length.nzw,
          mean.li.mean.int.length.c.nzw,
          mean.wz.prop.below.0,
          mean.wz.prop.above.1,
          mean.wz.prop.zero.width,
          mean.wz.cov.rate,
          mean.wz.mean.int.length,
          mean.wz.mean.int.length.c,
          mean.wz.mean.int.length.nzw,
          mean.wz.mean.int.length.c.nzw,
          mean.wt.prop.below.0,
          mean.wt.prop.above.1,
          mean.wt.prop.zero.width,
          mean.wt.cov.rate,
          mean.wt.mean.int.length,
          mean.wt.mean.int.length.c,
          mean.wt.mean.int.length.nzw,
          mean.wt.mean.int.length.c.nzw,
          mean.prop.zero.um.mcar,
          mean.prop.zero.bm.mcar,
          mean.prop.zero.bm.um.mcar,
          mean.wald.prop.below.0.mcar,
          mean.wald.prop.above.1.mcar,
          mean.wald.prop.zero.width.mcar,
          mean.wald.cov.rate.mcar,
          mean.wald.mean.int.length.mcar,
          mean.wald.mean.int.length.c.mcar,
          mean.wald.mean.int.length.mcar.nzw,
          mean.wald.mean.int.length.c.mcar.nzw,
          mean.wn.prop.below.0.mcar,
          mean.wn.prop.above.1.mcar,
          mean.wn.prop.zero.width.mcar,
          mean.wn.cov.rate.mcar,
          mean.wn.mean.int.length.mcar,
          mean.wn.mean.int.length.c.mcar,
          mean.wn.mean.int.length.mcar.nzw,
          mean.wn.mean.int.length.c.mcar.nzw,
          mean.li.prop.below.0.mcar,
          mean.li.prop.above.1.mcar,
          mean.li.prop.zero.width.mcar,
          mean.li.cov.rate.mcar,
          mean.li.mean.int.length.mcar,
          mean.li.mean.int.length.c.mcar,
          mean.li.mean.int.length.mcar.nzw,
          mean.li.mean.int.length.c.mcar.nzw,
          mean.wz.prop.below.0.mcar,
          mean.wz.prop.above.1.mcar,
          mean.wz.prop.zero.width.mcar,
          mean.wz.cov.rate.mcar,
          mean.wz.mean.int.length.mcar,
          mean.wz.mean.int.length.c.mcar,
          mean.wz.mean.int.length.mcar.nzw,
          mean.wz.mean.int.length.c.mcar.nzw,
          mean.wt.prop.below.0.mcar,
          mean.wt.prop.above.1.mcar,
          mean.wt.prop.zero.width.mcar,
          mean.wt.cov.rate.mcar,
          mean.wt.mean.int.length.mcar,
          mean.wt.mean.int.length.c.mcar,
          mean.wt.mean.int.length.mcar.nzw,
          mean.wt.mean.int.length.c.mcar.nzw)
variance = c(var.s.wald.prop.below.0,
             var.s.wald.prop.above.1,
             var.s.wald.prop.zero.width,
             var.s.wald.cov.rate,
             var.s.wald.mean.int.length,
             var.s.wald.mean.int.length.c,
             var.s.wald.mean.int.length.nzw,
             var.s.wald.mean.int.length.c.nzw,
             var.wilson.prop.below.0,
             var.wilson.prop.above.1,
             var.wilson.prop.zero.width,
             var.wilson.cov.rate,
             var.wilson.mean.int.length,
             var.wilson.mean.int.length.c,
             var.wilson.mean.int.length.nzw,
             var.wilson.mean.int.length.c.nzw,
             var.mean.prop.miss,
             var.prop.zero.um,
             var.prop.zero.bm,
             var.prop.zero.bm.um,
             var.wald.prop.below.0,
             var.wald.prop.above.1,
             var.wald.prop.zero.width,
             var.wald.cov.rate,
             var.wald.mean.int.length,
             var.wald.mean.int.length.c,
             var.wald.mean.int.length.nzw,
             var.wald.mean.int.length.c.nzw,
             var.wn.prop.below.0,
             var.wn.prop.above.1,
             var.wn.prop.zero.width,
             var.wn.cov.rate,
             var.wn.mean.int.length,
             var.wn.mean.int.length.c,
             var.wn.mean.int.length.nzw,
             var.wn.mean.int.length.c.nzw,
             var.li.prop.below.0,
             var.li.prop.above.1,
             var.li.prop.zero.width,
             var.li.cov.rate,
             var.li.mean.int.length,
             var.li.mean.int.length.c,
             var.li.mean.int.length.nzw,
             var.li.mean.int.length.c.nzw,
             var.wz.prop.below.0,
             var.wz.prop.above.1,
             var.wz.prop.zero.width,
             var.wz.cov.rate,
             var.wz.mean.int.length,
             var.wz.mean.int.length.c,
             var.wz.mean.int.length.nzw,
             var.wz.mean.int.length.c.nzw,
             var.wt.prop.below.0,
             var.wt.prop.above.1,
             var.wt.prop.zero.width,
             var.wt.cov.rate,
             var.wt.mean.int.length,
             var.wt.mean.int.length.c,
             var.wt.mean.int.length.nzw,
             var.wt.mean.int.length.c.nzw,
             var.prop.zero.um.mcar,
             var.prop.zero.bm.mcar,
             var.prop.zero.bm.um.mcar,
             var.wald.prop.below.0.mcar,
             var.wald.prop.above.1.mcar,
             var.wald.prop.zero.width.mcar,
             var.wald.cov.rate.mcar,
             var.wald.mean.int.length.mcar,
             var.wald.mean.int.length.c.mcar,
             var.wald.mean.int.length.mcar.nzw,
             var.wald.mean.int.length.c.mcar.nzw,
             var.wn.prop.below.0.mcar,
             var.wn.prop.above.1.mcar,
             var.wn.prop.zero.width.mcar,
             var.wn.cov.rate.mcar,
             var.wn.mean.int.length.mcar,
             var.wn.mean.int.length.c.mcar,
             var.wn.mean.int.length.mcar.nzw,
             var.wn.mean.int.length.c.mcar.nzw,
             var.li.prop.below.0.mcar,
             var.li.prop.above.1.mcar,
             var.li.prop.zero.width.mcar,
             var.li.cov.rate.mcar,
             var.li.mean.int.length.mcar,
             var.li.mean.int.length.c.mcar,
             var.li.mean.int.length.mcar.nzw,
             var.li.mean.int.length.c.mcar.nzw,
             var.wz.prop.below.0.mcar,
             var.wz.prop.above.1.mcar,
             var.wz.prop.zero.width.mcar,
             var.wz.cov.rate.mcar,
             var.wz.mean.int.length.mcar,
             var.wz.mean.int.length.c.mcar,
             var.wz.mean.int.length.mcar.nzw,
             var.wz.mean.int.length.c.mcar.nzw,
             var.wt.prop.below.0.mcar,
             var.wt.prop.above.1.mcar,
             var.wt.prop.zero.width.mcar,
             var.wt.cov.rate.mcar,
             var.wt.mean.int.length.mcar,
             var.wt.mean.int.length.c.mcar,
             var.wt.mean.int.length.mcar.nzw,
             var.wt.mean.int.length.c.mcar.nzw)
missingness = c(rep("no missing", 16), rep("MAR", 44),
                rep("MCAR", 43))
value.percent = value*100
variance.percent = variance*100
p = rep(0.05, 103)
proportion.missing = rep(0.3, 103)
n = rep(500, 103)

df = data.frame(p = p, n = n, proportion.missing = proportion.missing,
                missingness = missingness, designation = designation,
                value.percent = value.percent,
                variance.percent = variance.percent)
colnames(df) = c("p", "n", "proportion missing", "missingness category",
                 "designation", "value (in percent)", "variance (in percent)")
write.csv(df,
           file = "p05n500mp3_ind_div_by_3.csv")

#Function to compare the variance of Q_m with the mean of T_m
multimp_var2 = function(n, p, a, b, it, m, mp){
  #n is the number of binomial trials
  #p is the true proportion of successes
  #a is the first beta prior parameter
  #b is the second beta prior parameter
  #it is the number of iterations of multiple imputations
  #m is the number of imputed data sets
  #mp is the proportion of missing data for MCAR
  qm_vec = rep(NA,it)
  Tm_vec = rep(NA,it)
  qm_mcar_vec = rep(NA, it)
  Tm_mcar_vec = rep(NA, it)
  nonzero.width = rep(NA, it)
  for (i in 1:it){
    #Initialize x.d, generate y.d
    x.d = rep(NA, n)
    y.d = rbinom(n, 1, p)
    mcar.data = y.d
    no.missing.data = y.d
    #Generate x.d based on y.d
    for (u in 1:length(y.d)){
      if (y.d[u] ==1){
        x.d[u] = rbinom(1, 1, 0.2)
      }else {
        x.d[u] = rbinom(1,1, 0.6)
      }
    }
    
    #Create vector of indices indicating
    #missingness based on x.d
    miss = rep(NA, n)
    for (k in 1:length(x.d)){
      if (x.d[k] == 1){
        miss[k] = rbinom(1,1, 0.18) 
      }else{                       
        miss[k] = rbinom(1,1, 0.47)
      }
    }
    
    
    #Seperately compiling all observed y.d's 
    #corresponding to x.d = 1 and x.d = 0 (where
    #missingness index equals 0)
    y_x1 = c()
    y_x0 = c()
    for (j2 in 1:length(x.d)){
      if (miss[j2] == 0 & x.d[j2] ==1){
        y_x1 = c(y_x1, y.d[j2])
      }else if (miss[j2] == 0 & x.d[j2] == 0){
        y_x0 = c(y_x0, y.d[j2])
      }
    }
    
    #Determining r, the proportion of missingness
    count = 0 
    for (xi in 1:length(miss)){
      if (miss[xi] ==1){
        count = count + 1
      }
    }
    prop_miss = count/n
    
    #Initializing vectors for multiple imputation
    p_hat = rep(NA, m)
    var_p = rep(NA, m)
    #Below: creating m multiply imputed data sets,
    #each data set with missing y.d imputed based on 
    #whether x.d = 1 or x.d = 0 for that index,
    #with different imputation parameters for the 
    #two conditions
    an_x1 = a + sum(y_x1)
    bn_x1 = b + length(y_x1) - sum(y_x1)
    an_x0 = a + sum(y_x0)
    bn_x0 = b + length(y_x0) - sum(y_x0)
    for (j in 1:m) {
      theta_x1 = rbeta(1, an_x1, bn_x1)
      theta_x0 = rbeta(1, an_x0, bn_x0)
      for (k2 in 1:length(y.d)){
        if (miss[k2] ==1 & x.d[k2] == 1){
          y.d[k2] = rbinom(1, 1, theta_x1)
        }else if (miss[k2] == 1 & x.d[k2] == 0){
          y.d[k2] = rbinom(1,1, theta_x0)
        }
      }
      p_hat[j] = mean(y.d)
      var_p[j] = p_hat[j]*(1-p_hat[j])/n
    }
    #Rubin's rules for MAR data
    qm= sum(p_hat)/m
    um = sum(var_p)/m
    bm = var(p_hat)
    Tm = ((1 + (1/m))*bm) + um
    if (bm == 0){
      v = 1000000000
      rm = 0 #convert to rm = 0, for both bm = 0 and um = 0 and 
      #bm = 0 and um != 0
    }else{
      v = (m-1)*(1+(um/((1+ (1/m))*bm)))^2 #convert to normal if bm = 0 and um != 0
      rm = (1 + (1/m))*(bm/um)
    }
    
    #Multiple Imputation for MCAR data
    s = n*(1-mp)
    s_plus_1 = s + 1
    missing_vec = c(s_plus_1:n)
    obs = mcar.data[-missing_vec]
    an.mcar = a + sum(obs)
    bn.mcar = b + length(obs) - sum(obs)
    p.hat.mcar = rep(NA, m)
    var.p.mcar = rep(NA, m)
    for (ji in 1:m) {
      theta.mcar = rbeta(1, an.mcar, bn.mcar)
      imp.mcar = rbinom(n*mp, 1, theta.mcar)
      binom_data = c(obs, imp.mcar)
      p.hat.mcar[ji] = mean(binom_data)
      var.p.mcar[ji] = p.hat.mcar[ji]*(1-p.hat.mcar[ji])/n
    }
    #Rubin's rules for MCAR data
    qm.mcar= sum(p.hat.mcar)/m
    um.mcar = sum(var.p.mcar)/m
    bm.mcar = var(p.hat.mcar)
    Tm.mcar = ((1 + (1/m))*bm.mcar) + um.mcar
    if (bm.mcar == 0){
      v.mcar = 1000000000
      rm.mcar = 0
    }else{
      v.mcar = (m-1)*(1+(um.mcar/((1+ (1/m))*bm.mcar)))^2
      rm.mcar = (1 + (1/m))*(bm.mcar/um.mcar)
    }
    qm_vec[i] = qm
    Tm_vec[i] = Tm
    qm_mcar_vec[i] = qm.mcar
    Tm_mcar_vec[i] = Tm.mcar
    
  }
  return(data.frame(Tm = var(qm_vec),Tm_mean= mean(Tm_vec),
                    Tm.mcar = var(qm_mcar_vec), Tm.mcar.mean = mean(Tm_mcar_vec)))
}
#input your value of p repeated 1000 times
x= rep(0.5, 100)
y = sapply(x, multimp_var2,n=100, a=1, b=1, it = 1000, m=10)
#check the proportion of times that the mean of the Tms
#is greater than the variance of the qms.
#find the average proportion by which the mean of the Tms is 
#greater than the variance of the qms.
logical = rep(NA, length(x))
pgreater = rep(NA, length(x))
logical.mcar = rep(NA, length(x))
pgreater.mcar = rep(NA, length(x))
for (l in 1:length(x)){
  logical[l] = (y[[2,l]] > y[[1,l]])
  pgreater[l] = (y[[2,l]] - y[[1,l]])/y[[1,l]]
  logical.mcar[l] = (y[[4, l]] > y[[3, l]])
  pgreater.mcar[l] = (y[[4,l]] -y[[3, l]])/y[[3, l]]
}
varpgreater = var(pgreater)
base_percent = length(logical[logical==TRUE])/length(logical)
percent_greater = mean(pgreater)
percent_greater_25 = length(pgreater[pgreater> 0.25])/length(pgreater)

varpgreater_mcar = var(pgreater.mcar)
base_percent_mcar = length(logical.mcar[logical.mcar==TRUE])/length(logical.mcar)
percent_greater_mcar = mean(pgreater.mcar)
percent_greater_25_mcar = length(pgreater.mcar[pgreater.mcar> 0.25])/length(pgreater.mcar)