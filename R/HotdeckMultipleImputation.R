
#hotdeckMI: R function for making hotdeck multiple imputation data sets
#           as described in Sun, Qi, Heng, Gilbert (2020, JRSS-C)
#           Peter Gilbert May 18, 2023; reviewed by Michal Juraska May 19, 2023.

#Data set: Notation as described in Sun et al., a vaccine efficacy trial,
#data from the vaccine and placebo arm are included, where only a single baseline stratum K=1 is considered.
#This feature could be readily updated if needed.

#T = time to failure outcome (typically the time origin is the visit at which the
#    immune marker was measured or several days later.
#Delta = Indicator that T is observed (Delta=1) where Delta=0 means right-censored.
#X = Observed failure time (minimum of T and right-censoring time) (where X=T iff Delta=1)
#Z1discrete = Vector of covariates measured in everyone that are discrete for which complete matching is used for hotdeck MI 
#Z1scalar   = Vector of covariates measured in everyone that are continuous/scalar that are used for hotdeck MI 
#             through the distance function h(.,.)
#Z2 = The single covariate of interest that is measured in a subset of participants (e.g., the immune marker of interest) 
#     (the phase two covariate of interest)
#epsilonz = Indicator of whether Z2 is observed (1=observed, 0=missing).
#V = Mark variable (viral genotype or immunophenotype) only possibly measured in cases (that have Delta = 1).  
#    V could be binary or an integer or a scalar, an arbitrary type of univariable mark. V is NA if Delta = 0. 
#epsilonv = Indicator whether V is observed (1=observed, 0=missing).  For participants with Delta=0, epsilonv=1.
#Avdiscrete = Vector of 'auxiliary' mark features only possibly measured in cases (with Delta=1) that can be used to 
#             predict V and to help MAR for sequence missingness hold, restricting to variables that are discrete and 
#             for which complete matching is used for hotdeck MI.  Avdiscrete is a vector of all NAs if it is not used.
#Avscalar     Vector of 'auxiliary' mark features only possibly measured in cases (with Delta=1) that can be used to 
#             predict V, restricting to variables that are continuous/scalar that are used for hotdeck MI through the 
#             distance function h(.,.)
#M            Number of outputted data sets V^{(m)}_{ki}, m=1, \cdots, M
#L            Number of nearest neighbors
#
# Yanging Sun: 05/19/23: Verified that a step to bootstrap all cases is needed to make imputations proper
# (Murray, 2018, Statistical Science).  M bootstrapped data sets are generated, corresponding to the M outputted data sets.
# Note that a separate bootstrap data set is not generated for each case with a missing mark.
#
# For hotdeck MI to be valid, need to include in the matching variables all factors that impact whether the mark is observed.
# If a case-cohort/case-control correlates study is in view, and only a subset of cases is sampled, then the design variables
# influencing which subset is sampled must be included in the matching variables to yield proper imputations.
#
# Michal Juraska 05/18/23 QUESTION: IT SEEMS THAT CURRENTLY AVSCALAR MUST BE PROVIDED.
# IT WOULD BE MORE FLEXIBLE TO ALSO ALLOW NOT USING AVSCALAR (AS IS ALLOWED FOR AVDISCRETE).
# PG 05/19/23: For now leave this feature out.

#epsilona   = Indicator of whether Avdiscrete and Avscalar are available from a case (with Delta = 1).  
#             For participants with Delta=0, epsilona=1.

#########################################################################################################################
#The following summarizes some hotdeck multiple imputation notation from Section 2.1 of Sun et al. 2020
#(where k=1 always for the R script, so that k notation does not add information).
#Let i index study participant and k(=1) index baseline strata.
#The code considers 4 possible sets of available variables used for hotdeck imputation for a given case with 
#a missing mark, where in general all available variables are used:

#Case 1. H_{ki} = (T_{ki},Z1_{ki},Z2_{ki},A_{v,ki}) [Z2_{ki} is available and A_vi is available]
#Case 2. H_{ki} = (T_{ki},Z1_{ki},Z2_{ki})          [Z2_{ki} is available and A_vi is not available]
#Case 3. H_{ki} = (T_{ki},Z1_{ki},A_{v,ki})         [Z2_{ki} is not available and A_vi is available]
#Case 4. H_{ki} = (T_{ki},Z1_{ki})                  [Z2_{ki} is not available and A_vi is not available] 

# d(H_{ki},H_{kj}) is a distance between the two vectors H_{ki} and H_{kj}, where we use z-scores of Euclidean distance.

# Distances are only calculated for continuous/scalar H_{ki}; for discrete unordered H_{ki} exact matches are used.
# Exact matches of all discrete variables are required, otherwise only the continuous/scalar H_{ki} are used for 
# defining nearest neighbors.  There is one exception: In Case 1, if a perfect match on the discrete components of Z1_{ki}
# and A_{v,ki} is not possible, then a second chance is taken to seek a perfect match on the discrete components of
# A_{v,ki} only.

# Sun et al. describes use of the methhods based on a fixed L (amounting to variable bandwidh) or to fixed B (amounting to
#fixed bandwidth).  This code considers a fixed L only.

######################################################################
#Inputs into the function hotdeckMI:

#1. (X,Delta,Z1discrete,Z1scalar,Z2,epsilonz,V,epsilonv,Avdiscrete,Avscalar,epsilonv) 
#   where Z2, V, and Av have missing values satisfying
#   epsilonz==0 iff Z2 is missing;
#   epsilonv==0 iff (Delta==1 & V is missing);
#   Av is missing iff Delta==0 | (Delta==1 & epsilona==0);
#   if Delta==1 then Avscalar is missing iff epsilona==0.

#2. M = number of outputted data sets V^{(m)}_{ki}, m=1, \cdots, M.
#   Default value M=10.

#3. L = A counting number, the number of nearest neighbors, where L must be less than the number of 
#       cases with observed marks.
#   Default value L=5.

#####################################################
#Outputs of the function hotdeckMI:
#
#Vmat: This is the notation V^{(m)}_i in the paper.  Vmat is a matrix with NN rows and M columns, where 
#      NN is the number of participants (e.g., number of per-protocol vaccine recipients observed to be at-risk at the 
#      visit used for sampling of the immune marker plus number of per-protocol placebo recipients at-risk by the same time point).
#      A genotype/mark V is determined for 100% of 
#      cases (Delta=1), where cases with missing genotypes have M marks V imputed from nearest-neighbor cases with 
#      observed V values, and cases with observed V values simply have their observed genotype included M times.

#Zscoredistsnearest: a matrix with NN rows and L columns, which for each participant lists the L nearest distance values
#                    h(.,.).  Note that particicipants with an observed mark V simply have a vector of L NA values.
#                    This output is for rough diagnostic and code testing purposes.  If fewer than L nearest neighbors are
#                    available, then the number of columns with non-missing distances is the number available.

#####################################################
# Notes on the hotdeckMI function:
#
# 1. The hotdeck imputed marks are specific to a given immune marker studied as a correlate, because
#    that immune marker is included in the hotdeck, which is important for MAR being valid.  For multivariable marker analysis, 
#    an alternative approach would do the hotdeck including all of the markers in Avscalar could be advantageous, but that is 
#    more complicated if the markers have differnt missing data patterns and is not included.
# 
# 2. The function only works if there is at least one case with an observed mark in each of the 4 categories
#    defined by the cross-classification of epsilonz & epsilona (Case 1, 2, 3, 4 above).
#
# 3. The method requires 100% matches of neighbors on Z1discrete and Avdiscrete, except if not matches are available
#    it will also accept 100% matches on Avdiscrete but not on Z1discrete.  This feature is included because Avdiscrete
#    may be an important variable such as geographic region.
#
# 4. If wish to NOT require a perfect match on Z1discrete covariates, enter Z1discrete as a vector of all constant values,
#    e.g., Z1discrete <- matrix(rep(1,nrow(dat)),ncol=1).  This can be important if Z1discrete is 'weakly valuable'
#    for matching whereas in contrast other values are very important (like calendar date of case, geography of case)
#
# 5. If wish to NOT have X be used in the hotdeck distances to neighbors, enter X as a vector of constant values,
#    e.g., X <- rep(1,nrow(dat)).  
#
# 6. Similarly if wish all of Z1scalar, or elements of, to NOT be used in the hotdeck distances to neighbor, set the
#    relevant components to be a vector of constant values, e.g., Z1scalar <- matrix(rep(1,nrow(dat)),ncol=1)
#


hotdeckMI <- function(X,Delta,Z1discrete,Z1scalar,Z2,epsilonz,V,epsilonv,Avdiscrete,Avscalar,epsilona,M=10,L=5) {
# Input parameter requirements:
# X, Delta, Z2, epsilonz,V,epsilonv, epsilona are vectors, with no missing values allowed for X, Delta, epsilonz, epsilonv, epsilona.
# Z2 is missing iff epsilonz==0, V is missing iff epsilonv==0.
# Z1discrete, Z1scalar, Avdiscrete, Avscalar are matrices.  Z1discrete and Z1scalar have no missing values.
# Avdiscrete and Avscalar are missing iff epsilona==0.  Note the 'all-or-nothing' missing data pattern of sets of variables.

# Input data set checks
if (length(Z2[epsilonz==1 & is.na(Z2)])>0) { print("Problem 1")}
if (length(Z2[epsilonv==1 & Delta==1 & is.na(V)])>0) { print("Problem 2")}

# Make a distance function which for an individual subject i, calculates d(H_{ki},H_{kj}) for all j!=i (relative to the fixed i), 
# where j indexes all eligible neighbor cases that have the mark V observed.  H_{ki} and H_{kj} are vectors with p elements.
# The function distsi inputs the vector x=H_{ki} (1 by p vector for subject i) and a matrix Y of H_{kj} values Y (m by p) 
# where m is the number of eligible neighbor cases indexed to subject i.

distsi <- function(x,ymat) {
ymat <- ymat
m <- nrow(ymat)
p <- ncol(ymat)
xmat <- t(as.matrix(x))

# First remove columns with no variability across subjects pooling over x and ymat:
kpcols <- rep(TRUE,p)
for (j in 1:p) {
if (var(c(xmat[1,j],ymat[,j]))==0) {
kpcols[j] <- FALSE }}

p <- length(kpcols[kpcols])
if (p==0) {print("Problem in distsi with p=0") }

xmat <- t(as.matrix(x[kpcols]))
ymat <- ymat[,kpcols]
if (m==1) {ymat <- matrix(ymat,ncol=p) }
if (p==1) { ymat <- matrix(ymat,ncol=p) }

# First transform each of the p variables to be z-scores:
xZscore <- xmat
ymatZscore <- ymat

for (j in 1:p) {
mn <- mean(c(xmat[1,j],ymat[,j]))
sd <- sqrt(var(c(xmat[1,j],ymat[,j])))
xZscore[j] <- (xmat[1,j] - mn)/sd
ymatZscore[,j] <- (ymat[,j]-mn)/sd }

# Second compute Euclidean distances:
Eucdists <- rep(NA,m) 
for (i in 1:m) {
ans <- 0
for (j in 1:p) {
ans <- ans + (xZscore[j]-ymatZscore[i,j])^2 }
Eucdists[i] <- sqrt(ans) }
return(Eucdists) }

# Helper function that inputs a matrix of TRUE/FALSE values and returns a vector indicating whether all columns are TRUE

alltrue <- function(r, matt){
  return(apply(matt, 1, function(x, y){ identical(x, y) }, y=matt[r, ]))
}

NN <- length(Delta) # Number of rows/participants in the data set
n <- sum(Delta) # Number of cases
nmissmark <- sum(Delta*(1-epsilonv)) # Number of cases with a missing mark
Vmat <- matrix(NA,nrow=NN,ncol=M)
Zscoredistsnearest <- matrix(rep(NA,NN*L),nrow=NN)

# Fill in Vmat for cases with an observed mark:
for (i in 1:NN) {
if (epsilonv[i]==1) { Vmat[i,] <- rep(V[i],M) }
}

# Take a bootstrap sample of all cases, for implementing hotdeck based on bootstrapped cases
set.seed(12)
inds <- sample(c(1:NN)[Delta==1],replace=T)
# Re-define the data set based on this bootstrapped sample
Xbt <- X 
Z1discretebt <- Z1discrete 
Z1scalarbt <- Z1scalar 
Z2bt <- Z2
epsilonzbt <- epsilonz
Vbt <- V
epsilonvbt <- epsilonv
Avdiscretebt <- Avdiscrete
Avscalarbt <- Avscalar
epsilonabt <- epsilona
Deltabt <- Delta

Xbt[Delta==1] <- X[inds]
Z1discretebt[Delta==1,] <- Z1discrete[inds,]
Z1scalarbt[Delta==1,]   <- Z1scalar[inds,]
Z2bt[Delta==1] <- Z2[inds]
epsilonzbt[Delta==1] <- epsilonz[inds]
Vbt[Delta==1] <- V[inds]
epsilonvbt[Delta==1] <- epsilonv[inds]
Avdiscretebt[Delta==1,] <- Avdiscrete[inds,]
Avscalarbt[Delta==1,] <- Avscalar[inds,]
epsilonabt[Delta==1] <- epsilona[inds]
Deltabt[Delta==1] <- Delta[inds]

# Loop to create the between-imputation component m=1,...,M (Ruben and Shenker, 1986, JASA)
for (m in 1:M) {

# Loop over all observed cases with a missing mark, the set for which hotdeck imputations are needed
for (i in c(1:length(Delta))[Delta==1 & epsilonv==0]) {
cat(paste("m,i=",c(m,i)),"\n")

################################################################################################
# Case 1: Z2 measured (epsilonz==1) and Avdiscrete, Avscalar are available (epsilona==1) 

if (epsilonz[i]==1 & epsilona[i]==1) { 
# Hotdeck matching variables = (X,Z1discrete,Z1scalar,Z2,Avdiscrete,Avscalar)
# First, define cases with observed marks completely matched on the discrete covariates:
if (length(Avdiscretebt[,1][is.na(Avdiscretebt[,1])])==length(Avdiscretebt[,1])) { 
# Avdiscrete is not used
availableneighborcases <- Deltabt==1 & epsilonzbt==1 & epsilonabt==1 & epsilonvbt==1 & alltrue(i, Z1discretebt) & c(1:NN)!=i
} else {
# Avdiscrete is used
availableneighborcases <- Deltabt==1 & epsilonzbt==1 & epsilonabt==1 & epsilonvbt==1 & alltrue(i, Z1discretebt) & alltrue(i, Avdiscretebt) & c(1:NN)!=i }

if (sum(availableneighborcases[availableneighborcases])==0) {
# Second, if could not find complete matches based on all of the discrete variables, drop the 
# requirement for a complete match on Z1discrete, first aiming for a complete match on Avdiscrete:
if (!length(Avdiscretebt[,1][is.na(Avdiscretebt[,1])])==length(Avdiscretebt[,1])) { 
availableneighborcases <- Deltabt==1 & epsilonzbt==1 & epsilonabt==1 & epsilonvbt==1 & alltrue(i, Avdiscretebt) & c(1:NN)!=i }
else {availableneighborcases <- Deltabt==1 & epsilonzbt==1 & epsilonabt==1 & epsilonvbt==1 & c(1:NN)!=i }

# Michal Juraska 05/18/23 PROBLEM: PERFECT MATCHING ON Z1DISCRETE ALONE WHEN AVDISCRETE IS AVAILABLE IS NOT CONSIDERED. SHOULD IT BE IF
# NO MATCHES ARE ACHIEVED FOR EITHER Z1DISCRETE AND AVDISCRETE COMBINED OR AVDISCRETE ALONE?
# PG 05/19/23: Added the feature to allow the case of complete match on Z1discrete when a complete match on Avdiscrete fails.
 
#Third, if could not find complete matches based on Avdiscrete alone, check if there is a complete match based on Z1discrete alone:
if (sum(availableneighborcases[availableneighborcases])==0) {
availableneighborcases <- Deltabt==1 & epsilonzbt==1 & epsilonabt==1 & epsilonvbt==1 & alltrue(i, Z1discretebt) & c(1:NN)!=i }

#Fourth, if could not find complete matches based on Z1discrete alone, consider all cases:
if (sum(availableneighborcases[availableneighborcases])==0) {
availableneighborcases <- Deltabt==1 & epsilonzbt==1 & epsilonabt==1 & epsilonvbt==1 & c(1:NN)!=i }
if (sum(availableneighborcases[availableneighborcases])==0) { print("Problem 3 Case 1") }
}


# Michal Juraska 05/18/23 QUESTION: AVSCALAR IS ASSUMED TO BE ALWAYS OBSERVED WHEN AVDISCRETE IS OBSERVED.
# IT WOULD BE MORE FLEXIBLE TO ALSO ALLOW FOR AVDISCRETE TO BE AVAILABLE WHILE AVSCALAR IS UNAVAILABLE.
# THIS APPLIES TO BOTH CASES 1 AND 3.
# PG 05/19/23: For now, do not add this feature. Depending on viral load availability for the Sanofi Pasteur study,
# this may be added later.
  
vecti <- c(X[i],Z1scalar[i,],Z2[i],Avscalar[i,])
mattminusi <- cbind(Xbt[availableneighborcases],Z1scalarbt[availableneighborcases,],Z2bt[availableneighborcases],Avscalarbt[availableneighborcases,])
if(length(availableneighborcases[availableneighborcases])==1) {
mattminusi <- matrix(c(Xbt[availableneighborcases],Z1scalarbt[availableneighborcases,],Z2bt[availableneighborcases],Avscalarbt[availableneighborcases,]),
              ncol=length(vecti)) }
indsminusi <- c(1:NN)[availableneighborcases]
#print("Case 1 vecti:")
#cat(paste(vecti),"\n")
#print("Case 1 mattminusi:")
#cat(paste(mattminusi),"\n")
distsavailableneighborcases <- distsi(vecti,mattminusi)
# Identify the indices of the L smallest distances
Lnearestneighbors <- indsminusi[order(distsavailableneighborcases)][1:L]
# If length(indsminusi) < L, then NAs are introduced into Lnearestneighbors.  So only sample from the
# nonmissing elements of Lnearestneighbors
Lnearestneighbors <- Lnearestneighbors[!is.na(Lnearestneighbors)]
Vmat[i,m] <- sample(Vbt[Lnearestneighbors],size=1,replace=T)

# Youyi troubleshoot
# Are Lnearestneighbors ever missing?
print("Lnearestneighbors")
cat(paste(Vbt[Lnearestneighbors]))

# Troubleshoot an idisyncrosy with cbind:
#if (i==753) {
#cat(paste("length(vecti)=",length(vecti)),"\n")
#cat(paste("dim(mattminusi)=",dim(mattminusi)),"\n")
#cat(paste("indsminui=",indsminusi),"\n")
#cat(paste("table(availableneighborcases)=",table(availableneighborcases)),"\n")
#cat(paste("Avscalarbt[availableneighborcases,]=",Avscalarbt[availableneighborcases,]),"\n")
#cat(paste("Xbt[availableneighborcases,]=",Xbt[availableneighborcases]),"\n")
#cat(paste("Z1scalarbt[availableneighborcases,]=",Z1scalarbt[availableneighborcases]),"\n")
#cat(paste("Z2bt[availableneighborcases,]=",Z2bt[availableneighborcases]),"\n")
#}

# Printing for testing the function
cat(paste("Z2 observed and Av observed i:"),"\n")
cat(paste("Subject i X, Z1scalar, Z2, Avscalar = "),"\n")
print(vecti)
cat("\n")

cat(paste("Nearest neighbors X, Z1scalar, Z2, Avscalar = "),"\n")
print(mattminusi[order(distsavailableneighborcases)[1:L],])
cat(paste("\n"))

cat(paste("Farthest neighbors X, Z1scalar, Z2, Avscalar = "),"\n")
print(mattminusi[order(-distsavailableneighborcases)[1:L],])

if (length(vecti)!=ncol(mattminusi)) { print("Problem 7") 
}

#Note: NAs are not problematic; they merely indicate that there are fewer than L eligible cases
}

###############################################################################################
# Case 2: Z2 measured (epsilonz==1) and Avdiscrete, Avscalar are not available (epsilona==0)

if (epsilonz[i]==1 & epsilona[i]==0) { 
# Matching variables = (X,Z1discrete,Z1scalar,Z2)
# First define cases with observed marks completely matched on the discrete covariates:
availableneighborcases <- Deltabt==1 & epsilonzbt==1 & epsilonabt==0 & epsilonvbt==1 & alltrue(i, Z1discretebt) & c(1:NN)!=i  

if (sum(availableneighborcases[availableneighborcases])==0) {
# Second, if could not find complete matches based on the discrete covariates, consider all cases:
# Repeats the identical code as above:
availableneighborcases <- Deltabt==1 & epsilonzbt==1 & epsilonabt==0 & epsilonvbt==1 & c(1:NN)!=i
if (sum(availableneighborcases[availableneighborcases])==0) { print("Problem 3 Case 2") }
}

vecti <- c(X[i],Z1scalar[i,],Z2[i])
mattminusi <- cbind(Xbt[availableneighborcases],Z1scalarbt[availableneighborcases,],Z2bt[availableneighborcases])
if(length(availableneighborcases[availableneighborcases])==1) {
mattminusi <- matrix(c(Xbt[availableneighborcases],Z1scalarbt[availableneighborcases,],Z2bt[availableneighborcases]),
              ncol=length(vecti)) }
indsminusi <- c(1:NN)[availableneighborcases]
distsavailableneighborcases <- distsi(vecti,mattminusi)
# Identify the indices of the L smallest distances
Lnearestneighbors <- indsminusi[order(distsavailableneighborcases)][1:L]
# If length(indsminusi) < L, then NAs are introduced into Lnearestneighbors.  So only sample from the
# nonmissing elements of Lnearestneighbors
Lnearestneighbors <- Lnearestneighbors[!is.na(Lnearestneighbors)]
Vmat[i,m] <- sample(Vbt[Lnearestneighbors],size=1,replace=T)

cat(paste("Z2 observed and Av not observed i:"),"\n")
cat(paste("Subject i X, Z1scalar, Z2 = "),"\n")
print(vecti)
cat("\n")

cat(paste("Nearest neighbors X, Z1scalar, Z2 = "),"\n")
print(mattminusi[order(distsavailableneighborcases)[1:L],])
cat(paste("\n"))

cat(paste("Farthest neighbors X, Z1scalar, Z2 = "),"\n")
print(mattminusi[order(-distsavailableneighborcases)[1:L],])

}

# Case 3: Z2 not measured (epsilonz==0) and has the auxiliaries Avdiscrete, Avscalar measured (epsilona==1)

if (epsilonz[i]==0 & epsilona[i]==1) { 
# Matching variables = (X,Z1discrete,Z1scalar,Avdiscrete,Avscalar)
# First define cases with observed marks completely matched on the discrete covariates:
if (length(Avdiscretebt[,1][is.na(Avdiscretebt[,1])])==length(Avdiscretebt[,1])) {
# Avdiscrete is not used
availableneighborcases <- Deltabt==1 & epsilonzbt==0 & epsilonabt==1 & epsilonvbt==1 & alltrue(i, Z1discretebt) & c(1:NN)!=i
} else {
# Avdiscrete is used
availableneighborcases <- Deltabt==1 & epsilonzbt==0 & epsilonabt==1 & epsilonvbt==1 & alltrue(i, Z1discretebt) & alltrue(i, Avdiscretebt) & c(1:NN)!=i  }

if (sum(availableneighborcases[availableneighborcases])==0) {
# Second, if could not find complete matches based on the discrete covariates, consider all cases:
# Repeats the identical code as above:
availableneighborcases <- Deltabt==1 & epsilonzbt==0 & epsilonabt==1 & epsilonvbt==1 & c(1:NN)!=i
if (sum(availableneighborcases[availableneighborcases])==0) { print("Problem 3 Case 3") }
}

vecti <- c(X[i],Z1scalar[i,],Avscalar[i,])
mattminusi <- cbind(Xbt[availableneighborcases],Z1scalarbt[availableneighborcases,],Avscalarbt[availableneighborcases,])
if(length(availableneighborcases[availableneighborcases])==1) {
mattminusi <- matrix(c(Xbt[availableneighborcases],Z1scalarbt[availableneighborcases,],Avscalarbt[availableneighborcases,]),
              ncol=length(vecti)) }
indsminusi <- c(1:NN)[availableneighborcases]
distsavailableneighborcases <- distsi(vecti,mattminusi)
# Identify the indices of the L smallest distances
Lnearestneighbors <- indsminusi[order(distsavailableneighborcases)][1:L]
# If length(indsminusi) < L, then NAs are introduced into Lnearestneighbors.  So only sample from the
# nonmissing elements of Lnearestneighbors
Lnearestneighbors <- Lnearestneighbors[!is.na(Lnearestneighbors)]
# If length(indsminusi) < L, then NAs are introduced into Lnearestneighbors.  So only sample from the
# nonmissing elements of Lnearestneighbors
Lnearestneighbors <- Lnearestneighbors[!is.na(Lnearestneighbors)]
Vmat[i,m] <- sample(Vbt[Lnearestneighbors],size=1,replace=T)

cat(paste("Z2 not observed and Av observed i:"),"\n")
cat(paste("Subject i X, Z1scalar, Avscalar = "),"\n")
print(vecti)
cat("\n")

cat(paste("Nearest neighbors X, Z1scalar, Avscalar = "),"\n")
print(mattminusi[order(distsavailableneighborcases)[1:L],])
cat(paste("\n"))

cat(paste("Farthest neighbors X, Z1scalar, Avscalar = "),"\n")
print(mattminusi[order(-distsavailableneighborcases)[1:L],])

}

# Case 4: Z2 not measured (epsilonz==0) and Avdiscrete, Avscalar are not available (epsilona==0)

if (epsilonz[i]==0 & epsilona[i]==0) { 
# Matching variables = (X,Z1discrete,Z1scalar)
# First define cases with observed marks completely matched on the discrete covariates:
availableneighborcases <- Deltabt==1 & epsilonzbt==0 & epsilonabt==0 & epsilonvbt==1 & alltrue(i, Z1discretebt) & c(1:NN)!=i  

if (sum(availableneighborcases[availableneighborcases])==0) {
# Second, if could not find complete matches based on the discrete covariates, consider all cases:
# Repeats the identical code as above:
availableneighborcases <- Deltabt==1 & epsilonzbt==0 & epsilonabt==0 & epsilonvbt==1 & c(1:NN)!=i
if (sum(availableneighborcases[availableneighborcases])==0) { print("Problem 3 Case 4") }
}

vecti <- c(X[i],Z1scalar[i,])
mattminusi <- cbind(Xbt[availableneighborcases],Z1scalarbt[availableneighborcases,])
if(length(availableneighborcases[availableneighborcases])==1) {
mattminusi <- matrix(c(Xbt[availableneighborcases],Z1scalarbt[availableneighborcases,]),
              ncol=length(vecti)) }
indsminusi <- c(1:NN)[availableneighborcases]
distsavailableneighborcases <- distsi(vecti,mattminusi)
# Identify the indices of the L smallest distances
Lnearestneighbors <- indsminusi[order(distsavailableneighborcases)][1:L]
# If length(indsminusi) < L, then NAs are introduced into Lnearestneighbors.  So only sample from the
# nonmissing elements of Lnearestneighbors
Lnearestneighbors <- Lnearestneighbors[!is.na(Lnearestneighbors)]
Vmat[i,m] <- sample(Vbt[Lnearestneighbors],size=1,replace=T)

cat(paste("Z2 not observed and Av not observed i:"),"\n")
cat(paste("Subject i X, Z1scalar = "),"\n")
print(vecti)
cat("\n")

cat(paste("Nearest neighbors X, Z1scalar = "),"\n")
print(mattminusi[order(distsavailableneighborcases)[1:L],])
cat(paste("\n"))

cat(paste("Farthest neighbors X, Z1scalar = "),"\n")
print(mattminusi[order(-distsavailableneighborcases)[1:L],])

}
# Done filling in Vmat for subject i

Zscoredistsnearest[i,] <- distsavailableneighborcases[order(distsavailableneighborcases)][1:L]

}

}

ans <- list(Vmat,Zscoredistsnearest)
return(ans)

}


# Test data

#set.seed(32)
## Need to enter Z1discrete, Z1scalar, Avscalar as matrices, even if only one element
## X, Delta, Z2, V are univariate vectors
#X <- runif(800)
#Delta <- rbinom(800,1,0.5)
#Z1discrete <- matrix(rbinom(1600,1,0.5),ncol=2)
#Z1scalar <- matrix(rnorm(1600),ncol=2)
#Z2 <- rnorm(800)
#epsilonz <- rbinom(800,1,0.3)
#Z2[epsilonz==0] <- NA
#V <- runif(800)
#V[Delta==0] <- NA
#epsilonv <- rep(1,800)
#epsilonv[Delta==1] <- rbinom(length(Delta[Delta==1]),1,0.7)
#V[epsilonv==0] <- NA
#Avdiscrete <- matrix(rep(NA,800),ncol=1)
#Avdiscrete <- matrix(rbinom(1600,1,.5),ncol=2)
#Avdiscrete[Delta==0,1] <- NA
#Avdiscrete[Delta==0,2] <- NA
#Avscalar <- matrix(rnorm(1600)+3,ncol=2)
#Avscalar[Delta==0,1] <- NA
#Avscalar[Delta==0,2] <- NA
#epsilona <- rep(1,800) # Case where auxiliaries of cases are always measured
#M <- 10
#L <- 5
#
#ans <- hotdeckMI(X,Delta,Z1discrete,Z1scalar,Z2,epsilonz,V,epsilonv,Avdiscrete,Avscalar,epsilona,M,L)
