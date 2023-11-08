#
#Nov 2023
#


#
# F: fishing mortality 1/yr
# etaF: the coefficient determining the fish size with 50% fishing selectivity
setFishing = function(p, F, etaF=0.05 ) {
  
  for (iGroup in 1:p$nGroups) {
    ix = p$ix[[iGroup]]
    mFishing = etaF*max(p$mUpper[ix]) # selectivity at 0.05 of maximum size
    psi = ( 1 + (p$mc[ix]/mFishing)^(-3) )^(-1) # Standard trawl selectivity from Andersen (2019) Fig 5.2
    p$mortF[ix] = psi*F
  }
  return(p)
}