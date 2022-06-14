#Utility functions for statistical analyses
# Author: Binyamin Knisbacher

#increment odds is an option to add 0.5 to each cell (aka Haldane-Anscombe correction)
runChisq <- function(a, a2=NULL, b=NULL, b2=NULL, fisher=F, alternative="two.sided", retVal=0, print=0, increment_odds=0){
  if(is.null(a2)){ #matrix input
    chisqMat = a
    a = chisqMat[1,1]; a2 = chisqMat[2,1]; b = chisqMat[1,2]; b2 = chisqMat[2,2];
  } else{ #single value input
    chisqMat = matrix(c(a, a2, b, b2), nrow=2)
  }
  if(!fisher){
    chisqRes = chisq.test(chisqMat, correct=F)
  } else {
    chisqRes = fisher.test(chisqMat, alternative=alternative)
  }

  oddsR = ((a+increment_odds) / (a2+increment_odds)) / ((b+increment_odds) / (b2+increment_odds))
  if(print){
    print(chisqMat)
    print(chisqRes)
    print(paste("OR =",oddsR))
  }
  if(retVal==0){
    return(chisqRes$p.value)
  } else if(retVal==2){
    return(c(chisqRes$p.value, oddsR))
  } else{
    return(c(a, a2, b, b2, chisqRes$p.value, oddsR))
  }

}


runChisq2 <- function(a, a2=NULL, b=NULL, b2=NULL){
  if(is.null(a2)){ #matrix input
    chisqMat = a
    a = chisqMat[1,1]; a2 = chisqMat[2,1]; b = chisqMat[1,2]; b2 = chisqMat[2,2];
  } else{ #single value input
    chisqMat = matrix(c(a, a2, b, b2), nrow=2)
  }
  chisqRes = chisq.test(chisqMat, correct=F)
  oddsR = (a / a2) / (b / b2)
  percentA = a/(a+a2)*100
  print(paste(a, a2, b, b2, chisqRes$p.value, oddsR, percentA))
}
