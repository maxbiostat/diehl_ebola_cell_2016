getOdds <- function(mod, hier = TRUE){
  ## "inspired" by http://stackoverflow.com/questions/26417005/odds-ratio-and-confidence-intervals-from-glmer-output
  if(hier){
    cc <- confint(mod, parm = "beta_")  
    ctab <- cbind(est = fixef(mod), cc)
  }else{
    cc <- confint(mod)  
    ctab <- cbind(est = coef(mod), cc)
  }
  exp(ctab)
}
standz <- function(x) (x-mean(na.omit(x)))/(2*sd(na.omit(x)))
getSummary <- function(x, alpha = .95){
  return(c(
    lwr = quantile(x, probs = (1 - alpha)/2),
    mean = mean(x),
    upr = quantile(x, probs = (1 + alpha)/2)
  ))
}