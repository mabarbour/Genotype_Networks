pqlglmm.pois.model.offset <- function(y, offset.trans, groups, link, returnR=TRUE) {
  mod     <-  glmmPQL(y ~ 1 + offset(offset.trans),random=~1|groups,  family=quasipoisson(link=eval(link)), verbose=FALSE, control = lmeControl(opt = "optim")) 
  VarComp <- nlme::VarCorr(mod)
  beta0   <- as.numeric(mod$coefficients$fixed)
  omega   <- (as.numeric(VarComp[2,1]))
  var.a   <- (as.numeric(VarComp[1,1]))
  if (link=="log") {
    R.link  <- var.a/(var.a + omega*log(1/exp(beta0)+1))
    EY   	<- exp(beta0+var.a/2)
    R.org 	<- EY*(exp(var.a)-1)/(EY*(exp(var.a)-1)+omega) 
  }
  if (link=="sqrt") {
    R.link  <- var.a/(var.a + omega*0.25)
    R.org 	<- NA 
  }	
  if(returnR) return(list(R.link=R.link, R.org=R.org))
  else return(list(beta0=beta0, omega=omega, var.a=var.a))
}