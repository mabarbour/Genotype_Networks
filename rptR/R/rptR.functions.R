#### Functions taken from rptR package

# general count model function
pqlglmm.pois.model <- function(y, # response variable
                               groups, # grouping variable
                               link, # link function
                               returnR=TRUE) {
  mod     <-  glmmPQL(y ~ 1,random=~1|groups,  family=quasipoisson(link=eval(link)), verbose=FALSE) 
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

# general count model with offset option. Most parameters same as for general count model
pqlglmm.pois.model.offset <- function(y, 
                                      offset.trans, # variable to use for offsetting that has already been transformed according to the link function used.
                                      groups, 
                                      link, 
                                      returnR=TRUE) {
  mod     <-  glmmPQL(y ~ 1 + offset(offset.trans),random=~1|groups,  family=quasipoisson(link=eval(link)), verbose=FALSE) 
  VarComp <- nlme::VarCorr(mod)
  beta0   <- as.numeric(mod$coefficients$fixed)
  omega   <- (as.numeric(VarComp[2,1]))
  var.a   <- (as.numeric(VarComp[1,1]))
  if (link=="log") {
    R.link  <- var.a/(var.a + omega*log(1/exp(beta0)+1))
    EY     <- exp(beta0+var.a/2)
    R.org 	<- EY*(exp(var.a)-1)/(EY*(exp(var.a)-1)+omega) 
  }
  if (link=="sqrt") {
    R.link  <- var.a/(var.a + omega*0.25)
    R.org 	<- NA 
  }	
  if(returnR) return(list(R.link=R.link, R.org=R.org))
  else return(list(beta0=beta0, omega=omega, var.a=var.a))
}

# confidence interval estimation by parametric bootstrapping
bootstr <- function(y, 
                    offset.trans,
                    groups, 
                    k, # number of levels in groups
                    N, # total sample size
                    beta0, # coefficient of model intercept
                    var.a, # Variance explained by grouping factor
                    omega, # Residual variance
                    link) { # appropriate link function for count model. 
  groupMeans <- rnorm(k, 0, sqrt(var.a))
  if(link=="log")  mu <- exp(beta0 + groupMeans[groups])
  if(link=="sqrt") mu <- (beta0 + groupMeans[groups])^2
  if (omega<=1)    y.boot <- rpois(N, mu)
  else         y.boot <- rnbinom(N, size=(mu/(omega-1)), mu=mu)
  pqlglmm.pois.model.offset(y.boot, offset.trans, groups, link) # no offset necessary because y.boot is based off the model with the offset already in it.
}

# significance test by randomization
permut   <- function(y, offset.trans, groups, N, link) {
  samp <- sample(1:N, N)
  pqlglmm.pois.model.offset(y, offset.trans, groups[samp], link)
}

