## Ensemble Code
## Exact Calculation of Freq Coverage Ratio
## Truncated Densities
wkdir <- "C:\\Users\\Certified Copy\\Workspace\\Reference Prior\\Code\\Case1"
setwd(wkdir)

## Load libraries
library(msm)
library(truncdist)
library(invgamma)
library(statmod)

## Constants
sigma <- 3
beta.tgam <- 1
beta.tinvgam <- 1
beta.tinvgauss <- 1
beta.tweibull  <- 5
xmin <- 1
sc   <-1
nn <- 20



## Truncated Normal
rf.prior.tnorm <- function(th){
  
		r.lm  <- th^2
		l.lm  <- th
		gr.lm <- 2*th
		gl.lm <- 1

   num2  <- dnorm(r.lm,th,sigma)*gr.lm
   num1  <- dnorm(l.lm,th,sigma)*gl.lm
   
   num <- num2 - num1
   denom <- pnorm(r.lm, th, sigma) - pnorm(l.lm, th, sigma)	
	
   lam  <- num/denom -th/sigma^2
 
   b1 <- abs(lam)*denom/num1
   b2 <- abs(lam)*denom/num2

   expr <- log(abs(lam))
   gam  <- -digamma(1)
   
   expr.add <- gam -b2*as.numeric(lam > 0) -b1*as.numeric(lam < 0)
   expr.add <- expr.add + b1 + b2	
   expr.add <- expr.add + (1/(b1-b2))*(b1*digamma(1/b1) - b2*digamma(1/b2) )  
 
   expr     <- expr + expr.add		
         
   expr <- exp(expr)  
	
expr
}

rf.prior.jnorm <- function(th){
expr <- 1

expr
}

rf.post.tnorm <- function(th,x){
  
   if(sc == 1){
		r.lm <- 3*th 
		l.lm <- 2*th
		gr.lm <- 3 
		gl.lm <- 2
	}
   
   if( sc == 2){   
		r.lm  <- th^(2)
		l.lm  <- th^(3/2)
		gr.lm <- 2*th
		gl.lm <- (3/2)*sqrt(th)
    }
	
    if( sc == 3){
		r.lm  <- th^2/2
		l.lm  <- th^2/3
		gr.lm <- th
		gl.lm <- (2/3)*th
     }
	
   if( sc == 4){
		r.lm  <- th^2
		l.lm  <- th
		gr.lm <- 2*th
		gl.lm <- 1
     }	

   num2  <- dnorm(r.lm,th,sigma)*gr.lm
   num1  <- dnorm(l.lm,th,sigma)*gl.lm
   
   num <- num2 - num1
   denom <- pnorm(r.lm, th, sigma) - pnorm(l.lm, th, sigma)	
	
   lam  <- num/denom -th/sigma^2
 
   b1 <- abs(lam)*denom/num1
   b2 <- abs(lam)*denom/num2

   expr <- log(abs(lam))
   gam  <- -digamma(1)
   
   expr.add <- gam -b2*as.numeric(lam > 0) -b1*as.numeric(lam < 0)
   expr.add <- expr.add + b1 + b2	
   expr.add <- expr.add + (1/(b1-b2))*(b1*digamma(1/b1) - b2*digamma(1/b2) )  
 
   expr     <- expr + expr.add		
    
   n <- length(x)
   m.vec      <- rep(th,n)
   l.lm.vec   <- rep(l.lm,n)
   r.lm.vec   <- rep(r.lm,n)
  
   expr <- expr + sum(dnorm(x, m.vec, sigma,log=TRUE)) 
   expr <- expr - sum(pnorm(r.lm.vec, m.vec, sigma, log = TRUE) - pnorm(l.lm.vec, m.vec, sigma , log = TRUE))	
     
   expr <- exp(expr)  
	
expr
}

rf.jeff.tnorm <- function(th,x){

   expr <- log(1/sigma) 
   

   if(sc == 1){
		r.lm <- 3*th 
		l.lm <- 2*th
		gr.lm <- 3 
		gl.lm <- 2
	}
   
   if( sc == 2){   
		r.lm  <- th^(2)
		l.lm  <- th^(3/2)
		gr.lm <- 2*th
		gl.lm <- (3/2)*sqrt(th)
    }
	
    if( sc == 3){
		r.lm  <- th^2/2
		l.lm  <- th^2/3
		gr.lm <- th
		gl.lm <- (2/3)*th
     }
	
   if( sc == 4){
		r.lm  <- th^2
		l.lm  <- th
		gr.lm <- 2*th
		gl.lm <- 1
     }
   n    <- length(x)
   m.vec      <- rep(th,n)
   l.lm.vec   <- rep(l.lm,n)
   r.lm.vec   <- rep(r.lm,n)
  
   expr <- expr + sum(dnorm(x, m.vec, sigma,log=TRUE))
   expr <- expr - sum(pnorm(r.lm.vec, m.vec, sigma, log = TRUE) - pnorm(l.lm.vec, m.vec, sigma , log = TRUE))	
   
   expr <- exp(expr)  
	
expr
}

rf.post.texp <- function(th,x){


   if(sc == 1){
		r.lm <- 3*th 
		l.lm <- 2*th
		gr.lm <- 3 
		gl.lm <- 2
	}
   
   if( sc == 2){   
		r.lm  <- th^(2)
		l.lm  <- th^(3/2)
		gr.lm <- 2*th
		gl.lm <- (3/2)*sqrt(th)
    }
	
    if( sc == 3){
		r.lm  <- th^2/2
		l.lm  <- th^2/3
		gr.lm <- th
		gl.lm <- (2/3)*th
     }
	
   if( sc == 4){
		r.lm  <- th^2
		l.lm  <- th
		gr.lm <- 2*th
		gl.lm <- 1
     }	
	
   num2  <- dexp(r.lm,rate=1/th)*(gr.lm)
   num1  <- dexp(l.lm,rate=1/th)*(gl.lm)	
   
   num <- num2 - num1
   
   denom <- pexp(r.lm, rate= 1/th) - pexp(l.lm, rate=1/th)	
	
   lam  <- num/denom + 1/th
 
   b1 <- abs(lam)*denom/num1
   b2 <- abs(lam)*denom/num2

   expr <- log(abs(lam))
   gam  <- -digamma(1)
   
   expr.add <- gam -b2*as.numeric(lam > 0) -b1*as.numeric(lam < 0)
   expr.add <- expr.add + b1 + b2	
   expr.add <- expr.add + (1/(b1-b2))*(b1*digamma(1/b1) - b2*digamma(1/b2) )  
 
   expr     <- expr + expr.add		
    
   n <- length(x)
   m.vec      <- rep(th,n)
   l.lm.vec   <- rep(l.lm,n)
   r.lm.vec   <- rep(r.lm,n)
  
   expr <- expr + sum(dexp(x,rate=1/m.vec,log=TRUE))
   expr <- expr - sum( pexp(r.lm.vec, rate = 1/m.vec, log = TRUE) - pexp(l.lm.vec, rate = 1/m.vec, log = TRUE))	
   
   expr <- exp(expr)
   
expr
}


rf.jeff.texp <- function(th,x){

   expr <- log(1/th)
       

   if(sc == 1){
		r.lm <- 3*th 
		l.lm <- 2*th
		gr.lm <- 3 
		gl.lm <- 2
	}
   
   if( sc == 2){   
		r.lm  <- th^(2)
		l.lm  <- th^(3/2)
		gr.lm <- 2*th
		gl.lm <- (3/2)*sqrt(th)
    }
	
    if( sc == 3){
		r.lm  <- th^2/2
		l.lm  <- th^2/3
		gr.lm <- th
		gl.lm <- (2/3)*th
     }
	
   if( sc == 4){
		r.lm  <- th^2
		l.lm  <- th
		gr.lm <- 2*th
		gl.lm <- 1
     }
   
   n <- length(x)
   n <- length(x)
   m.vec      <- rep(th,n)
   l.lm.vec   <- rep(l.lm,n)
   r.lm.vec   <- rep(r.lm,n)
  
   expr <- expr + sum(dexp(x,rate=1/m.vec,log=TRUE))
   expr <- expr - sum( pexp(r.lm.vec, rate = 1/m.vec, log = TRUE) - pexp(l.lm.vec, rate = 1/m.vec, log = TRUE))
   #expr.norm <- log(integrate(Vectorize(rf.j.texp),l.lm,r.lm, subdivisions=5000, rel.tol=.Machine$double.eps^.05)$value)
   #expr <- expr - expr.norm
   expr <- exp(expr)	
	
expr
}


rf.post.tgam <- function(th,x){
   

   if(sc == 1){
		r.lm <- 3*th 
		l.lm <- 2*th
		gr.lm <- 3 
		gl.lm <- 2
	}
   
   if( sc == 2){   
		r.lm  <- th^(2)
		l.lm  <- th^(3/2)
		gr.lm <- 2*th
		gl.lm <- (3/2)*sqrt(th)
    }
	
    if( sc == 3){
		r.lm  <- th^2/2
		l.lm  <- th^2/3
		gr.lm <- th
		gl.lm <- (2/3)*th
     }
	
   if( sc == 4){
		r.lm  <- th^2
		l.lm  <- th
		gr.lm <- 2*th
		gl.lm <- 1
     }		

   num2 <- dgamma(r.lm, shape = th, scale = beta.tgam)*gr.lm
   num1 <- dgamma(l.lm, shape = th, scale = beta.tgam)*gl.lm
   
   num <- num2 - num1
   
   denom <- pgamma(r.lm, shape= th, scale = beta.tgam) - pgamma(l.lm, shape=th, scale = beta.tgam)	
   lam  <- num/denom -((th - 1)*trigamma(th)-1)
 
   b1 <- abs(lam)*denom/num1
   b2 <- abs(lam)*denom/num2

   expr <- log(abs(lam))
   gam  <- -digamma(1)
   
   expr.add <- gam -b2*as.numeric(lam > 0) -b1*as.numeric(lam < 0)
   expr.add <- expr.add + b1 + b2	
   expr.add <- expr.add + (1/(b1-b2))*(b1*digamma(1/b1) - b2*digamma(1/b2) )  
 
   expr     <- expr + expr.add		
    
   n <- length(x)
   m.vec      <- rep(th,n)
   l.lm.vec   <- rep(l.lm,n)
   r.lm.vec   <- rep(r.lm,n)
  
   expr <- expr + sum(dgamma(x,shape=m.vec,scale= beta.tgam,log=TRUE))
   expr <- expr - sum(pgamma(r.lm.vec, shape = m.vec, scale = beta.tgam,log = TRUE) - pgamma(l.lm.vec, shape =m.vec, scale = beta.tgam, log = TRUE))	

   expr <- exp(expr)
   
expr
}

rf.jeff.tgam <- function(th,x){

   expr <- .5*log( trigamma(th) ) 


   if(sc == 1){
		r.lm <- 3*th 
		l.lm <- 2*th
		gr.lm <- 3 
		gl.lm <- 2
	}
   
   if( sc == 2){   
		r.lm  <- th^(2)
		l.lm  <- th^(3/2)
		gr.lm <- 2*th
		gl.lm <- (3/2)*sqrt(th)
    }
	
    if( sc == 3){
		r.lm  <- th^2/2
		l.lm  <- th^2/3
		gr.lm <- th
		gl.lm <- (2/3)*th
     }
	
   if( sc == 4){
		r.lm  <- th^2
		l.lm  <- th
		gr.lm <- 2*th
		gl.lm <- 1
     }
   
   n <- length(x)
   m.vec      <- rep(th,n)
   l.lm.vec   <- rep(l.lm,n)
   r.lm.vec   <- rep(r.lm,n)
  
   expr <- expr + sum(dgamma(x,shape=m.vec,scale= beta.tgam,log=TRUE))
   expr <- expr - sum(pgamma(r.lm.vec, shape = m.vec, scale = beta.tgam,log = TRUE) - pgamma(l.lm.vec, shape =m.vec, scale = beta.tgam, log = TRUE))	
   
   expr <- exp(expr)
	
expr
}
	
rf.post.tlnorm <- function(th,x){
   

   if(sc == 1){
		r.lm <- 3*th 
		l.lm <- 2*th
		gr.lm <- 3 
		gl.lm <- 2
	}
   
   if( sc == 2){   
		r.lm  <- th^(2)
		l.lm  <- th^(3/2)
		gr.lm <- 2*th
		gl.lm <- (3/2)*sqrt(th)
    }
	
    if( sc == 3){
		r.lm  <- th^2/2
		l.lm  <- th^2/3
		gr.lm <- th
		gl.lm <- (2/3)*th
     }
	
   if( sc == 4){
		r.lm  <- th^2
		l.lm  <- th
		gr.lm <- 2*th
		gl.lm <- 1
     }	

   num2  <- dlnorm(r.lm,th,sigma)*gr.lm
   num1  <- dlnorm(l.lm,th,sigma)*gl.lm
   
   num <- num2 - num1
   denom <- plnorm(r.lm, th, sigma) - plnorm(l.lm, th, sigma)	
	
   lam  <- num/denom
 
   b1 <- abs(lam)*denom/num1
   b2 <- abs(lam)*denom/num2

   expr <- log(abs(lam))
   gam  <- -digamma(1)
   
   expr.add <- gam -b2*as.numeric(lam > 0) -b1*as.numeric(lam < 0)
   expr.add <- expr.add + b1 + b2	
   expr.add <- expr.add + (1/(b1-b2))*(b1*digamma(1/b1) - b2*digamma(1/b2) )  
 
   expr     <- expr + expr.add		
    
   n <- length(x)
   m.vec      <- rep(th,n)
   l.lm.vec   <- rep(l.lm,n)
   r.lm.vec   <- rep(r.lm,n)
  
   expr <- expr + sum(dlnorm(x, m.vec, sigma,log=TRUE))
   expr <- expr - sum(plnorm(r.lm.vec, m.vec, sigma, log = TRUE) - plnorm(l.lm.vec, m.vec, sigma , log = TRUE))	
   
   expr <- exp(expr)
	
expr
}

rf.jeff.tlnorm <- function(th,x){

   expr <- log(1/sigma) 
      
	  

   if(sc == 1){
		r.lm <- 3*th 
		l.lm <- 2*th
		gr.lm <- 3 
		gl.lm <- 2
	}
   
   if( sc == 2){   
		r.lm  <- th^(2)
		l.lm  <- th^(3/2)
		gr.lm <- 2*th
		gl.lm <- (3/2)*sqrt(th)
    }
	
    if( sc == 3){
		r.lm  <- th^2/2
		l.lm  <- th^2/3
		gr.lm <- th
		gl.lm <- (2/3)*th
     }
	
   if( sc == 4){
		r.lm  <- th^2
		l.lm  <- th
		gr.lm <- 2*th
		gl.lm <- 1
     }

   n <- length(x)
   m.vec      <- rep(th,n)
   l.lm.vec   <- rep(l.lm,n)
   r.lm.vec   <- rep(r.lm,n)
  
   expr <- expr + sum(dlnorm(x, m.vec, sigma,log=TRUE))
   expr <- expr - sum(plnorm(r.lm.vec, m.vec, sigma, log = TRUE) - plnorm(l.lm.vec, m.vec, sigma , log = TRUE))	
   
   expr <- exp(expr)   
   
expr
}

rf.post.tinvgam <- function(th,x){
   

   if(sc == 1){
		r.lm <- 3*th 
		l.lm <- 2*th
		gr.lm <- 3 
		gl.lm <- 2
	}
   
   if( sc == 2){   
		r.lm  <- th^(2)
		l.lm  <- th^(3/2)
		gr.lm <- 2*th
		gl.lm <- (3/2)*sqrt(th)
    }
	
    if( sc == 3){
		r.lm  <- th^2/2
		l.lm  <- th^2/3
		gr.lm <- th
		gl.lm <- (2/3)*th
     }
	
   if( sc == 4){
		r.lm  <- th^2
		l.lm  <- th
		gr.lm <- 2*th
		gl.lm <- 1
     }		

   num2 <- dinvgamma(1/l.lm, shape = th, scale = beta.tinvgam)*(-gl.lm/l.lm^2)
   num1 <- dinvgamma(1/r.lm, shape = th, scale = beta.tinvgam)*(-gr.lm/r.lm^2)
	
   num <-  num2 - num1
   
   denom <- pinvgamma(r.lm, shape= th, scale = beta.tinvgam) - pinvgamma(l.lm, shape=th, scale = beta.tinvgam)	
	
   lam  <- num/denom - ((th + 1)*trigamma(th)-1)
 
   b1 <- abs(lam)*denom/num1
   b2 <- abs(lam)*denom/num2

   expr <- log(abs(lam))
   gam  <- -digamma(1)
   
   expr.add <- gam -b2*as.numeric(lam > 0) -b1*as.numeric(lam < 0)
   expr.add <- expr.add + b1 + b2	
   expr.add <- expr.add + (1/(b1-b2))*(b1*digamma(1/b1) - b2*digamma(1/b2) )  
 
   expr     <- expr + expr.add		
    
   n <- length(x)
   m.vec      <- rep(th,n)
   l.lm.vec   <- rep(l.lm,n)
   r.lm.vec   <- rep(r.lm,n)
  
   expr <- expr + sum(dinvgamma(x,shape=m.vec,scale= beta.tinvgam,log=TRUE)) 
   expr <- expr - sum(pinvgamma(r.lm.vec, shape = m.vec, scale = beta.tinvgam,log = TRUE) - pinvgamma(l.lm.vec, shape =m.vec, scale = beta.tinvgam, log = TRUE))	


   expr <- exp(expr)
   #expr <- exp(expr)	
   
expr
}

rf.jeff.tinvgam <- function(th,x){

   expr <- .5*log( trigamma(th) ) 
    

   if(sc == 1){
		r.lm <- 3*th 
		l.lm <- 2*th
		gr.lm <- 3 
		gl.lm <- 2
	}
   
   if( sc == 2){   
		r.lm  <- th^(2)
		l.lm  <- th^(3/2)
		gr.lm <- 2*th
		gl.lm <- (3/2)*sqrt(th)
    }
	
    if( sc == 3){
		r.lm  <- th^2/2
		l.lm  <- th^2/3
		gr.lm <- th
		gl.lm <- (2/3)*th
     }
	
   if( sc == 4){
		r.lm  <- th^2
		l.lm  <- th
		gr.lm <- 2*th
		gl.lm <- 1
     }	


   n <- length(x)
   n <- length(x)
   m.vec      <- rep(th,n)
   l.lm.vec   <- rep(l.lm,n)
   r.lm.vec   <- rep(r.lm,n)
  
   expr <- expr + sum(dinvgamma(x,shape=m.vec,scale= beta.tinvgam,log=TRUE))
   expr <- expr - sum(pinvgamma(r.lm.vec, shape = m.vec, scale = beta.tinvgam,log = TRUE) - pinvgamma(l.lm.vec, shape =m.vec, scale = beta.tinvgam, log = TRUE))	
   
   expr <- exp(expr)	
      
expr
}

rf.post.tinvgauss <- function(th,x){
   

   if(sc == 1){
		r.lm <- 3*th 
		l.lm <- 2*th
		gr.lm <- 3 
		gl.lm <- 2
	}
   
   if( sc == 2){   
		r.lm  <- th^(2)
		l.lm  <- th^(3/2)
		gr.lm <- 2*th
		gl.lm <- (3/2)*sqrt(th)
    }
	
    if( sc == 3){
		r.lm  <- th^2/2
		l.lm  <- th^2/3
		gr.lm <- th
		gl.lm <- (2/3)*th
     }
	
   if( sc == 4){
		r.lm  <- th^2
		l.lm  <- th
		gr.lm <- 2*th
		gl.lm <- 1
     }	

   num2 <- dinvgauss(1/l.lm, th, beta.tinvgam)*(-gl.lm/l.lm^2)
   num1 <- dinvgauss(1/r.lm, th, beta.tinvgam)*(-gr.lm/r.lm^2)
	
   num <-  num2 - num1
   
  
   denom <- pinvgauss(r.lm,th,beta.tinvgauss) - pinvgauss(l.lm,th,beta.tinvgauss)	
	
   lam  <- num/denom
 
   b1 <- abs(lam)*denom/num1
   b2 <- abs(lam)*denom/num2

   expr <- log(abs(lam))
   gam  <- -digamma(1)
   
   expr.add <- gam -b2*as.numeric(lam > 0) -b1*as.numeric(lam < 0)
   expr.add <- expr.add + b1 + b2	
   expr.add <- expr.add + (1/(b1-b2))*(b1*digamma(1/b1) - b2*digamma(1/b2) )  
 
   expr     <- expr + expr.add		
    
   n <- length(x)
   m.vec      <- rep(th,n)
   l.lm.vec   <- rep(l.lm,n)
   r.lm.vec   <- rep(r.lm,n)
  
   expr <- expr + sum(dinvgauss(x, m.vec, beta.tinvgauss,log=TRUE))
   expr <- expr - sum(pinvgauss(r.lm.vec, m.vec, beta.tinvgauss,log = TRUE) - pinvgauss(l.lm.vec, m.vec,beta.tinvgauss, log = TRUE))	
  
   expr <- exp(expr)
	
expr
}

rf.jeff.tinvgauss <- function(th,x){

   expr <- .5*log(max( th^{-3}*beta.tinvgauss -2,0.001) )
      

   if(sc == 1){
		r.lm <- 3*th 
		l.lm <- 2*th
		gr.lm <- 3 
		gl.lm <- 2
	}
   
   if( sc == 2){   
		r.lm  <- th^(2)
		l.lm  <- th^(3/2)
		gr.lm <- 2*th
		gl.lm <- (3/2)*sqrt(th)
    }
	
    if( sc == 3){
		r.lm  <- th^2/2
		l.lm  <- th^2/3
		gr.lm <- th
		gl.lm <- (2/3)*th
     }
	
   if( sc == 4){
		r.lm  <- th^2
		l.lm  <- th
		gr.lm <- 2*th
		gl.lm <- 1
     }		

   n <- length(x)
   n <- length(x)
   m.vec      <- rep(th,n)
   l.lm.vec   <- rep(l.lm,n)
   r.lm.vec   <- rep(r.lm,n)
  
   expr <- expr + sum(dinvgauss(x, m.vec, beta.tinvgauss,log=TRUE))
   expr <- expr - sum(pinvgauss(r.lm.vec, m.vec, beta.tinvgauss,log = TRUE) - pinvgauss(l.lm.vec, m.vec,beta.tinvgauss, log = TRUE))	
  
   expr <- exp(expr)
  	
expr
}


rf.post.tweibull <- function(th,x){
  

   if(sc == 1){
		r.lm <- 3*th 
		l.lm <- 2*th
		gr.lm <- 3 
		gl.lm <- 2
	}
   
   if( sc == 2){   
		r.lm  <- th^(2)
		l.lm  <- th^(3/2)
		gr.lm <- 2*th
		gl.lm <- (3/2)*sqrt(th)
    }
	
    if( sc == 3){
		r.lm  <- th^2/2
		l.lm  <- th^2/3
		gr.lm <- th
		gl.lm <- (2/3)*th
     }
	
   if( sc == 4){
		r.lm  <- th^2
		l.lm  <- th
		gr.lm <- 2*th
		gl.lm <- 1
     }	
   
   num2 <- dweibull(r.lm, shape = beta.tweibull, scale = th)*gr.lm
   num1 <- dweibull(l.lm, shape = beta.tweibull, scale = th)*gl.lm   
   
   num <- num2 - num1
   
   denom <- pweibull(r.lm, shape= beta.tweibull, scale = th) - pweibull(l.lm, shape=beta.tweibull, scale = th)	
	
   lam  <- num/denom -beta.tweibull/th^(2*beta.tweibull + 2)
 
   b1 <- abs(lam)*denom/num1
   b2 <- abs(lam)*denom/num2

   expr <- log(abs(lam))
   gam  <- -digamma(1)
   
   expr.add <- gam -b2*as.numeric(lam > 0) -b1*as.numeric(lam < 0)
   expr.add <- expr.add + b1 + b2	
   expr.add <- expr.add + (1/(b1-b2))*(b1*digamma(1/b1) - b2*digamma(1/b2) )  
 
   expr     <- expr + expr.add		
    
   n <- length(x)
   m.vec      <- rep(th,n)
   l.lm.vec   <- rep(l.lm,n)
   r.lm.vec   <- rep(r.lm,n)
  
   expr <- expr + sum(dweibull(x,shape=beta.tweibull,scale=m.vec,log=TRUE))
   expr <- expr - sum(pweibull(r.lm.vec, shape= beta.tweibull, scale = m.vec,log = TRUE) - pweibull(l.lm.vec,shape=beta.tweibull, scale =m.vec,log = TRUE))	
   
   expr <- exp(expr)  
	
expr
}


rf.jeff.tweibull <- function(th,x){

   expr <- sqrt( beta.tweibull/th^2 + beta.tweibull*(beta.tweibull+1)/th^(2*(beta.tweibull+1) ) ) 
   

   if(sc == 1){
		r.lm <- 3*th 
		l.lm <- 2*th
		gr.lm <- 3 
		gl.lm <- 2
	}
   
   if( sc == 2){   
		r.lm  <- th^(2)
		l.lm  <- th^(3/2)
		gr.lm <- 2*th
		gl.lm <- (3/2)*sqrt(th)
    }
	
    if( sc == 3){
		r.lm  <- th^2/2
		l.lm  <- th^2/3
		gr.lm <- th
		gl.lm <- (2/3)*th
     }
	
   if( sc == 4){
		r.lm  <- th^2
		l.lm  <- th
		gr.lm <- 2*th
		gl.lm <- 1
     }
  
  n <- length(x)
   n <- length(x)
   m.vec      <- rep(th,n)
   l.lm.vec   <- rep(l.lm,n)
   r.lm.vec   <- rep(r.lm,n)
  
   expr <- expr + sum(dweibull(x,shape=beta.tweibull,scale=m.vec,log=TRUE))
   expr <- expr - sum(pweibull(r.lm.vec, shape= beta.tweibull, scale = m.vec,log = TRUE) - pweibull(l.lm.vec,shape=beta.tweibull, scale =m.vec,log = TRUE))	
   
   expr <- exp(expr)	
	
expr
}
 
#####################################################################
## Secondary Function
## Mean Absolute Deviation
## Posterior Comparison
#####################################################################
#Parallel Computation
###################################


tcomp <- function(r){
#######################################################################
library(msm)
library(truncdist)
library(invgamma)
library(statmod)

m1 <- 1000
m2 <- 1000

## Constants
sigma <- 3
beta.tgam <- 1
beta.tinvgam <- 1
beta.tinvgauss <- 1
beta.tweibull  <- 2
xmin <- 1
nn <- 20


th.vec <- c(5)
n.vec  <- c(3, 10, 30)

th.est1 <- th.est2 <- th.est3 <- rep(NA, m1 + 1)
th.est1.m <- th.est2.m <- the.est3.m <- rep(NA,m1/2)
rho1 <- rho2 <- rho3  <- rep(NA, m2)
rho <- matrix(NA,m2,6)
m.1 <- m.2 <- m.3 <-m.4 <- rep(NA, m2)

th.est   <- matrix(NA,m1/2,6)
md<- mse      <- matrix(NA,9,6)
th.hist  <- Out.result <- matrix(NA,m2,54)

#Start of new counter
md <-mse<-th.hist<- Out.result <-c()
count <- 0
	
for( k in 1:  length(n.vec) ){
		for( l in 1: length(th.vec) ){
			count <- count + 1
			#Start of new scenario
			m<-ms<-th.est <- rho<-c()
	
			for(sc in c(4)){
	
				for(i in 1:m2){

#################################
####Simulating X
#################################
	n <- n.vec[k]
	theta_0 <- th.vec[l]
   	spec.list <- c("norm","exp","gamma","lnorm", "invgamma","invgauss","weibull")
	spec <- spec.list[r]	
	if( r== 1){

   if(sc == 1){
		r.lm <- 3*theta_0 
		l.lm <- 2*theta_0
		gr.lm <- 3 
		gl.lm <- 2
	}
   
   if( sc == 2){   
		r.lm  <- theta_0^(2)
		l.lm  <- theta_0^(3/2)
		gr.lm <- 2*theta_0
		gl.lm <- (3/2)*sqrt(theta_0)
    }
	
    if( sc == 3){
		r.lm  <- theta_0^2/2
		l.lm  <- theta_0^2/3
		gr.lm <- theta_0
		gl.lm <- (2/3)*theta_0
     }
	
   if( sc == 4){
		r.lm  <- theta_0^2
		l.lm  <- theta_0
		gr.lm <- 2*theta_0
		gl.lm <- 1
     }
		
		x <- rtnorm(n,theta_0,sigma, l.lm, r.lm )
		if(sc == 1){
			low.lim   <-  ((max(x))/3)
			high.lim  <-  (min(x)/2)
		}
		
		if(sc == 2){
			low.lim   <-  ( max(x) )^(1/2)
			high.lim  <-  (min(x))^(2/3)
		}
		
		if(sc == 3){
			low.lim   <-  sqrt(2*(max(x)))
			high.lim  <-  sqrt(3*min(x))
		}
		
		if(sc == 4){
			low.lim   <-  sqrt(max(x))
			high.lim  <-  (min(x))
		}
		
		}
				
	if( r== 2){	
   if(sc == 1){
		r.lm <- 3*theta_0 
		l.lm <- 2*theta_0
		gr.lm <- 3 
		gl.lm <- 2
	}
   
   if( sc == 2){   
		r.lm  <- theta_0^(2)
		l.lm  <- theta_0^(3/2)
		gr.lm <- 2*theta_0
		gl.lm <- (3/2)*sqrt(theta_0)
    }
	
    if( sc == 3){
		r.lm  <- theta_0^2/2
		l.lm  <- theta_0^2/3
		gr.lm <- theta_0
		gl.lm <- (2/3)*theta_0
     }
	
   if( sc == 4){
		r.lm  <- theta_0^2
		l.lm  <- theta_0
		gr.lm <- 2*theta_0
		gl.lm <- 1
     }				
		x <- rtrunc(n,spec="exp",rate=1/theta_0,l.lm,r.lm)
		if(sc == 1){
			low.lim   <-  ((max(x))/3)
			high.lim  <-  (min(x)/2)
		}
		
		if(sc == 2){
			low.lim   <-  ( max(x) )^(1/2)
			high.lim  <-  (min(x))^(2/3)
		}
		
		if(sc == 3){
			low.lim   <-  sqrt(2*(max(x)))
			high.lim  <-  sqrt(3*min(x))
		}
		
		if(sc == 4){
			low.lim   <-  sqrt(max(x))
			high.lim  <-  (min(x))
		}
	}
				
	if( r== 3){
	
	   if(sc == 1){
		r.lm <- 3*theta_0 
		l.lm <- 2*theta_0
		gr.lm <- 3 
		gl.lm <- 2
	}
   
   if( sc == 2){   
		r.lm  <- theta_0^(2)
		l.lm  <- theta_0^(3/2)
		gr.lm <- 2*theta_0
		gl.lm <- (3/2)*sqrt(theta_0)
    }
	
    if( sc == 3){
		r.lm  <- theta_0^2/2
		l.lm  <- theta_0^2/3
		gr.lm <- theta_0
		gl.lm <- (2/3)*theta_0
     }
	
   if( sc == 4){
		r.lm  <- theta_0^2
		l.lm  <- theta_0
		gr.lm <- 2*theta_0
		gl.lm <- 1
     }
		x <- rtrunc(n,spec="gamma",shape=theta_0,scale=beta.tgam, a=l.lm, b=r.lm )
					
		if(sc == 1){
			low.lim   <-  ((max(x))/3)
			high.lim  <-  (min(x)/2)
		}
		
		if(sc == 2){
			low.lim   <-  ( max(x) )^(1/2)
			high.lim  <-  (min(x))^(2/3)
		}
		
		if(sc == 3){
			low.lim   <-  sqrt(2*(max(x)))
			high.lim  <-  sqrt(3*min(x))
		}
		
		if(sc == 4){
			low.lim   <-  sqrt(max(x))
			high.lim  <-  (min(x))
		}
	}
				
	if( r== 4){
	
   if(sc == 1){
		r.lm <- 3*theta_0 
		l.lm <- 2*theta_0
		gr.lm <- 3 
		gl.lm <- 2
	}
   
   if( sc == 2){   
		r.lm  <- theta_0^(2)
		l.lm  <- theta_0^(3/2)
		gr.lm <- 2*theta_0
		gl.lm <- (3/2)*sqrt(theta_0)
    }
	
    if( sc == 3){
		r.lm  <- theta_0^2/2
		l.lm  <- theta_0^2/3
		gr.lm <- theta_0
		gl.lm <- (2/3)*theta_0
     }
	
   if( sc == 4){
		r.lm  <- theta_0^2
		l.lm  <- theta_0
		gr.lm <- 2*theta_0
		gl.lm <- 1
     }
		x  <- rtrunc(n, spec="lnorm",meanlog=theta_0,sdlog=sigma, a=l.lm,b= r.lm )
		if(sc == 1){
			low.lim   <-  ((max(x))/3)
			high.lim  <-  (min(x)/2)
		}
		
		if(sc == 2){
			low.lim   <-  ( max(x) )^(1/2)
			high.lim  <-  (min(x))^(2/3)
		}
		
		if(sc == 3){
			low.lim   <-  sqrt(2*(max(x)))
			high.lim  <-  sqrt(3*min(x))
		}
		
		if(sc == 4){
			low.lim   <-  sqrt(max(x))
			high.lim  <-  (min(x))
		}
	}
				
	if ( r==5){
	
   if(sc == 1){
		r.lm <- 3*theta_0 
		l.lm <- 2*theta_0
		gr.lm <- 3 
		gl.lm <- 2
	}
   
   if( sc == 2){   
		r.lm  <- theta_0^(2)
		l.lm  <- theta_0^(3/2)
		gr.lm <- 2*theta_0
		gl.lm <- (3/2)*sqrt(theta_0)
    }
	
    if( sc == 3){
		r.lm  <- theta_0^2/2
		l.lm  <- theta_0^2/3
		gr.lm <- theta_0
		gl.lm <- (2/3)*theta_0
     }
	
   if( sc == 4){
		r.lm  <- theta_0^2
		l.lm  <- theta_0
		gr.lm <- 2*theta_0
		gl.lm <- 1
     }		
		x   <- rtrunc(n,spec="invgamma",shape=theta_0,scale=beta.tinvgam, a=l.lm, b=r.lm )
		if(sc == 1){
			low.lim   <-  ((max(x))/3)
			high.lim  <-  (min(x)/2)
		}
		
		if(sc == 2){
			low.lim   <-  ( max(x) )^(1/2)
			high.lim  <-  (min(x))^(2/3)
		}
		
		if(sc == 3){
			low.lim   <-  sqrt(2*(max(x)))
			high.lim  <-  sqrt(3*min(x))
		}
		
		if(sc == 4){
			low.lim   <-  sqrt(max(x))
			high.lim  <-  (min(x))
		}
	}
				
	if( r == 6){
	
   if(sc == 1){
		r.lm <- 3*theta_0 
		l.lm <- 2*theta_0
		gr.lm <- 3 
		gl.lm <- 2
	}
   
   if( sc == 2){   
		r.lm  <- theta_0^(2)
		l.lm  <- theta_0^(3/2)
		gr.lm <- 2*theta_0
		gl.lm <- (3/2)*sqrt(theta_0)
    }
	
    if( sc == 3){
		r.lm  <- theta_0^2/2
		l.lm  <- theta_0^2/3
		gr.lm <- theta_0
		gl.lm <- (2/3)*theta_0
     }
	
   if( sc == 4){
		r.lm  <- theta_0^2
		l.lm  <- theta_0
		gr.lm <- 2*theta_0
		gl.lm <- 1
     }		
		x   <- rtrunc(n,spec="invgauss",theta_0,beta.tinvgauss, a=l.lm, b=r.lm )
		
		if(sc == 1){
			low.lim   <-  ((max(x))/3)
			high.lim  <-  (min(x)/2)
		}
		
		if(sc == 2){
			low.lim   <-  ( max(x) )^(1/2)
			high.lim  <-  (min(x))^(2/3)
		}
		
		if(sc == 3){
			low.lim   <-  sqrt(2*(max(x)))
			high.lim  <-  sqrt(3*min(x))
		}
		
		if(sc == 4){
			low.lim   <-  sqrt(max(x))
			high.lim  <-  (min(x))
		}
	}
				
	if( r == 7){
	
   if(sc == 1){
		r.lm <- 3*theta_0 
		l.lm <- 2*theta_0
		gr.lm <- 3 
		gl.lm <- 2
	}
   
   if( sc == 2){   
		r.lm  <- theta_0^(2)
		l.lm  <- theta_0^(3/2)
		gr.lm <- 2*theta_0
		gl.lm <- (3/2)*sqrt(theta_0)
    }
	
    if( sc == 3){
		r.lm  <- theta_0^2/2
		l.lm  <- theta_0^2/3
		gr.lm <- theta_0
		gl.lm <- (2/3)*theta_0
     }
	
   if( sc == 4){
		r.lm  <- theta_0^2
		l.lm  <- theta_0
		gr.lm <- 2*theta_0
		gl.lm <- 1
     }		
		x     <- rtrunc(n,spec="weibull",shape=beta.tweibull,scale=theta_0, a=l.lm, b=r.lm )
		if(sc == 1){
			low.lim   <-  ((max(x))/3)
			high.lim  <-  (min(x)/2)
		}
		
		if(sc == 2){
			low.lim   <-  ( max(x) )^(1/2)
			high.lim  <-  (min(x))^(2/3)
		}
		
		if(sc == 3){
			low.lim   <-  sqrt(2*(max(x)))
			high.lim  <-  sqrt(3*min(x))
		}
		
		if(sc == 4){
			low.lim   <-  sqrt(max(x))
			high.lim  <-  (min(x))
		}				
	}
    

   th.est1[1] <- th.est2[1] <- th.est3[1]<- .5*(low.lim + high.lim)
   th1_0 <- th2_0 <- th3_0	<- th.est1[1]
		

		##Normal
		if( r== 1){
		counter <- 0
		M       <- optim(runif(1,low.lim,high.lim),rf.post.tnorm,method="L-BFGS-B",lower=low.lim,upper=high.lim,x=x)$value
		M       <- max(M,rf.post.tnorm(low.lim,x), rf.post.tnorm(high.lim,x))

			repeat{
				th 	<- runif(1,low.lim , high.lim )
				rat     <- exp(log(rf.post.tnorm(th,x)))/M
				
				u <- runif(1, 0, 1)			
			
				if( u < rat &!is.na(rat)){
					
					counter <- counter + 1
					th.est1[counter ]	  <- th
					th1_0             <- th
				}
				
				if(counter > m1){
					break
				}	
			}
		}

		
	
	
	rho1[i]    <- integrate(rf.post.tnorm,low.lim, theta_0,x)$value/integrate(rf.post.tnorm,low.lim, high.lim,x)$value
	#rho2[i]    <- integrate(Vectorize(rf.jeff.tnorm),low.lim, theta_0,x)$value/integrate(Vectorize(rf.jeff.tnorm),low.lim, high.lim,x)$value

		if(i %% 100 == 0 ){
			print( paste(Sys.time(), "Truncated Model, Iteration number:",i,k,l) )			
		}

}#End of loopi

	

	
	rho <- cbind(rho,rho1,rho2)
	th.est <- cbind(th.est, th.est1.m, th.est2.m)
	m <- c(m, mean(m.1),mean(m.2))
	
	ms <- c(ms, sqrt(mean(m.3)),sqrt(mean(m.4)))
	
    #th.est[, ll:hh]  <- cbind(th.est1.m,th.est2.m)
	#rho[,ll:hh]     <- cbind(rho1, rho2)
   	#md[count,ll:hh] <- c( mean(m.1), mean(m.2))
    #mse[count,ll:hh] <- c( mean(m.3), mean(m.4))

	#ll <- (count-1)*6+1
	#hh <- (count)*6
	
	
	#th.hist[,ll:hh]   <- th.est
	#Out.result[,ll:hh] <- rho	
}#End of loopsc
	Out.result <- cbind(Out.result,rho)
	th.hist <- cbind(th.hist,th.est)
    mse <- rbind(mse,ms)
    md <- rbind(md,m)


}#End of loopth

}
print(r)

Out <- list(Out.result, th.hist, md, mse)
#Out <- rho

Out
}#End of function


res1 <- tcomp(1)

res2 <- tcomp(2)

res3 <- tcomp(3)

res4 <- tcomp(4)

res5 <- tcomp(5)

r6 <- tcomp(6)



res2 <- lapply(c(5,6,7),tcomp)

res3 <- lapply(c(1,2,3),tcomp)


# Frequency Cov Comparison
fc.alpha <- function(alpha){

o <- length( which( res1[[1]][,j] <= alpha ) )/m2

o
}

i<- 0


par(mfrow=c(2,2))

for(c in 1:3){
	
	g <- th.vec[c]

	j <- i + (2*c-1) 
	f <- Vectorize(fc.alpha)
	curve(f,0, 1, lty=1,lwd=1, sub=
	bquote("True Value"~theta==.(g)),xlab=expression(alpha),ylab="Coverage",ylim=c(0,1),col="red",cex.lab=1.5)

	f <- Vectorize(fc.alpha)
	curve(f,0, 1, lty=2,lwd=1,xlab=expression(alpha),ylab="Frequency Coverage",add=TRUE,col="red")

	
	j <- i + (2*c-1) + 6
	f <- Vectorize(fc.alpha)
	curve(f,0, 1, lty=2,lwd=1,xlab=expression(alpha),ylab="Frequency Coverage",add=TRUE,col="red")


	j <- i + (2*c-1) + 12
	f <- Vectorize(fc.alpha)
	curve(f,0, 1, lty=3,lwd=1,xlab=expression(alpha),ylab="Frequency Coverage",add=TRUE,col="red")
	abline(0,1,col="blue")
		

}

fc.m <- matrix(NA,27,2)
for(c in 1:27){

j <- 2*(c-1) + 1
fc.m[c,1] <-fc.alpha(0.95)

j <- 2*c
fc.m[c,2] <-fc.alpha(0.95)

}

#######################################################################
library(msm)
library(truncdist)
library(invgamma)
library(statmod)

local_fn <-c("rf.jeff.texp","rf.jeff.tgam","rf.jeff.tinvgam","rf.jeff.tinvgauss","rf.jeff.tlnorm","rf.jeff.tnorm","rf.jeff.tweibull","rf.post.texp"     
 ,"rf.post.tgam","rf.post.tinvgam","rf.post.tinvgauss","rf.post.tlnorm","rf.post.tnorm","rf.post.tweibull")
 
no_cores <- detectCores()
cl <- makeCluster(no_cores)
registerDoSNOW(cl)

results<- foreach(r=1:7) %dopar% {
 tcomp(r)
 }
 
stopCluster(cl)
####################

#results<-system.time(parLapply(cl,1:7,tcomp))



#######################################################

fc.alpha <- function(alpha){
		
	o <- length( which( res1[[1]][,j] <= alpha ) )/m2

o
}


par(oma=c(0,3,0,0),mfrow=c(2,2))


	j  <-i <- 1 
		
	g <- Vectorize(fc.alpha)
	curve(g,0, 1, lty=1,lwd=1, main=""
	,xlab=expression(alpha),ylab="",ylim=c(0,1),col="red",cex=2,cex.lab=2)

	j <- i + 2
	g <- Vectorize(fc.alpha)
	curve(g,0, 1, lty=2,lwd=1,xlab=expression(alpha),ylab="",add=TRUE,col="red",cex=2,cex.lab=2)


	j <- i + 4
	g <- Vectorize(fc.alpha)
	curve(g,0, 1, lty=3,lwd=1,xlab=expression(alpha),ylab="",add=TRUE,col="red",cex=2,cex.lab=2)
	abline(0,1,col="blue")


	legend("topleft",legend=c("n=3","n=10","n=30"),lty=c(1,2,3),lwd=c(2,2),col="red")	
	mtext("Coverage",side=2,line=2,cex=1.5,cex.lab=1.5)


	j  <-i <- 2 
	
	
	g <- Vectorize(fc.alpha)
	curve(g,0, 1, lty=1,lwd=1, main=""
	,xlab=expression(alpha),ylab="",ylim=c(0,1),col="black",cex=2,cex.lab=2)

	j <- i + 2
	g <- Vectorize(fc.alpha)
	curve(g,0, 1, lty=2,lwd=1,xlab=expression(alpha),ylab="",add=TRUE,col="black",cex=2,cex.lab=2)
	legend("topleft",legend=c("n=3","n=10","n=30"),lty=c(1,2,3),lwd=c(2,2),col="red")


	j <- i + 4
	g <- Vectorize(fc.alpha)
	curve(g,0, 1, lty=3,lwd=1,xlab=expression(alpha),ylab="",add=TRUE,col="black",cex=2,cex.lab=2)
	abline(0,1,col="blue")

	mtext("Coverage",side=2,line=2,cex=1.5,cex.lab=1.5)		
	
	legend("topleft",legend=c("n=3","n=10","n=30"),lty=c(1,2,3),lwd=c(2,2),col="black")


	j <-i <- 1
	j <- i + 4
	plot(density(res1[[2]][,j]),col="red",lty=3,lwd=1,main="", xlab=expression(theta),ylab="", cex.lab=1.5,cex=1.5)

	j <- i + 2 
	lines(density(res1[[2]][,j]),col="red",lty=2,lwd=1, xlab=expression(theta),ylab="", cex.lab=1.5,cex=1.5)
	j <- i 
	lines(density(res1[[2]][,j]),col="red",lty=1,lwd=1, xlab=expression(theta),ylab="", cex.lab=1.5,cex=1.5)

	legend("topleft",legend=c("n=3","n=10","n=30"),lty=c(1,2,3),lwd=c(2,2),col="red")

	mtext(expression(paste(pi,"(", theta, "|", "x",")")),side=2,line=2,las=1,cex=1.5,cex.lab=1.5)

	j <-i <- 2
	j <-i + 4
	plot(density(res1[[2]][,j]),col="black",lty=3,lwd=1, main="",xlab=expression(theta),ylab="", cex.lab=1.5,cex=1.5)

	j <- i + 2 
	lines(density(res1[[2]][,j]),col="black",lty=2,lwd=1, xlab=expression(theta),ylab="", cex.lab=1.5,cex=1.5)
	j <- i
	lines(density(res1[[2]][,j]),col="black",lty=1,lwd=1, xlab=expression(theta),ylab="", cex.lab=1.5,cex=1.5)

	legend("topleft",legend=c("n=3","n=10","n=30"),lty=c(1,2,3),lwd=c(2,2),col="black")

	mtext(expression(paste(pi,"(", theta, "|", "x",")")),side=2,line=2,las=1,cex=1.5,cex.lab=1.5)
	





###################


fc.alpha <- function(alpha){
		
	o <- length( which( res1[[1]][,j] <= alpha ) )/m2

o
}


par(oma=c(0,3,0,0),mfrow=c(2,2))


	j  <-i <- i0 
		
	g <- Vectorize(fc.alpha)
	curve(g,0, 1, lty=1,lwd=1, main=""
	,xlab=expression(alpha),ylab="",ylim=c(0,1),col="red",cex=2,cex.lab=2)

	j <- i + 6
	g <- Vectorize(fc.alpha)
	curve(g,0, 1, lty=2,lwd=1,xlab=expression(alpha),ylab="",add=TRUE,col="red",cex=2,cex.lab=2)


	j <- i + 12
	g <- Vectorize(fc.alpha)
	curve(g,0, 1, lty=3,lwd=1,xlab=expression(alpha),ylab="",add=TRUE,col="red",cex=2,cex.lab=2)
	abline(0,1,col="blue")


	legend("topleft",legend=c("n=3","n=10","n=30"),lty=c(1,2,3),lwd=c(2,2),col="red")	
	mtext("Coverage",side=2,line=2,cex=1.5,cex.lab=1.5)


	j  <-i <- i0 +1 
	
	
	g <- Vectorize(fc.alpha)
	curve(g,0, 1, lty=1,lwd=1, main=""
	,xlab=expression(alpha),ylab="",ylim=c(0,1),col="black",cex=2,cex.lab=2)

	j <- i + 6
	g <- Vectorize(fc.alpha)
	curve(g,0, 1, lty=2,lwd=1,xlab=expression(alpha),ylab="",add=TRUE,col="black",cex=2,cex.lab=2)
	legend("topleft",legend=c("n=3","n=10","n=30"),lty=c(1,2,3),lwd=c(2,2),col="red")


	j <- i + 12
	g <- Vectorize(fc.alpha)
	curve(g,0, 1, lty=3,lwd=1,xlab=expression(alpha),ylab="",add=TRUE,col="black",cex=2,cex.lab=2)
	abline(0,1,col="blue")

	mtext("Coverage",side=2,line=2,cex=1.5,cex.lab=1.5)		
	
	legend("topleft",legend=c("n=3","n=10","n=30"),lty=c(1,2,3),lwd=c(2,2),col="black")


	j <-i <- i0
	j <- i + 12
	plot(density(res1[[2]][,j]),col="red",lty=3,lwd=1,main="", xlab=expression(theta),ylab="", cex.lab=1.5,cex=1.5)

	j <- i + 6
	lines(density(res1[[2]][,j]),col="red",lty=2,lwd=1, xlab=expression(theta),ylab="", cex.lab=1.5,cex=1.5)
	j <- i 
	lines(density(res1[[2]][,j]),col="red",lty=1,lwd=1, xlab=expression(theta),ylab="", cex.lab=1.5,cex=1.5)

	legend("topleft",legend=c("n=3","n=10","n=30"),lty=c(1,2,3),lwd=c(2,2),col="red")

	mtext(expression(paste(pi,"(", theta, "|", "x",")")),side=2,line=2,las=1,cex=1.5,cex.lab=1.5)

	j <-i <- i0 + 1
	j <-i + 12
	plot(density(res1[[2]][,j]),col="black",lty=3,lwd=1, main="",xlab=expression(theta),ylab="", cex.lab=1.5,cex=1.5)

	j <- i + 6 
	lines(density(res1[[2]][,j]),col="black",lty=2,lwd=1, xlab=expression(theta),ylab="", cex.lab=1.5,cex=1.5)
	j <- i
	lines(density(res1[[2]][,j]),col="black",lty=1,lwd=1, xlab=expression(theta),ylab="", cex.lab=1.5,cex=1.5)

	legend("topleft",legend=c("n=3","n=10","n=30"),lty=c(1,2,3),lwd=c(2,2),col="black")

	mtext(expression(paste(pi,"(", theta, "|", "x",")")),side=2,line=2,las=1,cex=1.5,cex.lab=1.5)
	


fc.m <- matrix(NA,3,2)
j <- 7
fc.m[1,1] <- fc.alpha(0.95)
j<- j + 1
fc.m[1,2] <- fc.alpha(0.95)

j <- 15

fc.m[2,1] <- fc.alpha(0.95)
j<- j + 1
fc.m[2,2] <- fc.alpha(0.95)

j <- 23
fc.m[3,1] <- fc.alpha(0.95)

j<- j + 1
fc.m[3,2] <- fc.alpha(0.95)



par(mfrow=c(2,2))

for(i in c(2,4,6)){


	j  = i
	
	
	g <- Vectorize(fc.alpha)
	curve(g,0, 1, lty=2,lwd=2, main=
	bquote("True Value"~theta== 5),xlab=expression(alpha),ylab="Coverage",ylim=c(0,1))

	j <- i + 6
	g <- Vectorize(fc.alpha)
	curve(g,0, 1, lty=3,lwd=2,xlab=expression(alpha),ylab="Frequency Coverage",add=TRUE)


	j <- i + 12
	g <- Vectorize(fc.alpha)
	curve(g,0, 1, lty=4,lwd=2,xlab=expression(alpha),ylab="Frequency Coverage",add=TRUE)
	abline(0,1,col="blue")

		
	
}

	legend("topleft",legend=c("3","10","30"),lty=c(2,3,4),lwd=c(2,2),pt.cex=1,cex=0.6)

