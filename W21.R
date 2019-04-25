##Model Truncated Weibull
##{scale, theta}
##{theta, scale}


rf.cond <- function(th1, th2, x){
	
     k <- length(x)	
     y   <- x	
     ans <- - (k )*log(th1) - sum((y - th2)/th1)

ans
}

rf.cond.unif <- function(th1, th2, x){
	
     k <- length(x)	
     y   <- x	
     ans <- - (k)*log(th1) - sum((y - th2)/th1)

ans
}


rf.marg.unif.lik <- function(th1, th2, x){

     k <- length(x)	
     y   <- x	
     t1  <- min(x)	

     ans <- (k-1)*log(th1)
    ans <- exp(ans)*(exp(-th1*sum(y-t1))- exp(-th1*sum(y)) )

ans
}

rf.marg.lik <- function(th1, th2, x){

     k <- length(x)	
     y   <- x	
     t1  <- min(x)	
    
     ans <- (k-2)*log(th1)
     ans <- exp(ans)*(exp(-th1*sum(y-t1))- exp(-th1*sum(y)) )	


ans
}




#########################
# Sampling from Proposal#
#########################
m1 <- 2000
m2 <- 2000

n.vec   <- c(3,10,30)
#rate
#truncation
th1.vec <- rep(c(1,2,4), each =3)
th2.vec <- rep(c(1,2,4), 3)

#phi <- 2

#n.vec   <- c(3,10)
#th1.vec <- c(2,8)
#th2.vec <- c(8,2)


th2.est   <- th1.est   <- rep(NA, m1 + 1)
th2.est.m <- th1.est.m <- rep(NA,m1/2)

th22.est   <- th12.est   <- rep(NA, m1 + 1)
th22.est.m <- th12.est.m <- rep(NA,m1/2)


rho11 <- rho12 <- rho21 <- rho22  <- rep(NA, m2)
mse <- matrix(NA,27,4)

m.4 <- m.3 <- m.1 <- m.2 <- rep(NA, m2)
x.list<- th2.list <- th1.list <- Out.result <- list()


count <- 0 
for( k in 1:length(n.vec) ){
for( l in 1:length(th1.vec) ){

   count <- count + 1		
   for(i in 1:m2){

   #################################
   # Simulate X
   #################################

    n <- n.vec[k]	
    theta_1 <- th1.vec[l]	
    theta_2 <- th2.vec[l]
    	
    x <- rtrunc(n,"exp", rate= th1.vec[l], th2.vec[l], Inf)
    x.list[[count]] <- x	
	
    th22 <- th2 <- runif(1, 1, min(x))	
    th12 <- th1 <- mean(x - th2)

    th1.low  <- 0
    th1.high <- min(x)


   for(j in 1:m1){

#########################
# rf prior
##########################

		j.p  <- j + 1
		
		t1 <- min(x)
		
		th1_n    <- rgamma(1,shape= length(x)-1,rate=sum(x - t1))
		th12_n   <- rgamma(1,shape= length(x),rate=sum(x - t1))

			
		M        <- factorial(length(x)-2)/(sum(x-t1))^{length(x)-1}
		rat1     <- rf.marg.lik(th1_n, th2, x)/(M*dgamma(th1_n, shape=length(x)-1, rate=sum(x-t1)))


		M        <- factorial(length(x)-1)/(sum(x-t1))^{length(x)}
		rat12     <- rf.marg.unif.lik(th12_n, th2, x)/(M*dgamma(th12_n, shape=length(x), rate=sum(x-t1)))

		u <- runif(1,0, 1)


		if( u < rat1 ){
			th1.est[j.p ] <- th1_n
			th1         <- th1_n
		}else{
			th1.est[j.p]  <- th1
		}

		if( u < rat12 ){
			th12.est[j.p ] <- th12_n
			th12           <- th12_n
		}else{
			th1.est[j.p]   <- th12
		}


		th2_n    <- rtnorm(1,t1,1,-Inf, t1)
		th22_n   <- rtnorm(1,t1,1,-Inf, t1)

		rat2    <- exp(rf.cond(th1, th2_n, x)  - rf.cond(th1, th2, x) )
		rat22   <- exp(rf.cond.unif(th12, th22_n, x)  - rf.cond.unif(th12, th22, x) )
		

		if( u < rat2 ){
			th2.est[j.p ] <- th2_n
			th2_0         <- th2_n
			th2           <- th2_0

		}else{
			th2.est[j.p]  <- th2.est[j]
		}		  		 
   

		if( u < rat22 ){
			th22.est[j.p ] <- th22_n
			th22            <- th22_n

		}else{
			th2.est[j.p]  <- th22
		}		

}


	ll        <- m1/2 + 1
	th1.est.m <- th1.est[ll:m1]
	th2.est.m <- th2.est[ll:m1]


	ll         <- m1/2 + 1
	th12.est.m <- th12.est[ll:m1]
	th22.est.m <- th22.est[ll:m1]


	mse.1     <- mean( abs(th1.est.m - theta_1))
	mse.2     <- mean( abs(th2.est.m - theta_2))

	mse.12     <- mean( abs(th12.est.m - theta_1))
	mse.22     <- mean( abs(th22.est.m - theta_2))

	m.1[i]    <- mse.1
	m.2[i]    <- mse.12
	m.3[i]    <- mse.2
	m.4[i]    <- mse.22

	rho11[i]    <- (2*length( which( th1.est.m <= theta_1) ) )/m1
	rho12[i]    <- (2*length( which( th12.est.m <= theta_1) ) )/m1

	rho21[i]    <- (2*length( which( th2.est.m <= theta_2) ) )/m1
	rho22[i]    <- (2*length( which( th22.est.m <= theta_2) ) )/m1

	if(i %% 100 == 0 ){
		print( paste(Sys.time(), "Iteration number:",i,k,l) )			
	}

}
	th1.list[[count]] <- cbind(th1.est.m, th12.est.m)
	th2.list[[count]] <- cbind(th2.est.m, th22.est.m)
	Out.result[[count]] <- cbind(rho11, rho12, rho21, rho22)
	mse[count,] <- c(mean(m.1), mean(m.2), mean(m.3),mean(m.4))

}
}


fc.alpha <- function(alpha){
		
	o <- length( which( Out.result[[i]][,j] <= alpha ) )/m2

o
}


