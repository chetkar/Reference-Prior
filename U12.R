##########################
# ARMS sampling
################################

library("HI")

rf.cond.lik <- function(th1,th2,x){
	
	k  <- length(x)
	t1 <- min(x)
	t2 <- max(x)


	ans <- -(k + 1)*log(th1)

ans
}

rf.cond.unif <- function(th1, th2, x){

	k  <- length(x)
	t1 <- min(x)
	t2 <- max(x)

	ans <- -(k)*log(th1*(th2- 1)) 
ans
}

rf.marg.unif.lik <- function(th2, x){

	k  <- length(x)
	t1 <- min(x)
	t2 <- max(x)

	ans <- -(k )*log(th2 - 1) + log( (th2/t2)^{k-1} - (1/t1)^{k-1})
	ans  <- exp(ans)
	
ans
}

cset <- function(x,y){
	
	t1 <- min(y)
	t2 <- max(y)
	(x > t2/t1)*(x < 50)

}

rf.marg <- function(x,y){


	k  <- length(y)
	t1 <- min(y)
	t2 <- max(y)

	ans <- -(k + 1)*log(x - 1) - log(x) + log( (x/t2)^k - (1/t1)^k)
	ans  <- ans -(x+1)/x
	
	
ans
}

rf.marg.unif <- function(x,y){


	k  <- length(y)
	t1 <- min(y)
	t2 <- max(y)

	ans <- -(k)*log(x - 1) - log(x) + log( (x/t2)^k - (1/t1)^k)
	ans  <- ans 
	
	
ans
}


rf.marg.lik <- function(th2,x){


	k  <- length(x)
	t1 <- min(x)
	t2 <- max(x)

	ans <- -(k + 1)*log(th2 - 1) - log(th2) 
	ans  <- exp(ans)*exp(-(th2+1)/th2)*( (th2/t2)^k - (1/t1)^k)
	
	
ans
}


rf.post <- function(th2,x){

	t1 <- min(x)
	t2 <- max(x)

	low <- t2/t1
	high <- Inf
	ans <- rf.marg.lik(th2,x)/integrate(rf.marg.lik,low,high,x)$value

ans
}


rf.post2 <- function(th2,x){
	t1 <- min(x)
	t2 <- max(x)

	low  <- t2/t1
	high <- 20
	ans  <- rf.marg.unif.lik(th2,x)/integrate(rf.marg.unif.lik,low,high,x)$value

ans
}

#########################
# Sampling from Proposal#
#########################
m1 <- 2000
m2 <- 2000

n.vec   <- c(3,10,30)
th1.vec <- rep(c(2,5,8), each =3)
th2.vec <- rep(c(2,5,8), 3)

th2.est   <- th1.est   <- rep(NA, m1 + 1)
th2.est.m <- th1.est.m <- rep(NA,m1/2)

th22.est   <- th12.est   <- rep(NA, m1 + 1)
th22.est.m <- th12.est.m <- rep(NA,m1/2)


rho11 <- rho12 <- rho21 <- rho22  <- rep(NA, m2)
mse <- matrix(NA,27,4)

m.4 <- m.3 <- m.1 <- m.2 <- rep(NA, m2)
x.list <- th2.list <- th1.list <- Out.result <- list()


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
    x <- runif(n, th1.vec[l], th1.vec[l]*th2.vec[l])
	
    t1 <- min(x)
    t2 <- max(x)
		
    y<- x
    th2.est  <- arms(1.2*t2/t1,rf.marg,     cset, m1,y)			
    th22.est <- arms(1.2*t2/t1,rf.marg.unif,cset, m1,y)


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

   rho11[i]    <- (2*length( which( th1.est.m <= theta_2) ) )/m1
   rho12[i]    <- (2*length( which( th12.est.m <= theta_2) ) )/m1

   rho21[i]    <- integrate(rf.post,t2/t1,theta_2,x)$value
   
   rho22[i]    <- integrate(rf.post2,t2/t1,theta_2,x)$value

 	
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


