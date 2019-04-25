##Model Unif(th1, th1*th2)
## {th1, th2 }


rf.cond <- function(th1, th2, x){
	
	k  <- length(x)
	t1 <- min(x)
	t2 <- max(x)

	ans <- - (k + 1)*log(th2 - 1) 

ans
}

rf.cond.unif <- function(th1, th2, x){
	
	k  <- length(x)
	t1 <- min(x)
	t2 <- max(x)

	ans <-  - (k )*log(th2 - 1)

ans
}

rf.marg.unif <- function(th1, x){

	k  <- length(x)
	t1 <- min(x)
	t2 <- max(x)
	
	
	ans <-  -(k-1)*log(t2 - th1)-log(th1)
	ans <- exp(ans)
	

ans
}



rf.marg.lik <- function(th1, x){

	k  <- length(x)
	t1 <- min(x)
	t2 <- max(x)
	
	ans <-  -k*log(t2 - th1)
	ans <- exp(ans)

ans
}

rf.post <- function(th1,x){

	t1 <- min(x)
	t2 <- max(x)

	low <- 0
	high <- t1
	ans <- rf.marg.lik(th1,x)/integrate(rf.marg.lik,low,high,x)$value

ans
}


rf.post2 <- function(th1,x){
	t1 <- min(x)
	t2 <- max(x)

	low <- 0
	high <- t1
	ans <- rf.marg.unif(th1,x)/integrate(rf.marg.unif,low,high,x)$value

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

#n.vec   <- c(10,30)
#th1.vec <- c(2,8)
#th2.vec <- c(8,2)

th2.est   <- th1.est   <- rep(NA, m1 + 1)
th2.est.m <- th1.est.m <- rep(NA,m1/2)

th22.est   <- th12.est   <- rep(NA, m1 + 1)
th22.est.m <- th12.est.m <- rep(NA,m1/2)


rho11 <- rho12 <- rho21 <- rho22  <- rep(NA, m2)
mse <- matrix(NA,27,4)

m.4 <- m.3 <- m.1 <- m.2 <- rep(NA, m2)
th2.list <- th1.list <- Out.result <- list()
x.list <-  list()


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
    x.list[[count]] <- x	
	
    th12 <- th1 <- min(x)
    th22 <- th2 <- max(x)/min(x)

    th1.low  <- 0
    th1.high <- min(x)


   for(j in 1:m1){

#########################
# rf prior
##########################

		j.p  <- j + 1
		t1 <- min(x)
		t2 <- max(x)

		th1_n    <- runif(1,0, t1)
		th12_n   <- runif(1,0, t1)


		rat1     <- rf.marg.lik(th1_n,x)*(t2 - t1)^{n}
		rat12    <- rf.marg.unif(th12_n,x)*(t2 - t1)^{n-1}

		u <- runif(1,0, 1)

		if( u < rat1 ){
			th1.est[j.p ] <- th1_n
			th1         <- th1_n
		}else{
			th1.est[j.p]  <- th1

		}

		if( u < rat12 ){
			th12.est[j.p ] <- th12_n
			th12         <- th12_n
		}else{
			th12.est[j.p]  <- th12

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


	rho11[i]    <- integrate(rf.post,0,theta_1,x)$value
	#rho12[i]    <- integrate(rf.post2,0,theta_1,x)$value

	rho21[i]    <- (2*length( which( th1.est.m <= theta_1) ) )/m1
	rho22[i]    <- (2*length( which( th12.est.m <= theta_1) ) )/m1


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

############ Plots#############33
###############################
par(mfrow=c(2,2))


i <- 1
j <-1
g <- Vectorize(fc.alpha)
curve(g,0, 1, lty=1,lwd=2, col ="red",xlab=expression(alpha),ylab="Coverage",main="Case 1 =(1,2,(n=3))")


i <- 1
j <-2
g <- Vectorize(fc.alpha)
curve(g,0, 1, lty=1,lwd=2 ,col ="black", add =TRUE)

abline(0,1, col="blue")

legend("topleft",legend=c("Prior1","Ref"),col=c("red","black"),lty=c(2,2),lwd=c(2,2))

i <- 2
j <-1
g <- Vectorize(fc.alpha)
curve(g,0, 1, lty=1,lwd=2 ,col ="red",xlab=expression(alpha),ylab="Coverage",main="Case 2 =(5,3,(n=3))")


i <- 2
j <-2
g <- Vectorize(fc.alpha)
curve(g,0, 1, lty=1,lwd=2 ,col ="black", add =TRUE)

abline(0,1, col="blue")

legend("topleft",legend=c("Prior1","Ref"),col=c("red","black"),lty=c(2,2),lwd=c(2,2))

i <- 3
j <-1
g <- Vectorize(fc.alpha)
curve(g,0, 1, lty=1,lwd=2 ,col ="red",xlab=expression(alpha),ylab="Coverage",main="Case 3 =(1,2,(n=10))")


i <- 3
j <-2
g <- Vectorize(fc.alpha)
curve(g,0, 1, lty=1,lwd=2 ,col ="black", add =TRUE)

abline(0,1, col="blue")

legend("topleft",legend=c("Prior1","Ref"),col=c("red","black"),lty=c(2,2),lwd=c(2,2))

i <- 4
j <-1
g <- Vectorize(fc.alpha)
curve(g,0, 1, lty=1,lwd=2 ,col ="red", xlab=expression(alpha),ylab="Coverage",main="Case 4 =(5,3,(n=10))")


i <- 4
j <-2
g <- Vectorize(fc.alpha)
curve(g,0, 1, lty=1,lwd=2 ,col ="black", add =TRUE)

abline(0,1, col="blue")

legend("topleft",legend=c("Prior1","Ref"),col=c("red","black"),lty=c(2,2),lwd=c(2,2))


