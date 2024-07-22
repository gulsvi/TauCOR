##################################################################
##    Conditinal GBA analysis: TauCOR calculates Tau and 
##       corrects LD matrix between SNPs in gene
##################################################################

## n   : sample size
## Z.g : vector of z-statistics for SNPs of gene
## Z.r : vector of z-statistics for SNPs of region around gene
## U.g : LD matrix between SNPs of gene
## U.r : LD matrix between SNPs of region around gene
## U.rg: LD matrix between SNPs of region around gene and SNPs of gene
## h2  : total trait heritability

 
####################################################
###  approx.estimation of tau by maxLh           ###
####################################################
TauCOR <- function(Z.r, U.r, U.rg, U.g,  n, eps = 0){

	LH <- function (Tau){ 
		Lam   <- eigvalU2 * Tau + eigvalU             
		logLH <-  sum(log(Lam) + ZrZr/Lam)            ## to minimize -log Likelihood
		return(logLH)
	}

	REF 	<- eigen(U.r, sym = TRUE)            
	eigvalU <- REF$val
	eigvecU <- REF$vec 
	eigvecU <- eigvecU[, eigvalU > eps] 	
	eigvalU <- eigvalU[  eigvalU > eps]
	rank.Ur <- length(eigvalU) 
	ndm     <- n/rank.Ur
 	Zr      <- t( t(Z.r) %*% eigvecU )

	eigvalU2  <- eigvalU^2  
	ZrZr      <- Zr^2
	B <- optimize(f = LH, interval = c(0, ndm),maximum = FALSE,tol=1e-22)
	Tau <- B$minimum 
	mat <-   crossprod(U.rg) * Tau  + U.g

	return(list(Tau=Tau,correctedU=mat))
}

	results  <- TauCOR(Z.r, U.r, U.rg, U.g,  n, eps = 0){
  
 
####################################################################
## SimSumStat: The function to simulate z-statistics and p-values 
## for the window (gene + region around gene)
## SimSumStat uses the Linear Regresion Model with Random SNP Effects   
####################################################################

SimSumStat <- function(rho, h2, Ncondi, U, M.r, M.g, n, mmm) {
	## rho    - the proportion of causal SNPs in gene among all causal SNPs in window. for example (0, 0.5 or 1)
	## h2     - the total trait heritability, for example (0.3, 0.5 or 0.7)
	## Ncondi - the total number of causal SNPs in window, for example (5 or 10)
	## U      - the LD matrix between SNPs in window (original or simulated LD matrix of size M x M)
	## M.r    - the total number of  SNPs in region
	## M.g    - the total number of  SNPs in gene

	M        <- M.g + M.r                     # the total number of SNPs in analysis
	Mp.g     <- round(rho *   Ncondi)         # the number of causal SNPs in gene
	Mp.r     <- Ncondi - Mp.g                 # the number of causal SNPs in region
	he2.reggen <- h2 * M/mmm                  # the heritability explained by region+gene	
	Tau      <- n * he2.reggen / Ncondi       # n x (heritability explained by one causal SNP)
	sqrt.Tau <- sqrt(Tau)     
	eigen.U  <- eigen(U, symmetric = TRUE)
	sqrt.eig.values <- sqrt(pmax(eigen.U$values,0))
	eig.vectors     <- eigen.U$vectors
	
	## simulating the causal statuses of SNPs
	
	if (Mp.r > 0) { c.r <- sample(c(rep(1, Mp.r), rep(0, M.r - Mp.r))) } else { c.r = rep(0,M.r) }
	if (Mp.g > 0) { c.g <- sample(c(rep(1, Mp.g), rep(0, M.g - Mp.g))) } else { c.g = rep(0,M.g) }
	cc  <- c(c.r, c.g)
    
	## simulating the joint z-statistics
	
	jz <- rep(0, M)             #  zeroing the vector of joint z-statistics
	X  <- as.vector(rnorm(Mp, mean=0, sd=1))
	#X  <- as.vector(X, center = TRUE, scale = TRUE))
	jz[which(cc == 1)] <- X 
	GWASz1 <- as.vector(U %*% jz)

	## simulating the GWAS z-statistics
	
    X <- as.vector(rnorm(M, mean = 0, sd = 1))
    #X <- as.vector(scale(X, center = TRUE, scale = TRUE))
	GWASz2 <- as.vector(eig.vectors %*% (sqrt.eig.values * X))
	GWASz <- GWASz1 * sqrt.Tau  + GWASz2 
	
	names(GWASz) <- colnames(U)
    GWASp <- as.vector(pchisq(GWASz^2, 1, lower.tail = FALSE))                    ## convert Z to Pval
	names(GWASp) <- colnames(U)
	return(list(Z = GWASz, P = GWASp))
}

