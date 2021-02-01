# Modified on Jan 14, 2019
# Located in ~/Desktop/Assortative_Mating/R/version_5/ folder

## Input
# A: matrix to be scaled
# r.sum: target row sum vector
# c.sum: target colum sum vector
# The algorithm can work for rectangular matrix but for now 
# assume square matrix input.

## Output
# result: scaled matrix

library('Rmosek')


convex.opt <- function(A, r.sum, c.sum) {
    # initialization 
    n <- dim(A)[1]
    A.vec <- c(A) # Flatten the matrix to 1d-array. It's done by "column"
    # A.vec <- A.vec/sum(A.vec)
    
    # generate linear row and column sum constraints matrix
    # print('entering constraint matrix generation routine')
    row.const <- Matrix(0, nc=length(A), nr=nrow(A)) # row sum mat
    row.const[cBind(as.vector(row(A)),1:length(A))] <- 1
    col.const <- Matrix(0, nc=length(A), nr=nrow(A)) # col sum mat 
    col.const[cBind(as.vector(col(A)),1:length(A))] <- 1
    C.mat <- rBind(row.const, col.const)
    C.vec <- c(r.sum, c.sum)
    # print('contraint matrix has been created')

    # set parameters for mosek
    sco1 <- list(sense = "min")
    sco1$c <- -log(A.vec)
    sco1$A <- C.mat
    sco1$bc <- rbind(blc = C.vec,
                     buc = C.vec)
    sco1$bx <- rbind(blx = rep(0,n^2),
                     bux = rep(1, n^2))
    NUMOPRO <- n^2; 
    opro <- matrix(list(), nrow=5, ncol= NUMOPRO)
    rownames(opro) <- c("type","j","f","g","h")
    for (i in 1:NUMOPRO) {
        opro[,i] <- list("ENT", i, 1.0, NA, NA)
    }
    # OPRC not needed since there is no non-linear constraints
    sco1$scopt <- list(opro=opro)
    sco1$iparam <- list(presolve_use='MSK_PRESOLVE_MODE_OFF',
                        intpnt_basis='MSK_BI_NEVER',
                        intpnt_max_iterations=500)
    
    # run mosek
    print('entering convex optimization')
    r <- mosek(sco1, list(verbose=1)) #Errors=1, MOSEK=2, Warnings=3 and Info=4, default=10
    # matrix(r$sol$itr$xx, nrow=n)
    print('convex optimization completed')
    
    if (is.null(r$sol$itr$xx)) {
        result <- NULL
    } else {
        result <- matrix(r$sol$itr$xx, nrow=n)
    }
    return(result)
}


matrix.scale <- function(A, r.sum, c.sum, tol=1e-8, max.iter=5000) {
    # The problem is in the form of 
    # minimize: xlog(x) - xlog(a)
    # subject to: row.sum(x_ij) = r.sum_i, col.sum(x_ij) = c.sum_j
    # bound: 0 <= x_ij <= 1
    
    # print('entering matrix scaling routine')
    # input check
    stopifnot(sum(r.sum) == sum(c.sum))
    stopifnot(dim(A)[1] == dim(A)[2])
    
    result <- convex.opt(A, r.sum, c.sum)
    
    return(result)
}
