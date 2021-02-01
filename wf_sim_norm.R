# Modified on Jan 14, 2019
# Located in ~/Desktop/Assortative_Mating/R/version_5/ folder

# This is a code for simulating an assortative mating in DIPLOIDS 
# "norm" stands for normalized version. It still supports the parallel 
# version as well.

# ======================== Input Parameters ========================
# m.mode: Mating mode. One of 
#         'random' (random mating), 
#         'gen' (assortative mating by genotype), 
#         'phen' (assortative mating by phenotype)
# N.pop: number of finite population size in the simulation 
# N.sim: number of simulations to perform
# G: The number of generations
# s10: The initial fraction from population 1 
#      s20 is automatically evaluated as s20 = 1 - s10
# s1g: The constant contributions from population 1 at g > 1
# s2g: The constant contributions from population 2 at g > 1
#      hg (contributions from the hybrid population) is automatically 
#      determined by hg = 1 - s1g - s2g
# am.const: constant (c) in assortative mating function
#           currently exp(-c * abs(t1 - t2)), this can be changed at 
#           init.R file
# k: number of loci contributing to the trait
# x.prob: probability for x_i = 1 (binomial random) for each k-trait loci.
#         default=1.0
# test.mode: if true, the allele frequencies are fixed as ('1','0')=(1,0) 
#            across all loci. default=T
# pld: ploidy (extension from the haploid model)
# rand.seed <- random seed for randomly generating L and X values, defult=1
# n.core <- number of clusters for parallel computing, default=4
# plt <- if we want to plot P(H|X) on the fly, default=T

# ======================== Output Parameters ========================
# HA: admixture fraction values from the soure population A 
# HA.P: P(H|X), probability for a random individual to have a given H value
# E.HA : E[H|X], expectation of admixture fraction 
# Var.HA : Var[H|X], variance of admixture fraction 
# E.Tr: E[T|X], expectation of trait
# Var.Tr: Var[T|X], variance of trait
# Tr.D.HA: P(T|H,X), trait distribution per given admixture fraction. 
#       a matrix with a size of a different number of admixture fraction in 
#       the generation x number of different traits (k+1))
# Tr.D: P(T|X), trait distribution in the population (not conditioning on 
#       the admixture fraction)
# P.TH: P(T,H|X), joint distribution of trait and admixture fraction
# Cov.HT: covaraince of trait and admixture fraction
# Corr.HT: correlation of trait and admixture fraction
# Corr.mate.Tr: correlation of traits between mates
# Corr.mate.HA: correlation of admixture fraction between mates
# A.Mat: P(L_i=1|H=h, X=x), alpha.mat for each generation, A.Mat[[g]] has
#       a dimension of length(h.vec) * k
# AL.F.1: E(L_i=1|X), allele frequency of allele '1'
# AL.F.0: E(L_i=0|X), allele frequency of allele '0'
#
# ===================================================================
# An individual should be tagged with:
# adimixture fraction (col=1), trait (col=2), 
# genotype 1, (col=3:(k+2)), genotype 2 (col=(k+3):(2k+2))
# ===================================================================

library(plyr)
library(gRbase)
library(parallel)
library(microbenchmark)
library(foreach)
library(doParallel)
library(doSNOW)
library(data.table)
library(tcltk)


wf.single.norm <- function(m.mode, N.pop, G, k, s10, s1g=0, s2g=0,
                      am.const=0.5, pld=2, loci.par, allele.freq) {
    
    library(plyr)
    ## =========== Initialization ===========
    G <- 1:G
    s20 <- 1 - s10
    hg <- 1 - s1g - s2g
    
    # Initialize variables to save
    # result <- list()
    # The length of the vectors are G+1 as we also record relevant statistics for 
    # the founding populations (i.e. G=0)
    Var.HA <- rep(NA, times=length(G)+1) # variance of admixture fraction
    E.HA <- rep(NA, times=length(G)+1) # expectation of admixture fraction
    Var.Tr <- rep(NA, times=length(G)+1) # variance of trait
    E.Tr <- rep(NA, times=length(G)+1) # expectation of trait
    Cov.HT <- rep(NA, times=length(G)+1) # covariance between admixture fraction and trait
    Cor.HT <- rep(NA, times=length(G)+1) # correlation between admixture fraction and trait
    Cor.mate.Tr <- rep(NA, times=length(G)+1) # correlation of traits between mates
    Cor.mate.HA <- rep(NA, times=length(G)+1) # correlation of admixture fraction between mates
    Avg.Homo.AL <- rep(NA, times=length(G)+1) # average homozygosity based on allele freq across all loci
    Avg.Hetero.AL <- rep(NA, times=length(G)+1) # average heterozygosity based on allele freq across all loci
    Avg.Homo.Gen <- rep(NA, times=length(G)+1) # average homozygote frequency based on genotype across all loci
    Avg.Hetero.Gen <- rep(NA, times=length(G)+1) # average heterozygote frequency  based on genotype across all loci
    
    # par.pool is a matrix format in the following order:
    # col.1: admixture fraction
    # col.2: NA (col reserved for trait value)
    # col.3~(k+2): genotype 1 of a given individual
    # col.(k+3)~(2k+2): genotype 2 of a given individual
    par.pool.f <- matrix(NA, nrow=N.pop, ncol=pld*k+2) #female potental parent
    par.pool.m <- matrix(NA, nrow=N.pop, ncol=pld*k+2) #male potental parent
    
    
    colnames(par.pool.f) <- c('am', 'tr', c(sapply(1:pld, 
                                                   function(x) paste(x, c(1:k), sep='.'))))
    colnames(par.pool.m) <- c('am', 'tr', c(sapply(1:pld, 
                                                   function(x) paste(x, c(1:k), sep='.'))))
    
    # founding population i.e. g=0
    n.1 <- round(N.pop * s10) #number of initial population from pop 1
    n.2 <- N.pop - n.1 #number of initial population from pop 2
    par.pool.f[1:n.1,] <- gen.pop.sample('A', n.1, k, allele.freq, ploidy=pld) # female parent
    par.pool.f[(n.1+1):N.pop,] <- gen.pop.sample('B', n.2, k, allele.freq, ploidy=pld) # female parent
    par.pool.m[1:n.1,] <- gen.pop.sample('A', n.1, k, allele.freq, ploidy=pld) # male parent
    par.pool.m[(n.1+1):N.pop,] <- gen.pop.sample('B', n.2, k, allele.freq, ploidy=pld) # male parent
    
    # compute trait from genotype (need to implement X portion later)
    par.pool.f[,2] <- eval.trait(par.pool.f[ ,c(-1,-2), drop=F], loci.par) # female
    par.pool.m[,2] <- eval.trait(par.pool.m[ ,c(-1,-2), drop=F], loci.par) # male
    
    par.pool.tot <- rbind(par.pool.f, par.pool.m) # total parental pool combining both genders
    
    if (any(is.na(par.pool.f)) || any(is.na(par.pool.m))) {
        stop("something went wrong with initial population assignment")
    }
    
    # Compute HA and allele frequency
    # ha.dist <- count(par.pool$HA, vars=h.vec)
    
    Var.HA[1] <-  var(par.pool.tot[,1])
    E.HA[1] <- mean(par.pool.tot[,1]) 
    Var.Tr[1] <- var(par.pool.tot[,2]) 
    E.Tr[1] <- mean(par.pool.tot[,2]) 
    Cov.HT[1] <- cov(par.pool.tot[,1], par.pool.tot[,2])
    Cor.HT[1] <- cor(par.pool.tot[,1], par.pool.tot[,2])
    
    # the first index for the correlation between mates are meaningness as 
    # there was no mating but for the data structure of the results
    # we'll assign it as 0 for now. 
    Cor.mate.Tr[1] <- 0 
    Cor.mate.HA[1] <- 0 
    
    # Average homozygosity and heterozygosity in the population across all loci
    Avg.Homo.AL[1] <- eval.avg.homo(par.pool.tot[,c(-1,-2), drop=F], k)
    Avg.Hetero.AL[1] <- 1 - Avg.Homo.AL[1]
    Avg.Homo.Gen[1] <- eval.genotype.freq(par.pool.tot[,c(-1,-2), drop=F], k)
    Avg.Hetero.Gen[1] <- 1 - Avg.Homo.Gen[1]
    
    
    for (g in G) {
        print(paste("generation = ", g))
        child.f <- matrix(NA, nrow=N.pop, ncol=pld*k+2) #female children
        child.m <- matrix(NA, nrow=N.pop, ncol=pld*k+2) #male children
        
        # pick N.pop # of pairs for mating
        if (m.mode == 'random') {
            # print('random mating')
            #             mate.ind <- matrix(unlist(lapply(1:N.pop, function(x) 
            #             {sample.int(N.pop, 2, replace=F)})), ncol = 2, byrow=T)
            mate.ind <- sample.mates(N.pop=N.pop, N.sample=N.pop)
            mate.ind <- data.table(mate.ind)
            names(mate.ind) <- c('Mf', 'Mm')
            
        } else if (m.mode == 'gen') {
            # print('assortative mating by amdixture fraction')
            if (sd(par.pool.tot[,1]) == 0) {
                scale <- 1
            } else {
                scale <- sd(par.pool.tot[,1])
            }
            mate.ind <- gen.offspring.norm(N.pop, mate.data.f=par.pool.f[,1]/scale,
                                           mate.data.m=par.pool.m[,1]/scale,
                                           am.const=am.const)
            
        } else if (m.mode == 'phen') {
            # print('assortative mating by phenotype')
            # phen <- par.pool[,2]/k # scale by the number of QT loci
            if (sd(par.pool.tot[,2]) == 0) {
                scale <- 1
            } else {
                scale <- sd(par.pool.tot[,2])
            }
            mate.ind <- gen.offspring.norm(N.pop, mate.data.f=par.pool.f[,2]/scale,
                                           mate.data.m=par.pool.m[,2]/scale,
                                           am.const=am.const)
            
        } else {
            stop('error: unsupported mating type')
        }
        
        
        ## admixture fraction of the children
        child.f[,1] <- (par.pool.f[mate.ind$Mf, 1] + par.pool.m[mate.ind$Mm, 1]) / 2 
        child.m[,1] <- (par.pool.f[mate.ind$Mf, 1] + par.pool.m[mate.ind$Mm, 1]) / 2 
        
        if(!all.equal(child.f[,1], child.m[,1])) {
            stop("The sibilings must have the same admixture fraction")
        }
        
        
        ## Genotype of children. For now, it can only handle diploid.
        #  multiple ploidy implementation will follow later. Need to clarify
        #  the gamet selection process. 
        al.pool.f <- par.pool.f[mate.ind$Mf, 3:(pld*k+2)] #mom's genotype
        al.pool.m <- par.pool.m[mate.ind$Mm, 3:(pld*k+2)] #dad's genotype
        # al.pool <- rbind(c(al.pool.1), c(al.pool.2))
        # pick.ind <- rbinom(N.pop*k, 1, prob=0.5) + 1
        # al.pick <- al.pool[cbind(pick.ind, 1:(N.pop*k))]
        # child[,3:(k+2)] <- matrix(al.pick, nrow=N.pop)
        
        sample.gene <- function(gen.ind) {
            temp <- matrix(gen.ind, nrow=pld, ncol=k, byrow=T)
            return(apply(temp, 2, sample, size=1))
        }
        
        child.f[, 3:(k+2)] <- aaply(al.pool.f, 1, sample.gene)
        child.f[, (k+3):(2*k+2)] <- aaply(al.pool.m, 1, sample.gene)
        
        child.m[, 3:(k+2)] <- aaply(al.pool.f, 1, sample.gene)
        child.m[, (k+3):(2*k+2)] <- aaply(al.pool.m, 1, sample.gene)
        
        ## Compute trait from genotype (need to implement X portion later)
        child.f[,2] <- eval.trait(child.f[ ,c(-1,-2), drop=F], loci.par) 
        child.m[,2] <- eval.trait(child.m[ ,c(-1,-2), drop=F], loci.par)
        
        
        # ========== Compute statistics for the current generation ==========
        child.tot <- rbind(child.f, child.m)
        
        # Compute variance and expectation of HA and Tr
        Var.HA[g+1] <-  var(child.tot[,1])
        E.HA[g+1] <- mean(child.tot[,1]) # expectation of admixture fraction
        Var.Tr[g+1] <- var(child.tot[,2]) # variance of trait
        E.Tr[g+1] <- mean(child.tot[,2]) # expectation of trait
        
        # Compute covariance and correlation between HA and Tr
        Cov.HT[g+1] <- cov(child.tot[,1], child.tot[,2])
        Cor.HT[g+1] <- cor(child.tot[,1], child.tot[,2])
        
        # Compute correlation between mates based on HA and Tr
        Cor.mate.HA[g+1] <- cor(par.pool.f[mate.ind$Mf, 1], par.pool.m[mate.ind$Mm, 1])
        Cor.mate.Tr[g+1] <- cor(par.pool.f[mate.ind$Mf, 2], par.pool.m[mate.ind$Mm, 2])
        
        # Compute average homozygosity and heterozygosity in the population across all loci
        Avg.Homo.AL[g+1] <- eval.avg.homo(child.tot[,c(-1,-2)], k)
        Avg.Hetero.AL[g+1] <- 1 - Avg.Homo.AL[g+1]
        
        Avg.Homo.Gen[g+1] <- eval.genotype.freq(child.tot[,c(-1,-2), drop=F], k)
        Avg.Hetero.Gen[g+1] <- 1 - Avg.Homo.Gen[g+1]
        
        ## Create parental pool for the next generation
        par.pool.f <- create.par.pool(child.f, s1g, s2g, N.pop, 
                                      k, allele.freq, loci.par, pld)
        par.pool.m <- create.par.pool(child.m, s1g, s2g, N.pop, 
                                      k, allele.freq, loci.par, pld)
        par.pool.tot <- rbind(par.pool.f, par.pool.m)
        
        colnames(par.pool.f) <- c('am', 'tr', c(sapply(1:pld, 
                                                       function(x) paste(x, c(1:k), sep='.'))))
        colnames(par.pool.m) <- c('am', 'tr', c(sapply(1:pld, 
                                                       function(x) paste(x, c(1:k), sep='.'))))
        colnames(par.pool.tot) <- c('am', 'tr', c(sapply(1:pld, 
                                                         function(x) paste(x, c(1:k), sep='.'))))
        
    }
    
    
    result <- matrix(NA, nrow=12, ncol=length(G)+1)
    result[1,] <- Var.HA
    result[2,] <- E.HA
    result[3,] <- Var.Tr
    result[4,] <- E.Tr
    result[5,] <- Cov.HT
    result[6,] <- Cor.HT
    result[7,] <- Cor.mate.HA
    result[8,] <- Cor.mate.Tr
    result[9,] <- Avg.Homo.AL
    result[10,] <- Avg.Hetero.AL
    result[11,] <- Avg.Homo.Gen
    result[12,] <- Avg.Hetero.Gen
    
    rownames(result) <- c('Var.HA', 'E.HA', 'Var.Tr', 'E.Tr', 
                          'Cov.HT', 'Cor.HT', 'Cor.mate.HA', 'Cor.mate.Tr',
                          'Avg.Homo.AL', 'Avg.Hetero.AL', 
                          'Avg.Homo.Gen', 'Avg.Hetero.Gen')
    
    return(result)

}



wf.sim.norm <- function(m.mode, N.pop, N.sim, G, k, z, s10, s1g=0, s2g=0,
                       am.const=0.5, x.prob=1, test.mode=T, pld=2,
                       rand.seed=1, n.core=4){
    
    source.folder <- '~/Desktop/Assortative_Mating/R/version_5'
    source(file.path(source.folder, 'init_wf.R'))
    if (!(m.mode %in% c('random','gen','phen')))
        stop('The am.mode should be one of random, gen, phen')
    
    
    ## =========== Initialization ===========
    L <- 1:k  # L: a set of pointers to the loci contributing to the trait
    X <- generate.x.z(k, z=z, r.seed=rand.seed)
    allele.freq <- generate.pq(k, r.seed=rand.seed, test.mode=test.mode) 
    
    print(X)
    # write allele frequencies and Xi's of the current run
    curr.t <- format(Sys.time(), "%Y_%m_%d_%H%M%S")
    x.save.f <- paste0("x_", curr.t, ".txt")
    al.frq.save.f <- paste0("al_frq_", curr.t, ".txt")
    
    write.table(X, file=x.save.f, quote=F, row.names=F, col.names=F)
    write.table(allele.freq, file=al.frq.save.f, quote=F, row.names=F) 
    
    
    ## Need to verify why this step was necessary
    #  ===========
    if (sum(X) != k) {
        # if x.prob != 1, pick X_i's so that E[T|X,M=A] >= E[T|X,M=B]
        # by flipping 0 and 1's if E[T|X,M=A] < E[T|X,M=B]
        X.tf <- as.logical(X)
        ET.A <- sum(allele.freq$p[X.tf]) + sum(1 - allele.freq$p[!X.tf])
        ET.B <- sum(allele.freq$q[X.tf]) + sum(1 - allele.freq$q[!X.tf])
        if (ET.A < ET.B) {
            X <- 1 - X
        }
    }
    # ===========
    
    ## Partition Loci based on the X values
    loci.par <- split(L, X) #partition of loci based on xi value
    
    if (is.null(loci.par[['1']]))
        loci.par['1'] <- list(NULL)
    if (is.null(loci.par[['0']]))
        loci.par['0'] <- list(NULL)
    
    # extend the loci partition to accomodate multiple ploidy
    if (pld != 1) {
        loci.par.1 <- loci.par$'1'
        loci.par.0 <- loci.par$'0'
        
        for (p in 1:(pld-1)) {
            loci.par$'1' <- c(loci.par$'1', (loci.par.1 + (k*p)))
            loci.par$'0' <- c(loci.par$'0', (loci.par.0 + (k*p)))
        }
    }
    
    
    
    ## =========== Do computation for N.sim times using parallel ========
    cl <- makeCluster(n.core, oufile="")
    registerDoParallel(cl)
    # registerDoSNOW(cl)
    # pb <- txtProgressBar(0, N.sim, style=3)
    
    result <- foreach(i=1:N.sim, .combine=rbind, .verbose=F,
                      .export='wf.single.norm') %dopar% {
        source(file.path(source.folder, 'init_wf.R'))
        source(file.path(source.folder, 'matrix_scaling.R'))
        
        # sink(file=log, append=TRUE, type="message")
        # cat(paste("Starting iteration", i, "\n"))
        wf.single.norm(m.mode=m.mode, N.pop=N.pop, G=G, k=k, s10=s10, 
                       s1g=s1g, s2g=s2g, am.const=am.const, pld=pld,
                       loci.par=loci.par, allele.freq=allele.freq)
    }

    stopCluster(cl)

    return(result)
}

