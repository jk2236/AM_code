# Modified on Jan 14, 2019
# This is the version for randomizing X_i's
# Located in ~/Desktop/Assortative_Mating/R/version_5/ folder

# This is a collection of library functions for simulating 
# an assortative mating in diploids. 
# for diploid with sex simulation 

library(data.table)
source('~/Desktop/Assortative_Mating/R/version_5/gen_pq_doc.R') #for gen.pq in random pq generation mode

rand <- function(r.seed) {
    # set random seed for generate functions
    set.seed(r.seed)
}

# generate.x <- function(k, p.x, r.seed) {
#     # Generate random X values for the trait loci, either 0 or 1
#     # k: number of trait loci
#     # p.x: probability for a trait loci to have value 1
#     rand(r.seed)
#     return(c(rbinom(k, 1, p.x)))
# }


generate.x.z <- function(k, z, r.seed) {
    # This function generates X values under prescribed sum.
    # Each locus has equal probability of being 1 or 0 (i.e., neutral model).
    # 
    # k: number of loci
    # z: sum of X values across all k loci
    # r.seed: random seed 
    rand(r.seed)
    x.z <- rep(0, k)
    x.z[sample.int(k, size=z)] <- 1
    return(x.z)
}

generate.pq <- function(k, r.seed, test.mode=T, rand.mode='sim') {
    # modified just for X_i run
    
    # Generate random L values for the trait loci. The values are sampled
    # from the uniform distribution [0,1] and p >= q
    # If test.mode=F, there're two possible modes to generate the pq frequencies
    # rand.mode='unif': draw p,q from uniform[0,1] distribution
    # rand.mode='sim': draw p,q from simulated random drift 
    rand(r.seed)
    
    if (test.mode) {
        prob <- data.frame(p=rep(1.0, k), q=rep(0.0, k))
    } else {
            stop('unsupported rand.mode in generate.pq')
    }
    
    return(prob)
}

generate.comb <- function(seq.1, seq.2) {
    # this function generates all possible combinations of seq1 and seq2 
    # with repetition. Not being used in the current version of simulation code.
    cbind(rep.int(seq.1, length(seq.2)),
          c(t(matrix(rep.int(seq.2, length(seq.1)), nrow=length(seq.2)))))
}

assort.func <-function(assort.const, t1, t2) {
    am <- exp(-assort.const * abs(t1-t2))
    # am <- 1 / (1 + assort.const * abs(t1-t2))
    return(am)
}


gen.pop.sample <- function(group.id, n.sample, k, allele.freq, ploidy=2) {
    # This function generates n number of random population samples of 
    # a particular group. Input:
    # group.id: either 'A' or 'B'
    # n.sample: number of random samples to generate
    # k: number of trait loci
    # allele.freq: allele frequency at each loci
    # ploidy: default is 2 for diploid 
    
    # The output is a matrix format in the following order:
    # col.1: admixture fraction
    # col.2: NA (col reserved for trait value)
    # col.3~(ploidy*k+2): genotype of a given individual 
    
    pop.sample <- matrix(NA, nrow=n.sample, ncol=ploidy*k+2)
    if(group.id == 'A') {
        pop.sample[,1] <- rep(1, n.sample)
        pop.sample[,3:(ploidy*k+2)] <- sapply(1:(ploidy*k), 
                                              function(x) 
                                                  rbinom(n.sample, 1, 
                                                         rep(allele.freq$p, ploidy)[x]))
    } else if (group.id == 'B') {
        pop.sample[,1] <- rep(0, n.sample)
        # pair.IL <- expand.grid(i=1:n.sample, l=1:k)
        pop.sample[,3:(ploidy*k+2)] <- sapply(1:(ploidy*k), 
                                              function(x) 
                                                  rbinom(n.sample, 1, 
                                                         rep(allele.freq$q, ploidy)[x]))
    } else {
        stop('unsupported group.id')
    }
    return(pop.sample)
}


eval.trait <- function(genotype, loci.par) {
    # This function computes trait from the given genotype and a given
    # X_i values (loci.par) across the trait loci.
#     if (!is.matrix(genotype)) {
#         genotype <- matrix(genotype)
#     }
    return(rowSums(genotype[ ,loci.par$`1`, drop=F]) 
           + rowSums(1 - genotype[ ,loci.par$`0`, drop=F]))
}




eval.allele.freq <- function(genotype, k) {
    # This function computes allele frequencies in the population.
    # It assumes diploid and biallelic case with allele "1" and allele "0".
    
    # Input
    # genotype: N.pop x 2k matrix with col.1-k first genotype
    #           and col.(k+1)-2k second genotype
    # k: number of loci
    
    # Output
    # Allele frequency for allele "1" and "0" for each loci
    # The output is in 2 x k matrix with the first row for "1" allele
    # and the second row for "0" allele frequencies.
    
    if (dim(genotype)[2] != (2*k))
        stop("check input genotype format")
    
    N <- dim(genotype)[1] # number of individuals in the population
    al.freq.1 <- (colSums(genotype[,1:k, drop=F]) + 
                      colSums(genotype[,(k+1):(2*k), drop=F])) / (2*N)
    al.freq.0 <- 1 - al.freq.1
    al.freq <- data.frame('A1'=al.freq.1, 'A0'=al.freq.0, row.names=1:k)
    
    return(al.freq)
}



eval.genotype.freq <- function(genotype, k) {
    # This function computes homozygote and heterozygote genotype frequencies 
    # in the population. It assumes diploid and biallei with allele "1" and "0".
    
    # Input
    # genotype: N.pop x 2k matrix with col.1-k first genoytpe
    #           and col.(k+1)-2k second genotype
    # k: number of loci
    
    # Output
    # Homozygote genotype frequencies average over all k loci. 
    
    if (dim(genotype)[2] != (2*k))
        stop("check input genotype format")
    
    N <- dim(genotype)[1] # number of individuals in the population
    homo.freq <- colSums(genotype[,1:k, drop=F] == 
                             genotype[,(k+1):(2*k), drop=F]) / N
    avg.homo.freq <- mean(homo.freq)
    
    return(avg.homo.freq)
    
}


eval.avg.homo <- function(genotype, k) {
    # This function computes average homozygosity 
    # across all k loci. Considering biallelic case, let A_1,i and A_0,i be
    # the frequency of "1" allele and "0" allele at i-th locus, repectively.
    # Avg.Homo = (1/k)*sum(i=1 to k)(A_1,i^2 + A_0,i^2)
    # Avg.Hetero = (1/k)*sum(i=1 to k)(1 - A_1,i^2 - A_0,i^2)
    #            = 1 - Avg.Homo
    
    # Input:
    # genotype: N.pop x 2k matrix with col.1-k first genotype
    #           and col.(k+1)-2k second genotype
    # k: number of loci
    
    # Output:
    # Avg.Homo: average homozygosity
    
    al.freq <- eval.allele.freq(genotype, k)
    
    Avg.Homo <- mean(rowSums(al.freq^2))
    return(Avg.Homo)
}



create.par.pool <- function(child, s1g, s2g, n.pop, 
                            k, allele.freq, loci.par, pld=2) {
    # create parental pool from combining source populations and the
    # hybrid population. 
    # child: matrix of the child created from the mating
    # n.pop: population size, constant across all generations
    # k: number of trait loci
    # allele.freq: allele frequency at each loci
    # loci.par: X_i values across the trait loci
    # pld: ploidy (default for diploids)
    
    # The input child and all the output is the following format:
    # col.1: admixture fraction
    # col.2: NA (col reserved for trait value)
    # col.3~(k+2): genotype of a given individual 
    
    n.1 <- round(n.pop * s1g) #number of people from population A
    n.2 <- round(n.pop * s2g) #number of people from population B
    n.h <- n.pop - n.1 - n.2 #number of hybrid in the new parental population
    
    p.pool <- matrix(NA, nrow=n.pop, ncol=pld*k+2)
    if (n.1 != 0) {
        p.pool[1:n.1,] <- gen.pop.sample('A', n.1, k, allele.freq, ploidy=pld)
    }
    if (n.2 != 0) {
        p.pool[(n.1+1):(n.1+n.2),] <- gen.pop.sample('B', n.2, k, 
                                                     allele.freq, ploidy=pld)
    }
    if (n.h != 0) {
        p.pool[(n.1+n.2+1):n.pop,] <- child[sample.int(n=n.pop, size=n.h),]
    }
    
    p.pool[,2] <- eval.trait(p.pool[ , c(-1,-2), drop=F], loci.par)
    
    if (any(is.na(p.pool))) {
        stop('error occurred during the parental pool generation')
    }
    
    return(p.pool)
}



sample.mates <- function(N.pop, N.sample) {
    # This function samples mates for random mating scheme.
    # For a given N male and N female, it samples, it randomly generates
    # N pairs of male and female with or without replacement.
    mate.ind <- cbind(sample.int(N.pop, N.sample, replace=T), 
                      sample.int(N.pop, N.sample, replace=T))
    return(mate.ind)
}


gen.offspring.norm <- function(N.pop, mate.data.f, mate.data.m, am.const) {
    # generate all possible mate combinations for normalization
    # This function was modified for the new version including male and female
    # specific populations. Sept 6, 2017
    
    mate.comb <- CJ(1:N.pop, 1:N.pop) 
    names(mate.comb) <- c('Mf', 'Mm')
    mate.comb.1 <- mate.comb$Mf
    mate.comb.2 <- mate.comb$Mm
    
    p.mate <- assort.func(am.const, mate.data.f[mate.comb.1], mate.data.m[mate.comb.2])
    p.mate.mat <- matrix(p.mate, nrow=N.pop, ncol=N.pop)
    r.sum <- rep(1, N.pop)
    p.mate.mat.norm <- matrix.scale(p.mate.mat, r.sum, r.sum)
    mate.pick <- sample.int(N.pop^2, size=N.pop, replace=T, prob=c(p.mate.mat.norm))
    # sample with replacement or sample without replacement?
    mate.ind <- mate.comb[mate.pick, ]
    
    return(mate.ind)
}


# gen.offspring <- function(N.pop, mate.data, am.const) {
#     counter <- 0
#     mate.ind <- matrix(NA, nrow=N.pop, ncol=2)
#     while(TRUE) {
#         # print(counter)
#         N.sample <- N.pop
#         mate.try <- sample.mates(N.pop=N.pop, N.sample=N.sample)
#         p.mate <- assort.func(am.const, mate.data[mate.try[,1]], mate.data[mate.try[,2]])
#         p.cut <- runif(N.sample)
#         success <- (p.mate >= p.cut)
#         mate.result <- mate.try[success,]
#         n.success <- sum(success)
#         # print(n.success)
#         # print((counter + n.success) >= N.pop)
#         if ((counter + n.success) >= N.pop) {
#             mate.ind[(counter+1):N.pop,] <- mate.result[1:(N.pop-counter),]
#             break   
#         } else {
#             mate.ind[(counter+1):(counter+n.success),] <- mate.result
#             counter <- counter + n.success
#         }
#     }
#     return(mate.ind)
# }





error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    # this function draws error bar on a existing plot. 
    # x, y: original data points in the existing plot (mean value)
    # upper: standard deviation
    
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
        stop("vectors must be same length")
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


# combine.single <- function(x1, x2) {
#     return(rbindlist(list(x1, x2), use.names=T))
# }


cal.fst <- function(prob, ploidy=1) {
    # This function calculates multilocus Fst in two alleles & two populations model.
    # prob: data.frame with $p and $q column. Each row represents different locus.
    k <- dim(prob)[1]
    del.mean <- mean(prob$q - prob$p)
    s.del.sq <- sum((prob$q - prob$p - del.mean)^2) / k
    del.sq.mean <- sum((prob$q - prob$p)^2) / k
    
    p.mean <- mean(prob$p)
    q.mean <- mean(prob$q)
    
    s.p.sq <- sum((prob$p - p.mean)^2) / k
    s.q.sq <- sum((prob$q - q.mean)^2) / k
    
    Fst.k <- ploidy*del.sq.mean / (2*(p.mean*(1-p.mean) - s.p.sq +
                                   q.mean *(1-q.mean ) - s.q.sq) + ploidy*del.sq.mean)
    return(Fst.k)
    
}




