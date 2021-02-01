# Modified on Jan 14, 2019
# Located in ~/Desktop/Assortative_Mating/R/version_5/ folder

######################################################
#Simulations under Nicholson, Smith et ... Donnelly (2002) J Royal Stat Soc B.
#Model is
#1) Generate ancestral allele frequencies pi_i according to neutral SFS.
#2) add "drift variables" that are normal(0, c_j*pi_i*(1-pi_i)), where 
#c_j is like a population-specific Fst. Drift variables are independent
#for each locus and each population.
#3) any allele frequencies that end up less than 0 or greater than 1 are set to 0 or 1.

#We are ignoring their SNP ascertainment step.

# n.loci <- 1000
# n.anc <- 20000 #Ancestral n, for generating ancestral allele frequencies.
# c_A <- 0.3 # "c' in pop A.
# c_B <- 0.3 # 'c' in pop B. c_A and c_B are related to Fst.

#generate ancestral allele frequences for n.loci loci under neutral sfs with Ne = n.anc.
gen.neut.sfs <- function(n.loci, n.anc){
    probs.uw <- 1/(1:(n.anc-1))
    probs <- probs.uw/sum(probs.uw)
    cdf <- cumsum(probs)
    percs <- runif(n.loci, 0 , 1)
    get.cdfind <- function(perc){
        (sum(cdf <= perc) + 1)/n.anc 
    }
    sapply(percs, get.cdfind)
}

#supply a vector of allele frequncies and a drift parameter c, and one realization
#of the post-drift allele frequencies is provided.
drift.afs <- function(afs, c){
    drift <- rnorm(length(afs), 0, sqrt(c*afs*(1-afs)))
    drifted <- afs + drift
    drifted[drifted < 0] <- 0
    drifted[drifted > 1] <- 1
    return(drifted)
}

#Simulates p and q by drawing ancestral frequencies from a neutral SFS
#then adding drift according to Nicholson et al. n.anc determines 
#allele frequency distribution in ancestral pop. c_A and c_B are 
#drift parameters. rm.fixed.same=TRUE removes loci where the same
#allele is fixed in both populations.
sim.pq <- function(n.loci, n.anc, c_A, c_B, rm.fixed.same=T){
    anc.af <- gen.neut.sfs(n.loci, n.anc)
    drift.A <- drift.afs(anc.af, c_A)
    drift.B <- drift.afs(anc.af, c_B)
    # in the current am.sim, p > q, might need to fix later to make it 
    # consistent with Doc's paper
    # deris1 <- drift.B > drift.A # <- uncomment this later
    deris1 <- drift.B < drift.A # <- comment this later
    p <- numeric(n.loci)
    q <- numeric(n.loci)
    p[deris1 == 1] <- drift.A[deris1 == 1]
    q[deris1 == 1] <- drift.B[deris1 == 1]
    p[deris1 == 0] <- 1 - drift.A[deris1 == 0]
    q[deris1 == 0] <- 1 - drift.B[deris1 == 0]
    random1 <- drift.B == drift.A
    x <- rbinom(sum(random1), 1, 1/2)
    #For loci where, after drift, allele freqs are equal, decide "1" allele randomly.
    p[random1] <- p[random1] + x * (1 - 2 * p[random1])
    q[random1] <- q[random1] + x * (1 - 2 * q[random1])
    if(rm.fixed.same){
        p <- p[random1 == 0]
        q <- q[random1 == 0]
    }
    return(cbind(p, q))
}