# Modified on Jan 14, 2019
# Located in ~/Desktop/Assortative_Mating/R/version_5/ folder

# This is only for the simulated pi, qi mode with other parameters being 
# the ones of the base case. 

rm(list=ls())

library(parallel)

run.sim <- function(N.pop, N.sim, G, k, z, s10, s1g=0, s2g=0,
                    am.const=0.5, x.prob=1, test.mode=T, pld=2,
                    rand.seed=1, n.core=8) {
    
    source.folder <- '~/Desktop/Assortative_Mating/R/version_5'
    source(file.path(source.folder, 'wf_sim_norm.R'))
    
    result <- list()
    
    # ========= random ==========
    print('random')
    print(Sys.time())
    t1 <- proc.time()
    result.rand <- wf.sim.norm(m.mode='random', N.pop=N.pop, N.sim=N.sim, G=G, k=k, z=z,
                               s10=s10, s1g=s1g, s2g=s2g, am.const=am.const, 
                               x.prob=x.prob, test.mode=test.mode, pld=pld,
                               rand.seed=rand.seed, n.core=n.core)
    print(proc.time() - t1)
    # print(Sys.time())
    
    # ========= genotype ==========
    print('genotype')
    print(Sys.time())
    t1 <- proc.time()
    result.gen <- wf.sim.norm(m.mode='gen', N.pop=N.pop, N.sim=N.sim, G=G, k=k, z=z,
                               s10=s10, s1g=s1g, s2g=s2g, am.const=am.const, 
                               x.prob=x.prob, test.mode=test.mode, pld=pld,
                               rand.seed=rand.seed, n.core=n.core)
    print(proc.time() - t1)
    # print(Sys.time())
    
    
    # ========= phenotype ==========
    print('phenotype')
    print(Sys.time())
    t1 <- proc.time()
    result.phen <- wf.sim.norm(m.mode='phen', N.pop=N.pop, N.sim=N.sim, G=G, k=k, z=z,
                               s10=s10, s1g=s1g, s2g=s2g, am.const=am.const, 
                               x.prob=x.prob, test.mode=test.mode, pld=pld,
                               rand.seed=rand.seed, n.core=n.core)
    print(proc.time() - t1)
    print(Sys.time())
    
    
    result$rand <- result.rand
    result$gen <- result.gen
    result$phen <- result.phen
    
    return(result)
}



# ================== Run simulation ==================
save.folder <- '~/Desktop/Assortative_Mating/data_ver_5_xi/'

if(!dir.exists(save.folder)) {
    dir.create(save.folder)
}

setwd(save.folder)

num.loci <- 10
r.seed.vec <- c(634, 879, 410, 745, 663, 613, 925, 15, 192, 978)

for (rand.seed in r.seed.vec) {
    
    print(paste("random seed:", rand.seed))
    
    # N.pop=1000; N.sim=100; G=40; k=num.loci; s10=0.5;
    # s1g=0.0; s2g=0.0; am.const=0.5; x.prob=0.5; test.mode=T; pld=2;
    
    N.pop=100; N.sim=10; G=40; k=num.loci; z=5; s10=0.5;
    s1g=0.0; s2g=0.0; am.const=0.5; x.prob=0.5; test.mode=T; pld=2;
    
    n.core=detectCores();
    
    if (test.mode == F) {
        print('Random pq generation')
    }
    if (x.prob != 1.0) {
        print('Randomizing X!')
    }
    cat(paste("\nk = ", k, '\n'))
    
    
    result <- run.sim(N.pop=N.pop, N.sim=N.sim, G=G, k=k, z=z, s10=s10,
                      s1g=s1g, s2g=s2g, am.const=am.const,
                      x.prob=x.prob, test.mode=test.mode, pld=pld,
                      rand.seed=rand.seed, n.core=n.core)
    
    
    
    # ================== Process the result and save ==================
    random <- list()
    gen <-  list()
    phen <- list()
    
    n.par <- 12
    
    # Var.HA, E.HA, Var.Tr, E.Tr, Cov.HT, Cor.HT, Cor.mate.HA, Cor.mate.Tr,
    # Avg.Homo.AL, Avg.Hetero.AL, Avg.Homo.Gen, Avg.Hetero.Gen, order
    for (i in 1:n.par) {
        random[[i]] <- matrix(NA, nrow=N.sim, ncol=G+1)
        gen[[i]] <- matrix(NA, nrow=N.sim, ncol=G+1)
        phen[[i]] <- matrix(NA, nrow=N.sim, ncol=G+1)
    }
    
    for (r in 1:n.par) {
        for (i in 1:N.sim) {
            random[[r]][i,] <- result$rand[n.par*(i-1)+r, ]
            gen[[r]][i,] <- result$gen[n.par*(i-1)+r, ]
            phen[[r]][i,] <- result$phen[n.par*(i-1)+r, ]
        }
    }
    
    
    for (i in 1:n.par) {
        random[[i+n.par]] <- colMeans(random[[i]])
        gen[[i+n.par]] <- colMeans(gen[[i]])
        phen[[i+n.par]] <- colMeans(phen[[i]])
    }
    

    names.list <- c('Var.HA', 'E.HA', 'Var.Tr', 'E.Tr', 'Cov.HT', 'Cor.HT',
                    'Cor.mate.HA', 'Cor.mate.Tr', 'Avg.Homo.AL', 'Avg.Hetero.AL',
                    'Avg.Homo.Gen', 'Avg.Hetero.Gen',
                    'Var.HA.Avg', 'E.HA.Avg', 'Var.Tr.Avg', 'E.Tr.Avg',
                    'Cov.HT.Avg', 'Cor.HT.Avg',
                    'Cor.mate.HA.Avg', 'Cor.mate.Tr.Avg',
                    'Avg.Homo.AL.Avg', 'Avg.Hetero.AL.Avg',
                    'Avg.Homo.Gen.Avg', 'Avg.Hetero.Gen.Avg')
    
    names(random) <- names.list
    names(gen) <- names.list
    names(phen) <- names.list
    
    save.name <- sprintf('npop_%i_nsim_%i_k_%i_s10_%1.2f_test_%s_xp_%1.1f_rseed_%i.RData', 
                         N.pop, N.sim, k, s10, test.mode, x.prob, rand.seed)
    save.list <- c('result', 'random', 'gen', 'phen',
                   'N.pop', 'N.sim', 'G', 'k', 'z', 's10', 's1g', 's2g', 'pld',
                   'am.const', 'x.prob', 'test.mode', 'rand.seed', 'n.core')
    save(list=save.list, file=paste(save.folder, save.name, sep=''))
    cat("\n data has been saved \n")
    
    
    # ================== Plot ==================
    
    plt.save.name.1 <- sprintf('npop_%i_nsim_%i_k_%i_s10_%1.2f_test_%s_xp_%1.2f_rseed_%i.png', 
                               N.pop, N.sim, k, s10, test.mode, x.prob, rand.seed)
    
    png(filename=paste(save.folder, plt.save.name.1, sep=''), 
        width=2000, height=900, pointsize=15)
    par(mfrow=c(2,5), lwd=2, cex=1.2, mar=c(6, 4, 2, 2))
    plot(0:G, random$E.HA.Avg, type='b', col='red', xlab='Generations', 
         ylab='Avg E[HA]', main='Mean E[HA]', ylim=c(0,1.0), pch=1)
    lines(0:G, gen$E.HA.Avg, type='b', col='blue', pch=2)
    lines(0:G, phen$E.HA.Avg, type='b', col='green', pch=3)
    
    plot(0:G, random$Var.HA.Avg, type='b', col='red', xlab='Generations', 
         ylab='Avg Var[HA]', main='Mean Var[HA]', ylim=c(0, max(gen$Var.HA.Avg)), pch=1)
    lines(0:G, gen$Var.HA.Avg, type='b', col='blue', pch=2)
    lines(0:G, phen$Var.HA.Avg, type='b', col='green', pch=3)
    
    plot(0:G, random$Cov.HT.Avg, type='b', col='red', xlab='Generations', 
         ylab='Avg Cov[H,T]', main='Mean Covariance', ylim=c(0, random$Cov.HT.Avg[1]), pch=1)
    lines(0:G, gen$Cov.HT.Avg, type='b', col='blue', pch=2)
    lines(0:G, phen$Cov.HT.Avg, type='b', col='green', pch=3)
    
    plot(1:G, random$Cor.mate.HA.Avg[-1], type='b', col='red', xlab='Generations', 
         ylab='Avg Cor.mate.HA', main='Mean Cor.mate.HA', 
         ylim=c(-0.2, 1), pch=1)
    lines(1:G, gen$Cor.mate.HA.Avg[-1], type='b', col='blue', pch=2)
    lines(1:G, phen$Cor.mate.HA.Avg[-1], type='b', col='green', pch=3)
    
    plot(0:G, random$Avg.Homo.AL.Avg, type='b', col='red', xlab='Generations', 
         ylab='Avg Avg Homozygosity', main='Mean Avg Homozygosity', 
         ylim=c(0.45,0.6), pch=1)
    lines(0:G, gen$Avg.Homo.AL.Avg, type='b', col='blue', pch=2)
    lines(0:G, phen$Avg.Homo.AL.Avg, type='b', col='green', pch=3)
    
    
    
    plot(0:G, random$E.Tr.Avg, type='b', col='red', xlab='Generations', 
         ylab='Avg E[Tr]', main='Mean E[Tr]', ylim=c(0,2*k), pch=1)
    lines(0:G, gen$E.Tr.Avg, type='b', col='blue', pch=2)
    lines(0:G, phen$E.Tr.Avg, type='b', col='green', pch=3)
    
    plot(0:G, random$Var.Tr.Avg, type='b', col='red', xlab='Generations', 
         ylab='Avg Var[Tr]', main='Mean Var[Tr]', 
         ylim=c(0, max(phen$Var.Tr.Avg)), pch=1)
    lines(0:G, gen$Var.Tr.Avg, type='b', col='blue', pch=2)
    lines(0:G, phen$Var.Tr.Avg, type='b', col='green', pch=3)
    
    plot(0:G, random$Cor.HT.Avg, type='b', col='red', xlab='Generations', 
         ylab='Avg Cor[H,T]', main='Mean Correlation', ylim=c(0,1), pch=1)
    lines(0:G, gen$Cor.HT.Avg, type='b', col='blue', pch=2)
    lines(0:G, phen$Cor.HT.Avg, type='b', col='green', pch=3)
    
    plot(1:G, random$Cor.mate.Tr.Avg[-1], type='b', col='red', xlab='Generations', 
         ylab='Avg Cor.mate.Tr', main='Mean Cor.mate.Tr', 
         ylim=c(-0.2, 1), pch=1)
    lines(1:G, gen$Cor.mate.Tr.Avg[-1], type='b', col='blue', pch=2)
    lines(1:G, phen$Cor.mate.Tr.Avg[-1], type='b', col='green', pch=3)
    
    plot(0:G, random$Avg.Homo.Gen.Avg, type='b', col='red', xlab='Generations', 
         ylab='Avg Avg Homozygote Freq', main='Mean Avg Homozygote Freq', 
         ylim=c(0,1), pch=1)
    lines(0:G, gen$Avg.Homo.Gen.Avg, type='b', col='blue', pch=2)
    lines(0:G, phen$Avg.Homo.Gen.Avg, type='b', col='green', pch=3)
    
    
    par(fig=c(0,1,0,1), oma=c(0,0,0,0), mar=c(0,0,0,0), new=T)
    plot(0,0, type='n', bty='n', xaxt='n', yaxt='n')
    legend('bottom', legend=c('rand', 'gen', 'phen'), xpd=T, horiz=T, inset=c(0,0),
           bty='n', col=c('red', 'blue', 'green'), lty=c(1,1,1), pch=c(1,2,3))
    
    dev.off()
    
    
    
    # ======= individual plots ========
    if (N.sim >= 100) {
        n.plot <- 100
    } else {
        n.plot <- N.sim
    }
    
    ind.plot <- sample(1:N.sim, n.plot)
    col.map <- rainbow(n.plot)
    
    plt.save.name.2 <- sprintf('npop_%i_nsim_%i_k_%i_s10_%1.2f_test_%s_xp_%1.2f_rseed_%i_2.png', 
                               N.pop, N.sim, k, s10, test.mode, x.prob, rand.seed)
    
    png(filename=paste(save.folder, plt.save.name.2, sep=''), 
        width=3000, height=1000, pointsize=15)
    
    par(mfrow=c(3,10), lwd=2, cex=1.2, mar=c(6, 4, 2, 2))
    
    # ========= random plots
    plot(0:G, random$E.HA[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='E[HA]', main='rand E[HA]', ylim=c(0,1))
    for (i in 2:n.plot) {
        lines(0:G, random$E.HA[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, random$Var.HA[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Var[HA]', main='rand Var[HA]', 
         ylim=c(0, random$Var.HA.Avg[1]))
    for (i in 2:n.plot) {
        lines(0:G, random$Var.HA[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, random$E.Tr[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='E[Tr]', main='rand E[Tr]', ylim=c(0,2*k))
    for (i in 2:n.plot) {
        lines(0:G, random$E.Tr[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, random$Var.Tr[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Var[Tr]', main='rand Var[Tr]', 
         ylim=c(0, max(random$Var.Tr.Avg)))
    for (i in 2:n.plot) {
        lines(0:G, random$Var.Tr[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, random$Cov.HT[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Cov[H,T]', main='rand Cov[H,T]', 
         ylim=c(0, random$Cov.HT[1]))
    for (i in 2:n.plot) {
        lines(0:G, random$Cov.HT[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, random$Cor.HT[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Cor[H,T]', main='rand Cor[H,T]', 
         ylim=c(0, 1), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, random$Cor.HT[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, random$Cor.mate.HA[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Cor[mates HA]', main='rand Cor[mates HA]', 
         ylim=c(0, 1), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, random$Cor.mate.HA[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, random$Cor.mate.Tr[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Cor[mates Tr]', main='rand Cor[mates Tr]', 
         ylim=c(0, 1), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, random$Cor.mate.Tr[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, random$Avg.Homo.AL[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Avg.Homo.AL', main='rand Avg.Homo.AL', 
         ylim=c(0, 1), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, random$Avg.Homo.AL[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, random$Avg.Homo.Gen[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Avg.Homo.AL', main='rand Avg.Homo.Gen', 
         ylim=c(0, 1), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, random$Avg.Homo.Gen[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    
    # ======== genotype plots
    plot(0:G, gen$E.HA[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='E[HA]', main='gen E[HA]', ylim=c(0,1), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, gen$E.HA[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, gen$Var.HA[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Var[HA]', main='gen Var[HA]', 
         ylim=c(0, max(gen$Var.HA.Avg)), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, gen$Var.HA[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, gen$E.Tr[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='E[Tr]', main='gen E[Tr]', ylim=c(0,2*k), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, gen$E.Tr[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, gen$Var.Tr[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Var[Tr]', main='gen Var[Tr]', 
         ylim=c(0, max(gen$Var.Tr.Avg)), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, gen$Var.Tr[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, gen$Cov.HT[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Cov[H,T]', main='gen Cov[H,T]', 
         ylim=c(0, gen$Cov.HT[1]), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, gen$Cov.HT[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, gen$Cor.HT[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Cor[H,T]', main='gen Cor[H,T]', 
         ylim=c(0, 1), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, gen$Cor.HT[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, gen$Cor.mate.HA[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Cor[mates HA]', main='gen Cor[mates HA]', 
         ylim=c(0, 1), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, gen$Cor.mate.HA[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, gen$Cor.mate.Tr[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Cor[mates Tr]', main='gen Cor[mates Tr]', 
         ylim=c(0, 1), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, gen$Cor.mate.Tr[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, gen$Avg.Homo.AL[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Avg.Homo.AL', main='gen Avg.Homo.AL', 
         ylim=c(0, 1), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, gen$Avg.Homo.AL[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, gen$Avg.Homo.Gen[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Avg.Homo.AL', main='gen Avg.Homo.Gen', 
         ylim=c(0, 1), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, gen$Avg.Homo.Gen[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    
    
    
    # ======== phenotype plots
    plot(0:G, phen$E.HA[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='E[HA]', main='phen E[HA]', ylim=c(0,1), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, phen$E.HA[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, phen$Var.HA[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Var[HA]', main='phen Var[HA]', 
         ylim=c(0, phen$Var.HA.Avg[1]), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, phen$Var.HA[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, phen$E.Tr[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='E[Tr]', main='phen E[Tr]', ylim=c(0,2*k), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, phen$E.Tr[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, phen$Var.Tr[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Var[Tr]', main='phen Var[Tr]', 
         ylim=c(0, max(phen$Var.Tr.Avg)), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, phen$Var.Tr[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, phen$Cov.HT[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Cov[H,T]', main='phen Cov[H,T]', 
         ylim=c(0, phen$Cov.HT[1]), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, phen$Cov.HT[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, phen$Cor.HT[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Cor[H,T]', main='phen Cor[H,T]', 
         ylim=c(0, 1), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, phen$Cor.HT[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, phen$Cor.mate.HA[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Cor[mates HA]', main='phen Cor[mates HA]', 
         ylim=c(0, 1), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, phen$Cor.mate.HA[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, phen$Cor.mate.Tr[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Cor[mates Tr]', main='phen Cor[mates Tr]', 
         ylim=c(0, 1), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, phen$Cor.mate.Tr[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, phen$Avg.Homo.AL[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Avg.Homo.AL', main='phen Avg.Homo.AL', 
         ylim=c(0, 1), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, phen$Avg.Homo.AL[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    plot(0:G, phen$Avg.Homo.Gen[ind.plot[1], ], type='l', col=col.map[1], 
         xlab='Generations', ylab='Avg.Homo.AL', main='phen Avg.Homo.Gen', 
         ylim=c(0, 1), pch=1)
    for (i in 2:n.plot) {
        lines(0:G, phen$Avg.Homo.Gen[ind.plot[i], ], type='l', col=col.map[i])
    }
    
    
    
    
    par(fig=c(0,1,0,1), oma=c(0,0,0,0), mar=c(0,0,0,0), new=T)
    plot(0,0, type='n', bty='n', xaxt='n', yaxt='n')
    legend('bottom', legend=sprintf('npop=%i, nsim=%i, k=%i, s10=%1.2f, test=%s, xp=%1.1f', 
                                    N.pop, N.sim, k, s10, test.mode, x.prob))
    
    dev.off()
    
}













