library(ggplot2)
library(assertr)
library(data.table)
library(reshape2)
library(dplyr)

# assign nucleotide biases, based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2648205/
nuc.prob <- c(19.745, 19.745, 30.255, 30.255) / 100
read.length <- 200
base.vec <- c("A", "T", "G", "C")
# set static parameters
n.sim <- 100
#read.length <- 1000

read.lengths <- c(100, 200, 300, 400, 600, 800, 1000)
# generate all mutation rates and kmer values
mu.rates <- seq(0, 20, by=1) / 100
kmer.sizes <- seq(3, 31, by=2)

# generate genomes to sample from
num.target <- 100
length.target <- 20000
num.contaminent <- 1
length.contaminent <- 2200000
generate.contaminant <- FALSE
generate.target <- FALSE


target.mat <- matrix(nrow=num.target, ncol=length.target)
contam.mat <- matrix(nrow=num.contaminent, ncol=length.contaminent)

for (i in 1:num.target)
{
  target.mat[i,] <- sample(base.vec, length.target, prob = nuc.prob, replace = TRUE)
}

# write target sequences to file
for (i in 1:num.contaminent)
{
  contam.mat[i,] <- sample(base.vec, length.contaminent, prob = nuc.prob, replace = TRUE)
}


# only use when generating new k-mer hashes
if (generate.target == TRUE)
{
  write.csv(target.mat, paste("sim_target_seq.csv", sep = ""), row.names=FALSE)
  for (i in 1:length(kmer.sizes)) 
  {
    kmer.length <- kmer.sizes[[i]]
    
    if (nrow(target.mat) == 1)
    {
      data_list <- unique(sapply(seq(from=1, to=(ncol(target.mat) - kmer.length + 1)), function(j) paste(target.mat[, j:(j + kmer.length - 1)], collapse="")))
    } else {
      data_list <- unique(c(unlist(sapply(seq(from=1, to=(ncol(target.mat) - kmer.length + 1)), function(j) unique(apply(target.mat[, j:(j + kmer.length - 1)], 1, paste, collapse=""))))))
    }
    
    present <- data.frame(x = data_list)
    
    write.csv(present, paste("sim_target_kmers_", kmer.length, ".csv", sep = ""), row.names=FALSE)
  }
}

# only use when generating new k-mer hashes
if (generate.contaminant == TRUE)
{
  for (i in 1:length(kmer.sizes)) 
  {
    kmer.length <- kmer.sizes[[i]]
    
    if (nrow(contam.mat) == 1)
    {
      data_list <- unique(c(sapply(seq(from=1, to=(ncol(contam.mat) - kmer.length + 1)), function(j) paste(contam.mat[, j:(j + kmer.length - 1)], collapse=""))))
    } else {
      data_list <- unique(c(sapply(seq(from=1, to=(ncol(contam.mat) - kmer.length + 1)), function(j) unique(apply(contam.mat[, j:(j + kmer.length - 1)], 1, paste, collapse="")))))
    }
    
    present <- data.frame(x = data_list)
    
    write.csv(present, paste("sim_contam_kmers_", kmer.length, ".csv", sep = ""), row.names=FALSE)
  }
}


target.mat <- as.matrix(read.csv(file = paste("sim_target_seq.csv", sep = ""), header = TRUE))

#read.lengths <- c(200)
for (ri in 1:length(read.lengths))
{
  read.length <- read.lengths[[ri]]
  
  # generate dataframe to hold values
  all.sim <- data.frame(kmer = c(), mu = c(), prop.kmers = c(), prop.fp = c())
  # ref <- rep(TRUE, read.length)
  
  # generate reference sequences for each simulation
  ref.mat <- matrix(nrow=n.sim, ncol=read.length)
  nonref.mat <- matrix(nrow=n.sim, ncol=read.length)
  
  for (sim in 1:n.sim)
  {
    # get random slice from targets
    ref.seq <- sample(1:num.target, 1)
    ref.pos <- sample(1:(length.target - read.length), 1)
    ref <- target.mat[ref.seq, ref.pos:(ref.pos + read.length - 1)]
    ref.mat[sim,] <- ref
    
    # repeat for contaminents
    ref.seq <- sample(1:num.contaminent, 1)
    ref.pos <- sample(1:(length.contaminent - read.length), 1)
    ref <- contam.mat[ref.seq, ref.pos:(ref.pos + read.length - 1)]
    nonref.mat[sim,] <- ref
  }
  
  query.mat <- array(dim=c(n.sim, read.length, length(mu.rates)))
  nontarget.mat <- array(dim=c(n.sim, read.length, length(mu.rates)))
  
  for (mi in 1:length(mu.rates))
  {
    mur <- mu.rates[[mi]]
    query.mat[,,mi] <- ref.mat
    nontarget.mat[,,mi] <- nonref.mat
    
    # run for query
    if (mur == 0)
    {
      next
    } else {
      nmu <- rpois(n=n.sim, lambda = (read.length * mur))
    }
    
    for (sim in 1:n.sim)
    {
      pos.mu <- sort(sample(seq(1, read.length), nmu[sim], replace = FALSE))
      to_mutate <- ref.mat[sim, pos.mu]
      to_mutate[ref.mat[sim, pos.mu] == "A"] <- sample(base.vec[base.vec != "A"], 1, replace = TRUE)
      to_mutate[ref.mat[sim, pos.mu]== "T"] <- sample(base.vec[base.vec != "T"], 1, replace = TRUE)
      to_mutate[ref.mat[sim, pos.mu] == "G"] <- sample(base.vec[base.vec != "G"], 1, replace = TRUE)
      to_mutate[ref.mat[sim, pos.mu] == "C"] <- sample(base.vec[base.vec != "C"], 1, replace = TRUE)
      query.mat[sim, pos.mu, mi] <- to_mutate
    }
    
    # and for contaminent
    nmu <- rpois(n=n.sim, lambda = (read.length * mur))
    for (sim in 1:n.sim)
    {
      pos.mu <- sort(sample(seq(1, read.length), nmu[sim], replace = FALSE))
      to_mutate <- nonref.mat[sim, pos.mu]
      to_mutate[nonref.mat[sim, pos.mu] == "A"] <- sample(base.vec[base.vec != "A"], 1, replace = TRUE)
      to_mutate[nonref.mat[sim, pos.mu]== "T"] <- sample(base.vec[base.vec != "T"], 1, replace = TRUE)
      to_mutate[nonref.mat[sim, pos.mu] == "G"] <- sample(base.vec[base.vec != "G"], 1, replace = TRUE)
      to_mutate[nonref.mat[sim, pos.mu] == "C"] <- sample(base.vec[base.vec != "C"], 1, replace = TRUE)
      nontarget.mat[sim, pos.mu, mi] <- to_mutate
    }
    
  }
    
  for (k in 1:length(kmer.sizes))
  {
    kmer.length <- kmer.sizes[[k]]
    
    # get k-mers of contaminating sequences
    k.mat <- read.csv(file = paste("sim_target_kmers_", kmer.length, ".csv", sep = ""), header = TRUE)
    target.kmers <- setkey(as.data.table(k.mat))
    num.kmers <- (read.length - kmer.length + 1)
    
    query.slice.list <- lapply(seq_len(num.kmers), function(i) query.mat[, i:(i + kmer.length - 1),])
    j <- 1
    match.count <- sapply(seq_len(num.kmers), function(i) sapply(seq_len(length(mu.rates)), function(j) col_concat(query.slice.list[[i]][,,j]) %in% as.matrix(target.kmers[.(col_concat(query.slice.list[[i]][,,j])), nomatch = 0L])), simplify = 'array')
    num.matches <- sapply(seq_len(length(mu.rates)), function(j) rowSums(match.count[,j,]))
    #ref.slice.list <- lapply(seq_len(num.kmers), function(i) ref.mat[, i:(i + kmer.length - 1)])
    #match.list <- lapply(seq_len(num.kmers), function(i) lapply(seq_len(length(mu.rates)), function(j) query.slice.list[[i]][,,j] == ref.slice.list[[i]]))
    #list.match <- sapply(seq_len(num.kmers), function(i) sapply(seq_len(length(mu.rates)), function(j) rowSums(match.list[[i]][[j]]) == kmer.length), simplify = 'array')
    #num.matches <- sapply(seq_len(length(mu.rates)), function(j) rowSums(list.match[,j,]))

    # check for reads from contaminent
    contam.slice.list <- lapply(seq_len(num.kmers), function(i) nontarget.mat[, i:(i + kmer.length - 1),])
    false.match.count <- sapply(seq_len(num.kmers), function(i) sapply(seq_len(length(mu.rates)), function(j) col_concat(contam.slice.list[[i]][,,j]) %in% as.matrix(target.kmers[.(col_concat(contam.slice.list[[i]][,,j])), nomatch = 0L])), simplify = 'array')
    num.false.pos <- sapply(seq_len(length(mu.rates)), function(j) rowSums(false.match.count[,j,]))
    # query.slice.list <- lapply(seq_len(num.kmers), function(i) nonref.mat[, i:(i + kmer.length - 1)])
    # false.match.count <- sapply(seq_len(num.kmers), function(i) col_concat(query.slice.list[[i]]) %in% as.matrix(target.kmers[.(col_concat(query.slice.list[[i]])), nomatch = 0L]))
    # num.false.pos <- rowSums(false.match.count)
    # match.list <- lapply(seq_len(num.kmers), function(i) query.slice.list[[i]] == ref.slice.list[[i]])
    # list.match <- sapply(seq_len(num.kmers), function(i) rowSums(match.list[[i]]) == kmer.length)
    # num.false.pos <- rowSums(list.match)
    
    # stack columns
    colnames(num.matches) <- mu.rates
    colnames(num.false.pos) <- mu.rates
    
    num.matches <- melt(num.matches, id.vars=mu.rates,
                     value.name = "Matches")
    colnames(num.matches) <- c("Var1", "Mu", "Matches")
    
    num.false.pos <- melt(num.false.pos, id.vars=mu.rates,
                        value.name = "FP")
    colnames(num.false.pos) <- c("Var1", "Mu", "FP")
    concat.mat <- merge(num.matches, num.false.pos, by=c("Var1","Mu"))
    
    # update dataframe
    to.append <- data.frame(kmer = kmer.length, mu = concat.mat$Mu, prop.kmers = (concat.mat$Matches / num.kmers), prop.fp = (concat.mat$FP / num.kmers))
    all.sim <- rbind(all.sim, to.append)
  }
  
  
  
  # calculate overall matching for mash
  #tp <- subset(all.sim, type == "TP")
  num.kmers <- (read.length - all.sim$kmer + 1)
  match.kmers <- round(num.kmers * all.sim$prop.kmers)
  mismatch.kmers <- num.kmers - match.kmers
  jaccard <- match.kmers / (match.kmers + (2 * mismatch.kmers))
  mash.tp <- ifelse(jaccard == 0, 0, 1 - ((-1/all.sim$kmer) * log((2 * jaccard) / (1 + jaccard))))
  
  #to.append <- data.frame(kmer = tp$kmer, type = "TP.mash", mu = tp$mu, prop.kmers = mash.tp)
  all.sim$mash <- mash.tp
  
  # calculate false positive mash
  #fp <- subset(all.sim, type == "FP")
  match.kmers <- round(num.kmers * all.sim$prop.fp)
  mismatch.kmers <- num.kmers - match.kmers
  jaccard <- match.kmers / (match.kmers + (2 * mismatch.kmers))
  mash.fp <- ifelse(jaccard == 0, 0, 1 - ((-1/all.sim$kmer) * log((2 * jaccard) / (1 + jaccard))))
  # to.append <- data.frame(kmer = tp$kmer, type = "FP.mash", mu = tp$mu, prop.kmers = mash.fp)
  # all.sim <- rbind(all.sim, to.append)
  all.sim$mash.fp <- mash.fp
  
  
  p <- ggplot(all.sim, aes(x=as.factor(mu), y=prop.kmers)) + facet_grid(~ kmer) + geom_boxplot() + scale_x_discrete(breaks = function(x){x[c(TRUE, FALSE)]}) + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.position="none") + ylab("Proportion matching k-mers") + xlab("Mutation rate") + scale_y_continuous(limits=c(NA, 1))
  p
  
  ggsave(paste("kmer_mu_comp_readlen", read.length, ".png", sep=""), width=8000, height=2000, units="px")
  
  # subset to sensible kmer sizes
  
  # calculate proportions of false positives
  #all.prob <- data.frame(kmer = c(), no.matches = c())
  
  # based on poppunk paper, k = (log(l) + log(1 - p) - log(p))/log(4) so p = L/((4^k) + L)
  
  # for (i in 1:length(kmer.sizes))
  # {
  #   kmer.length <- kmer.sizes[[i]]
  #   for (j in 1:n.sim)
  #   {
  #     kmer.prob <- 1/4**kmer.length
  #     dist <- sample(seq(0, 1), (read.length - kmer.length + 1), prob = c((1 - kmer.prob), kmer.prob), replace = TRUE)
  #     to.append <- data.frame(kmer = kmer.length, no.matches = sum(dist), prop.kmers = sum(dist) / (read.length - kmer.length + 1))
  #     
  #     all.prob <- rbind(all.prob, to.append)
  #   }
  # }
  
  
  p <- ggplot(all.sim, aes(x=as.factor(kmer), y=prop.fp)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.position="none") + ylab("Proportion matching k-mers") + xlab("Kmer size")
  p
  ggsave(paste("kmer_random_match_readlen", read.length, ".png", sep=""))
  
  
  # subsample <- subset(all.sim, kmer > 1)
  # p <- ggplot(subsample, aes(x=as.factor(kmer), y=prop.fp)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.position="none") + ylab("Proportion matching k-mers") + xlab("Kmer size")
  # p
  # ggsave(paste("kmer_random_match_readlen", read.length, "_lessk1.png", sep=""))
  
  
  # calculate effects of cutoffs
  cut.offs <- seq(0, 100, by=2.5) / 100
  all.cutoffs <- data.frame(kmer = c(), mu = c(), cut.off = c(), prop.matches = c(), prop.fp = c(), tp = c(), fp = c(), prop.matches.mash = c(), prop.fp.mash = c(), tp.mash = c(), fp.mash = c())
  
  for (kmer.length in kmer.sizes)
  {
    for (mur in mu.rates) 
    {
      subsample <- subset(all.sim, kmer == kmer.length & mu == mur)
      subsample$true.rate <- subsample$prop.kmers - subsample$prop.fp
      subsample$true.rate.mash <- subsample$mash - subsample$mash.fp
      for (co in cut.offs)
      {
        above.cut.off.tp <- sum(subsample$prop.kmers >= co) / length(subsample$prop.kmers)
        above.cut.off.fp <- sum(subsample$prop.fp >= co) / length(subsample$prop.kmers)
        # true positives are those where false positives don't affect decision i.e. prop.kmers - false positives > cutoff
        true.positives <- subsample[,"true.rate"] >= co & subsample[,"prop.kmers"] >= co
        false.positives <- subsample[,"true.rate"] < co & subsample[,"prop.kmers"] >= co
        
        above.cut.off.mash <- sum(subsample$mash >= co) / length(subsample$prop.kmers)
        above.cut.off.mash.fp <- sum(subsample$mash.fp >= co) / length(subsample$prop.kmers)
        # true positives are those where false positives don't affect decision i.e. prop.kmers - false positives > cutoff
        true.positives.mash <- subsample[,"true.rate.mash"] >= co & subsample[,"mash"] >= co
        false.positives.mash <- subsample[,"true.rate.mash"] < co & subsample[,"mash"] >= co
        
        to.append <- data.frame(kmer = kmer.length, mu = mur, cut.off = co, prop.matches = above.cut.off, prop.fp = above.cut.off.fp, tp = sum(true.positives) / length(subsample$prop.kmers), fp=sum(false.positives) / length(subsample$prop.kmers), prop.matches.mash = above.cut.off.mash, prop.fp.mash = above.cut.off.mash.fp, tp.mash = sum(true.positives.mash) / length(subsample$prop.kmers), fp.mash=sum(false.positives.mash) / length(subsample$prop.kmers))
        all.cutoffs <- rbind(all.cutoffs, to.append)
      }
    }
  }
  
  
  # generate ROC
  #all.cutoffs <- all.cutoffs[order(all.cutoffs$cut.off, decreasing=TRUE),]
  
  mu.plot <- 0.02
  subsample <- subset(all.cutoffs, mu == mu.plot)
  p <- ggplot(subsample, aes(x=as.factor(cut.off), y=prop.matches, colour = as.factor(kmer), group = as.factor(kmer))) + facet_grid(~mu) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("Proportion reads assigned target") + xlab("Cut-off")
  p
  ggsave(paste("kmer_cutoff_mu", mu.plot, "_readlen", read.length, ".png", sep=""))
  # p <- ggplot(subsample, aes(x=fp, y=tp, colour = as.factor(kmer), group = as.factor(kmer))) + facet_grid(~mu) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("True positives") + xlab("False positives")
  # p
  # ggsave(paste("kmer_ROC_mu", mu.plot, "_readlen", read.length, ".png", sep=""))
  
  # mu.plot <- 0.05
  # subsample <- subset(all.cutoffs, mu == mu.plot)
  # p <- ggplot(subsample, aes(x=as.factor(cut.off), y=prop.matches, colour = as.factor(kmer), group = as.factor(kmer))) + facet_grid(~mu) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("Proportion reads assigned target") + xlab("Cut-off")
  # p
  # ggsave(paste("kmer_cutoff_mu", mu.plot, "_readlen", read.length, ".png", sep=""))
  # # p <- ggplot(subsample, aes(x=fp, y=tp, colour = as.factor(kmer), group = as.factor(kmer))) + facet_grid(~mu) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("True positives") + xlab("False positives")
  # # p
  # # ggsave(paste("kmer_ROC_mu", mu.plot, "_readlen", read.length, ".png", sep=""))
  
  mu.plot <- 0.1
  subsample <- subset(all.cutoffs, mu == mu.plot)
  p <- ggplot(subsample, aes(x=as.factor(cut.off), y=prop.matches, colour = as.factor(kmer), group = as.factor(kmer))) + facet_grid(~mu) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("Proportion reads assigned target") + xlab("Cut-off")
  p
  ggsave(paste("kmer_cutoff_mu", mu.plot, "_Q10_readlen", read.length, ".png", sep=""))
  # p <- ggplot(subsample, aes(x=fp, y=tp, colour = as.factor(kmer), group = as.factor(kmer))) + facet_grid(~mu) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("True positives") + xlab("False positives")
  # p
  # ggsave(paste("kmer_ROC_mu", mu.plot, "_readlen", read.length, ".png", sep=""))
  
  # mu.plot <- 0.12
  # subsample <- subset(all.cutoffs, mu == mu.plot)
  # p <- ggplot(subsample, aes(x=as.factor(cut.off), y=prop.matches, colour = as.factor(kmer), group = as.factor(kmer))) + facet_grid(~mu) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("Proportion reads assigned target") + xlab("Cut-off")
  # p
  # ggsave(paste("kmer_cutoff_mu", mu.plot, "_Q9_readlen", read.length, ".png", sep=""))
  # # p <- ggplot(subsample, aes(x=fp, y=tp, colour = as.factor(kmer), group = as.factor(kmer))) + facet_grid(~mu) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("True positives") + xlab("False positives")
  # # p
  # # ggsave(paste("kmer_ROC_mu", mu.plot, "_readlen", read.length, ".png", sep=""))
  # 
  # mu.plot <- 0.13
  # subsample <- subset(all.cutoffs, mu == mu.plot)
  # p <- ggplot(subsample, aes(x=as.factor(cut.off), y=prop.matches, colour = as.factor(kmer), group = as.factor(kmer))) + facet_grid(~mu) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("Proportion reads assigned target") + xlab("Cut-off")
  # p
  # ggsave(paste("kmer_cutoff_mu", mu.plot, "_Q9_readlen", read.length, ".png", sep=""))
  # # p <- ggplot(subsample, aes(x=fp, y=tp, colour = as.factor(kmer), group = as.factor(kmer))) + facet_grid(~mu) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("True positives") + xlab("False positives")
  # # p
  # # ggsave(paste("kmer_ROC_mu", mu.plot, "_readlen", read.length, ".png", sep=""))
  # 
  # mu.plot <- 0.15
  # subsample <- subset(all.cutoffs, mu == mu.plot)
  # p <- ggplot(subsample, aes(x=as.factor(cut.off), y=prop.matches, colour = as.factor(kmer), group = as.factor(kmer))) + facet_grid(~mu) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("Proportion reads assigned target") + xlab("Cut-off")
  # p
  # ggsave(paste("kmer_cutoff_mu", mu.plot, "_readlen", read.length, ".png", sep=""))
  # # p <- ggplot(subsample, aes(x=fp, y=tp, colour = as.factor(kmer), group = as.factor(kmer))) + facet_grid(~mu) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("True positives") + xlab("False positives")
  # # p
  # # ggsave(paste("kmer_ROC_mu", mu.plot, "_readlen", read.length, ".png", sep=""))
  # 
  # mu.plot <- 0.16
  # subsample <- subset(all.cutoffs, mu == mu.plot)
  # p <- ggplot(subsample, aes(x=as.factor(cut.off), y=prop.matches, colour = as.factor(kmer), group = as.factor(kmer))) + facet_grid(~mu) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("Proportion reads assigned target") + xlab("Cut-off")
  # p
  # ggsave(paste("kmer_cutoff_mu", mu.plot, "_Q8_readlen", read.length, ".png", sep=""))
  # # p <- ggplot(subsample, aes(x=fp, y=tp, colour = as.factor(kmer), group = as.factor(kmer))) + facet_grid(~mu) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("True positives") + xlab("False positives")
  # # p
  # # ggsave(paste("kmer_ROC_mu", mu.plot, "_readlen", read.length, ".png", sep=""))
  
  mu.plot <- 0.2
  subsample <- subset(all.cutoffs, mu == mu.plot)
  p <- ggplot(subsample, aes(x=as.factor(cut.off), y=prop.matches, colour = as.factor(kmer), group = as.factor(kmer))) + facet_grid(~mu) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("Proportion reads assigned target") + xlab("Cut-off")
  p
  ggsave(paste("kmer_cutoff_mu", mu.plot, "_readlen", read.length, ".png", sep=""))
  # p <- ggplot(subsample, aes(x=fp, y=tp, colour = as.factor(kmer), group = as.factor(kmer))) + facet_grid(~mu) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("True positives") + xlab("False positives")
  # p
  # ggsave(paste("kmer_ROC_mu", mu.plot, "_readlen", read.length, ".png", sep=""))
  
  # kmer.sample <- 11
  # subsample <- subset(all.cutoffs, kmer == kmer.sample)
  # p <- ggplot(subsample, aes(x=as.factor(cut.off), y=prop.matches, colour = as.factor(mu), group = as.factor(mu))) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("Proportion reads assigned target") + xlab("Cut-off")
  # p
  # ggsave(paste("kmer_cutoff_k", kmer.sample, "_readlen", read.length, ".png", sep=""))
  # 
  kmer.sample <- 13
  subsample <- subset(all.cutoffs, kmer == kmer.sample)
  p <- ggplot(subsample, aes(x=as.factor(cut.off), y=prop.matches, colour = as.factor(mu), group = as.factor(mu))) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("Proportion reads assigned target") + xlab("Cut-off")
  p
  ggsave(paste("kmer_cutoff_k", kmer.sample, "_readlen", read.length, ".png", sep=""))
  
#   kmer.sample <- 9
#   subsample <- subset(all.cutoffs, kmer == kmer.sample)
#   p <- ggplot(subsample, aes(x=as.factor(cut.off), y=prop.matches, colour = as.factor(mu), group = as.factor(mu))) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("Proportion reads assigned target") + xlab("Cut-off")
#   p
#   ggsave(paste("kmer_cutoff_k", kmer.sample, "_readlen", read.length, ".png", sep=""))
#   
   write.csv(all.sim, paste("all_simulations_readlen", read.length, ".csv", sep = ""), row.names=FALSE)
   write.csv(all.cutoffs, paste("all_cutoffs_readlen", read.length, ".csv", sep = ""), row.names=FALSE)
}

mu.plot <- 0.16
read.length.choice <- 200
read.length <- read.lengths[[1]]
all.sim <- read.csv(file = paste('kmer_simulations/new/all_simulations.csv', sep = ""), stringsAsFactors = FALSE)
all.cutoffs <- read.csv(file = paste('kmer_simulations/new/all_cutoffs.csv', sep=""), stringsAsFactors = FALSE)

# for (ri in 1:length(read.lengths))
# {
#   read.length <- read.lengths[[ri]]
#   temp.sim <- read.csv(file = paste('kmer_simulations/new/all_simulations_readlen', read.length, '.csv', sep = ""), stringsAsFactors = FALSE)
#   temp.sim$read.length <- read.length
# 
#   all.sim <- rbind(all.sim, temp.sim)
# }
# 
# write.csv(all.sim, paste("all_simulations.csv", sep = ""), row.names=FALSE)

# # recalculate cutoffs
# cut.offs <- seq(0, 100, by=2.5) / 100
# all.cutoffs <- data.frame(kmer = c(), mu = c(), cut.off = c(), prop.matches = c(), prop.fp = c(), tp = c(), fp = c(), prop.matches.mash = c(), prop.fp.mash = c(), tp.mash = c(), fp.mash = c(), read.length = c())
# 
# for (read.l in read.lengths)
#   {
#   for (kmer.length in kmer.sizes)
#   {
#     for (mur in mu.rates)
#     {
#       subsample <- subset(all.sim, kmer == kmer.length & mu == mur & read.length == read.l)
#       subsample$true.rate <- subsample$prop.kmers - subsample$prop.fp
#       subsample$true.rate.mash <- subsample$mash - subsample$mash.fp
#       for (co in cut.offs)
#       {
#         above.cut.off.tp <- sum(subsample$prop.kmers >= co) / length(subsample$prop.kmers)
#         above.cut.off.fp <- sum(subsample$prop.fp >= co) / length(subsample$prop.kmers)
#         # true positives are those where false positives don't affect decision i.e. prop.kmers - false positives > cutoff
#         true.positives <- subsample[,"true.rate"] >= co & subsample[,"prop.kmers"] >= co
#         false.positives <- subsample[,"true.rate"] < co & subsample[,"prop.kmers"] >= co
# 
#         above.cut.off.mash <- sum(subsample$mash >= co) / length(subsample$prop.kmers)
#         above.cut.off.mash.fp <- sum(subsample$mash.fp >= co) / length(subsample$prop.kmers)
#         # true positives are those where false positives don't affect decision i.e. prop.kmers - false positives > cutoff
#         true.positives.mash <- subsample[,"true.rate.mash"] >= co & subsample[,"mash"] >= co
#         false.positives.mash <- subsample[,"true.rate.mash"] < co & subsample[,"mash"] >= co
# 
#         to.append <- data.frame(kmer = kmer.length, mu = mur, cut.off = co, prop.matches = above.cut.off.tp, prop.fp = above.cut.off.fp, tp = sum(true.positives) / length(subsample$prop.kmers), fp=sum(false.positives) / length(subsample$prop.kmers), prop.matches.mash = above.cut.off.mash, prop.fp.mash = above.cut.off.mash.fp, tp.mash = sum(true.positives.mash) / length(subsample$prop.kmers), fp.mash=sum(false.positives.mash) / length(subsample$prop.kmers), read.length = read.l)
#         all.cutoffs <- rbind(all.cutoffs, to.append)
#       }
#     }
#   }
# }
# 
# write.csv(all.cutoffs, paste("all_cutoffs.csv", sep = ""), row.names=FALSE)

subsample <- subset(all.sim, (kmer == 7 | kmer == 11 | kmer == 15 | kmer == 19 | kmer == 23 | kmer == 27 | kmer == 31) & read.length == read.length.choice)
p <- ggplot(subsample, aes(x=as.factor(mu), y=prop.kmers)) + facet_grid(~ kmer) + geom_boxplot() + scale_x_discrete(breaks = function(x){x[c(TRUE, FALSE)]}) + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.position="none") + ylab("Proportion matching k-mers") + xlab("Mutation rate") + scale_y_continuous(limits=c(NA, 1))# + stat_summary(fun.y=median, colour="black", geom="text", show_guide = FALSE, vjust=-1, aes( label=round(after_stat(y), digits=2)))
p
ggsave(paste("kmer_mu_comp_readlen", read.length.choice, ".svg", sep=""))
p <- ggplot(subsample, aes(x=as.factor(kmer), y=prop.fp)) + geom_boxplot() + theme_light() + theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.position="none") + ylab("Prop. random matching k-mers") + xlab("K-mer size") + scale_y_continuous(breaks = seq(0,1,0.1))# + stat_summary(fun.y=median, colour="black", geom="text", show_guide = FALSE, vjust=-1, aes( label=round(after_stat(y), digits=2)))
p
ggsave(paste("kmer_fp_readlen", read.length.choice, ".svg", sep=""))

# read in previous analysis
subsample <- subset(all.cutoffs, (mu == 0 | mu == 0.01 | mu == 0.1 | mu == 0.12 | mu == 0.16) & read.length == read.length.choice & (kmer == 7 | kmer == 11 | kmer == 15 | kmer == 19 | kmer == 23 | kmer == 27 | kmer == 31))
p <- ggplot(subsample, aes(x=cut.off, y=prop.matches, colour = as.factor(kmer), group = as.factor(kmer))) + facet_grid(~mu) + geom_point() + geom_line() + theme_light() + theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 13), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("Prop. reads assigned target") + xlab("Cut-off") + labs(color = "K-mer size") + scale_x_continuous(breaks = seq(0,1,0.2))# + geom_vline(xintercept = 0.5, color="black", linetype = "longdash")
p
ggsave(paste("cutoff_v_mu", read.length.choice, ".svg", sep=""))

# draw ROC curve
subsample <- subset(all.cutoffs, (mu == 0 | mu == 0.01 | mu == 0.1 | mu == 0.12 | mu == 0.16) & read.length == read.length.choice & (kmer == 7 | kmer == 11 | kmer == 15 | kmer == 19 | kmer == 23 | kmer == 27 | kmer == 31))
# correct error at 0% cutoff
subsample$prop.matches[subsample$cut.off == 0] <- 1.0
subsample$prop.fp[subsample$cut.off == 0] <- 1.0
subsample$prop.matches[subsample$cut.off == 1] <- 0.0
subsample$prop.fp[subsample$cut.off == 1] <- 0.0

p <- ggplot(subsample[order(subsample$cut.off),], aes(x=prop.fp, y=prop.matches, colour = as.factor(kmer), group = as.factor(kmer))) + facet_grid(~mu) + geom_path(linewidth=1.2) + theme_light() + theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 13), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("True positive rate") + xlab("False positive rate") + labs(color = "K-mer size") + scale_x_continuous(breaks = seq(0,1,0.2)) + geom_abline(intercept = 0, slope = 1, color="black", linetype = "longdash")
p
ggsave(paste("ROC_len", read.length.choice, ".svg", sep=""))


subsample <- subset(all.sim, read.length == read.length.choice)
p <- ggplot(subsample, aes(x=as.factor(kmer), y=prop.fp)) + geom_boxplot() + theme_light() + theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.position="none") + ylab("Prop. random matching") + xlab("K-mer size") + scale_y_continuous(breaks = seq(0,1,0.1))
p
ggsave(paste("kmer_fp_readlen", read.length.choice, ".svg", sep=""))


# calculate RMSE of y=x
RMSE <- all.sim[all.sim$mu <= 0.2,] %>%
  group_by(kmer, read.length) %>%
  #sqrt(mean(((1 - mu) - mash)^2))
  summarise(RSME = sqrt(mean(((1 - mu) - mash)^2)))

subsample <- subset(all.cutoffs, (mu == 0 | mu == 0.01 | mu == 0.1 | mu == 0.12 | mu == 0.16) & read.length == read.length.choice & (kmer == 7 | kmer == 11 | kmer == 15 | kmer == 19 | kmer == 23 | kmer == 27 | kmer == 31))
p <- ggplot(subsample, aes(x=cut.off, y=prop.matches.mash, colour = as.factor(kmer), group = as.factor(kmer))) + facet_grid(~mu) + geom_point() + geom_line() + theme_light() + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("Prop. reads assigned target") + xlab("Cut-off") + labs(color = "Kmer size") + scale_x_continuous(breaks = seq(0,1,0.2))
p
ggsave(paste("mash_mu_comp_readlen", read.length.choice, ".svg", sep=""))

subsample <- subset(all.sim, read.length == read.length.choice & (kmer == 7 | kmer == 11 | kmer == 15 | kmer == 19 | kmer == 23 | kmer == 27 | kmer == 31))
p <- ggplot(subsample, aes(x=as.factor(kmer), y=mash.fp)) + geom_boxplot() + theme_light() + theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.position="none") + ylab("Mash-like index (S)") + xlab("K-mer size") + scale_y_continuous(breaks = seq(0,1,0.1))
p
ggsave(paste("mash_fp_readlen", read.length.choice, ".svg", sep=""))

subsample <- subset(all.sim, read.length == read.length.choice & (kmer == 7 | kmer == 11 | kmer == 15 | kmer == 19 | kmer == 23 | kmer == 27 | kmer == 31))
p <- ggplot(subsample, aes(x=1-mu, y=mash, colour = as.factor(kmer))) + facet_grid(~kmer, scales="free_y") + geom_point(show.legend = FALSE) + geom_smooth(method=lm, level=0.95, colour = "black") + geom_abline(slope = 1, intercept=0, linetype="dashed", color = "black") + theme_light() + theme(axis.text.x = element_text(size = 10, angle=45, vjust=1, hjust=1), axis.text.y = element_text(size = 10), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("Mash-like index (S)") + xlab("Real mutation rate") + labs(color = "Kmer size") + scale_x_continuous(breaks = seq(0,1,0.04)) + scale_y_continuous(limits = c(0, NA), breaks = seq(0,1,0.2)) # + stat_summary(fun.y=median, colour="black", geom="point", shape=18, size=2,show_guide = FALSE)
p
ggsave(file="mash_mu_concordance.svg", plot=p)

RMSE.subsample <- subset(RMSE, (kmer == 7 | kmer == 11 | kmer == 15 | kmer == 19 | kmer == 23 | kmer == 27 | kmer == 31))
p <- ggplot(RMSE.subsample, aes(x=as.factor(kmer), y=RSME, fill = as.factor(read.length))) + geom_col(position = "dodge") + theme_light() + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("RMSE") + xlab("K-mer size") + labs(fill = "Read length") + scale_y_continuous(breaks = seq(0,1,0.1))
p
ggsave(file="RMSE_mash.svg", plot=p)

# correct error at 0% cutoff
subsample <- subset(all.cutoffs, (mu == 0 | mu == 0.01 | mu == 0.1 | mu == 0.12 | mu == 0.16) & read.length == read.length.choice & (kmer == 7 | kmer == 11 | kmer == 15 | kmer == 19 | kmer == 23 | kmer == 27 | kmer == 31))
subsample$prop.matches.mash[subsample$cut.off == 0] <- 1.0
subsample$prop.fp.mash[subsample$cut.off == 0] <- 1.0
subsample$prop.matches.mash[subsample$cut.off == 1] <- 0.0
subsample$prop.fp.mash[subsample$cut.off == 1] <- 0.0
p <- ggplot(subsample[order(subsample$cut.off),], aes(x=prop.fp.mash, y=prop.matches.mash, colour = as.factor(kmer), group = as.factor(kmer))) + facet_grid(~mu) + geom_path(linewidth=1.2) + theme_light() + theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 13), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12)) + ylab("True positive rate") + xlab("False positive rate") + labs(color = "K-mer size") + scale_x_continuous(breaks = seq(0,1,0.2)) + geom_abline(intercept = 0, slope = 1, color="black", linetype = "longdash")
p
ggsave(paste("ROC_mash_len", read.length.choice, ".svg", sep=""))

