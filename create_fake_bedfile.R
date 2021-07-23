rm(list=ls())
library(bigsnpr)
library(bit64)

NROWS <- 14000
NCOLS <- 64000
prob_of_NAs <- 0.05
res_of_probs = (1-prob_of_NAs) / 3

fake <- snp_fake(NROWS, NCOLS)
G <- fake$genotypes
#G[] <- sample(c(0,1,2,3), size = length(G), replace = TRUE, prob = c(res_of_probs,res_of_probs,res_of_probs,prob_of_NAs))

#for(i in 1:NROWS){
 # G[i,] <- sample(c(0,1,2,3), size = NCOLS, replace = TRUE, prob = c(res_of_probs,res_of_probs,res_of_probs,prob_of_NAs))
#}

library(foreach)
library(doParallel)
cl = makeCluster(8)
registerDoParallel(cl)
print('parallel generating fake_snps')
foreach::foreach(i=1:NROWS, .combine=c, .packages = 'bigsnpr') %dopar%{
  G[i,] <- sample(c(0,1,2,3), size = NCOLS, replace = TRUE, prob = c(res_of_probs,res_of_probs,res_of_probs,prob_of_NAs))
}
stopCluster(cl)
print('generated...')
print('saving bed file...')

#G[1:NROWS,] <- apply(1:NROWS,sample(c(0,1,2,3), size = NCOLS, replace = TRUE, prob = c(res_of_probs,res_of_probs,res_of_probs,prob_of_NAs)))

tmp <- tempfile(fileext = ".bed")

save_path = '/home/bigsnpr_medicine/fake_snps/100_100_0.05.bed'
bed <- snp_writeBed(fake, save_path)
print(paste('done'))
