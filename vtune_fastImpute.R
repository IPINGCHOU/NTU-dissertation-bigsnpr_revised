# For vtune setting
options(digits=10)
args = commandArgs(trailingOnly=TRUE)

if (length(args)<3) {
  stop("Args not enough (row, column, cores)", call.=FALSE)
}

# init testing parameters
library(bigsnpr)
NROWS <- as.numeric(args[1])
NCOLS <- as.numeric(args[2])
NCORES <- as.numeric(args[3])
prob_of_NAs <- as.double(args[4])
code <- rep(NA_real_, 256)
code[1:3] <- c(0,1,2)

CODE_IMPUTE_LABEL <- c(0, 0.5, 1, rep(NA, 253))
CODE_IMPUTE_PRED  <- c(0, 1, 2, NA, 0, 1, 2, rep(NA, 249))
# Matrix size
code <- rep(NA_real_, 256)
code[1:3] <- c(0,1,2)
# Generating testing matrix
G_toy <- FBM.code256(NROWS, NCOLS, code, init = sample(as.raw(0:2), NROWS*NCOLS, replace = TRUE)) 

res_of_probs <- (1-prob_of_NAs) / 3
G_toy[] <- sample(c(0,1,2,3), size = length(G_toy), replace = TRUE, prob = c(res_of_probs,res_of_probs,res_of_probs,prob_of_NAs))
G_toy_CHR <- rep(1,NCOLS)
G_toy_X <- G_toy$copy(code = CODE_IMPUTE_LABEL)

print(paste('Matrix size: ', NROWS, ' x ', NCOLS))
print(paste('The total NAs: ', length(G_toy)*prob_of_NAs))
print(paste('Core: ', NCORES))
snp_fastImpute(G_toy_X, G_toy_CHR, ncores = NCORES)
print(paste('done'))