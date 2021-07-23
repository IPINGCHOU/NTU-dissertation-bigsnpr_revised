library(parallel)
library(bigsnpr)
library(doParallel)
library(ggplot2)


## Without this code, the MPI session may fail to shut down properly in case an
## error occurs, and then the script won't terminate until it hits the walltime
# options(error=quote(assign(".mpi.err", FALSE, env = .GlobalEnv)))

snp_fastImpute_GPU <- function(Gna, infos.chr,
                               alpha = 1e-4,
                               size = 200,
                               p.train = 0.8,
                               n.cor = nrow(Gna),
                               seed = NA,
                               ncores = 1,
                               cal_mode = 'gpu_hist') {

      imputeChr_GPU <- function(X, X2, infos.imp, ind.chr, alpha, size,
                            p.train, n.cor, seed, ncores, cal_mode) {

        old <- .Random.seed
        on.exit(.Random.seed <<- old, add = TRUE)

    # Do something only if there is still something to do
        if (any(is.na(infos.imp[1, ind.chr]))) {
          n <- nrow(X)

          if (!is.na(seed)) set.seed(seed)

          # correlation between variants
          print(paste('Corr computing...'))
          corr_start <- Sys.time()
          corr <- snp_cor(
            Gna       = X,
            ind.row   = sort(sample(n, size = n.cor)),
            ind.col   = ind.chr,
            size      = size,
            alpha     = alpha,
            fill.diag = FALSE,
            ncores    = ncores
          )
          print(paste('Corr computed...'))
          print(Sys.time() - corr_start)

          # imputation ===
          print(paste('Start imputation...'))
          imputation_start <- Sys.time()
          # for (i in seq_along(ind.chr)) {
          cl <- makeCluster(ncores)
          print(cl)
          registerDoParallel(cl)
          total_length <- as.numeric(length(seq_along(ind.chr)))

          foreach(i = 1:total_length, .packages='bigsnpr') %dopar%{

            if (!is.na(seed)) set.seed(seed + i)

            snp <- ind.chr[i]
            # Do something only if it wasn't done before
            if (is.na(infos.imp[1, snp])) {

              X.label <- X[, snp]

              nbna <- length(indNA <- which(is.na(X.label)))
              if (nbna > 0) {
                indNoNA <- setdiff(seq_len(n), indNA)
                ind.train <- sort(sample(indNoNA, size = p.train * length(indNoNA)))
                ind.val <- setdiff(indNoNA, ind.train)

                ind.col <- ind.chr[which(corr[, i] != 0)]
                if (length(ind.col) < 5L)
                  ind.col <- intersect(setdiff(-size:size + snp, snp), ind.chr)

                data.train <- xgboost::xgb.DMatrix(
                  label = X.label[ind.train],
                  data  = X2[ind.train, ind.col, drop = FALSE])

                bst.params <- list(
                  objective  = "binary:logistic",
                  max_depth  = 4,
                  base_score = min(max(1e-7, mean(X.label[ind.train])), 1 - 1e-7),
                  verbose    = 0,
                  nthread    = 1,
                  tree_method = cal_mode ###neo
                )

                bst <- xgboost::xgb.train(
                  data    = data.train,
                  params  = bst.params,
                  nrounds = 10
                )

                # error of validation
                pred2 <- stats::predict(bst, X2[ind.val, ind.col, drop = FALSE])
                infos.imp[2, snp] <- mean(round(2 * pred2) != (2 * X.label[ind.val]))
                # imputation
                pred <- stats::predict(bst, X2[indNA, ind.col, drop = FALSE])
                X2[indNA, snp] <- as.raw(round(2 * pred) + 4)

              }

              # this variant is done
              infos.imp[1, snp] <- nbna / n
            }
          }
          stopCluster(cl)
          print(paste('End imputation...'))
          print(Sys.time() - imputation_start)
        }

    invisible()
    }

    if (!requireNamespace("xgboost", quietly = TRUE))
    stop2("Please install package 'xgboost'.")

    print(paste('Preprocessing...'))
    bigsnpr:::check_args(infos.chr = "assert_lengths(infos.chr, cols_along(Gna))")

    CODE_IMPUTE_LABEL <- c(0, 0.5, 1, rep(NA, 253))
    CODE_IMPUTE_PRED  <- c(0, 1, 2, NA, 0, 1, 2, rep(NA, 249))

    X  <- Gna$copy(code = CODE_IMPUTE_LABEL)
    X2 <- Gna$copy(code = CODE_IMPUTE_PRED)

    infos.imp <- bigsnpr:::FBM_infos(Gna)
    
    # split CHR to 23 chromosomes
    # +1 for start the idx from 1
    idx <- seq_along(CHR)
    ind.chrs <- split(idx, cut_number(idx, 22))
    total_len <- as.numeric(length(ind.chrs))

    for (i in 1:total_len) {
        imputeChr_GPU(X, X2, infos.imp, ind.chrs[[i]], alpha, size, p.train, n.cor, seed, ncores, cal_mode)
    }

    infos.imp

}

#========================
# NCORES <- 16
# DEVICE <- 1
# file_path <- "./fake_snps/2k_2k_0.05.bed"
# file_path <- './med_data/sample4.bed'

args = commandArgs(trailingOnly=TRUE)
NCORES <- as.numeric(args[1])
DEVICE <- as.numeric(args[2])
file_path <- as.character(args[3])
temp_path <- as.character((args[4]))
 
set.seed(1024)
SEED <- 1024

# reading files
# if temp_path == 'tmp , take tempfile (default /tmp folder) as temp folder
if(temp_path == 'tmp'){temp_path = tempfile()}

rds <- paste0(temp_path, '.rds')

if(file.exists(rds)==FALSE){
  rds <- snp_readBed(file_path, backingfile = temp_path)
}
print(paste('rds path: ', rds))

obj<- snp_attach(rds)
G <- obj$genotypes
CHR <- obj$map$chromosome # ind.chr

NROWS <- dim(G)[1]
NCOLS <- dim(G)[2]
CHR_LEN <- length(CHR)
print(paste('CHR length: ', CHR_LEN))
print(paste('Matrix size: ', NROWS, ' x ', NCOLS))
print(paste('Core: ', NCORES))

default_nproc_blas <- function(){
   cl <- parallel::makePSOCKcluster(1)
   on.exit(parallel::stopCluster(cl), add = TRUE)
   parallel::clusterEvalQ(cl, RhpcBLASctl::blas_get_num_procs())[[1]]

}

if(default_nproc_blas() > 1){
   options(default.nproc.blas = 1)
}

if(DEVICE==1){
  a = Sys.time()
  print(paste('Now using GPU version snpfastImpute'))
  G2 <- snp_fastImpute_GPU(G, CHR, ncores = NCORES, seed = SEED, cal_mode = 'gpu_hist')
  print(Sys.time()-a)
}else if(DEVICE==2){
  print(paste('Now using CPU version snpfastImpute'))
  G2 <- snp_fastImpute_GPU(G, CHR, ncores = NCORES, seed = SEED, cal_mode = 'hist')
}else if(DEVICE==3){
  print(paste('Now using original version snpfastImpute with CPU'))
  G2 <- snp_fastImpute(G, CHR, ncores = NCORES, seed = SEED)
}


print(paste('done'))

#=== get imputed
print(paste('check imputed infos.imp and G'))
print(G2[1:2, 1:10])
big_counts(G, ind.col = 1:10)
# You need to change the code of G
# To make this permanent, you need to save (modify) the file on disk
obj$genotypes$code256 <- CODE_IMPUTE_PRED
# obj<- snp_save(obj)
big_counts(obj$genotypes, ind.col = 1:10)

print(paste('checked'))

