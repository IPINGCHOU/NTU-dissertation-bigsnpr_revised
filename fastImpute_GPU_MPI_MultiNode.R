library(bigsnpr)
library(doParallel)
library(Rmpi)
library(doMPI)
library(ggplot2)

## Without this code, the MPI session may fail to shut down properly in case an
## error occurs, and then the script won't terminate until it hits the walltime
# options(error=quote(assign(".mpi.err", FALSE, env = .GlobalEnv)))

FBM_infos_cur <- function(Gna) {

    base <- sub_bk(Gna$backingfile, "-infos-impute", stop_if_not_ext = FALSE)
    rds <- paste0(base, ".rds")

    if (file.exists(rds)) {
        cat("Attach\n")
        readRDS(rds)
    } else {
        cat("Create\n")
        FBM(2, ncol(Gna), backingfile = base, init = NA_real_)$save()
    }
}


snp_fastImpute_GPU <- function(Gna, infos.chr,
                               alpha = 1e-04,
                               size = 200,
                               p.train = 0.8,
                               n.cor = nrow(Gna),
                               seed = NA,
                               ncores = 1,
                               cal_mode = 'gpu_hist') {
    comb <- function(...){
        chr <- c()
        info <- c()
        for(i in list(...)){
            chr <- cbind(chr, i$impute_chromesome)
            info <- cbind(info, i$impute_info)
        }
        return(list(impute_chromesome=chr, impute_info=info))
    }

    
    FBM_infos_cur <- function(Gna) {

        base <- sub_bk(Gna$backingfile, "-infos-impute", stop_if_not_ext = FALSE)
        rds <- paste0(base, ".rds")

        if (file.exists(rds)) {
            cat("Attach\n")
            readRDS(rds)
        } else {
            cat("Create\n")
            FBM(2, ncol(Gna), backingfile = base, init = NA_real_)$save()
        }
    }
    

    imputeChr_GPU <- function(X, X2, infos.imp, ind.chr, alpha, size,
                            p.train, n.cor, seed, ncores, cal_mode){
    
        # get current slave info
        id <- mpi.comm.rank(comm = 0)
        np <- mpi.comm.size(comm = 0)
        hostname <- mpi.get.processor.name(short = FALSE)

        msg <- sprintf("Hello world from process %03d of %03d, on host %s\n", id, np, hostname)
        cat(msg)
        cat("\n")
        
        # === imputeChr_GPU start ===
        old <- .Random.seed
        on.exit(.Random.seed <<- old, add = TRUE)
        
        # Do somthing only if there is still something to do
        if(any(is.na(infos.imp[1, ind.chr]))){
            
            # prepare impute matrix
            impute_chromesome <- X2[, ind.chr]
            impute_info <- infos.imp[, ind.chr]
            offset <- ind.chr[1] -1
            
            n <- nrow(X)
            
            if(!is.na(seed)) set.seed(seed)
            
            # get correlation between variants
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
            
#             imputation ===
            print(paste('Start imputation...'))
            imputation_start <- Sys.time()
            cluster_spec <- rep(hostname, ncores) 
            incl <- makeCluster(cluster_spec)
            print('=======')
            print(incl)
            print('=======')
            registerDoParallel(incl)
            total_length <- as.numeric(length(seq_along(ind.chr)))

            impute_inner_chr_info <- foreach(i = 1:total_length, .packages = 'bigsnpr', .multicombine=TRUE, .combine=comb) %dopar% {
#             for(i in seq_along(ind.chr)){
                if(!is.na(seed)) set.seed(seed + i)

                snp <- ind.chr[i]
                # Do something only if it wasn't done before
                if(is.na(infos.imp[1, snp])){
                    X.label <- X[, snp]
                    nbna <- length(indNA <- which(is.na(X.label)))

                    if(nbna > 0){
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
                        #infos.imp[2, snp] <- mean(round(2 * pred2) != (2 * X.label[ind.val]))
                        impute_info[2, snp - offset] <- mean(round(2 * pred2) != (2 * X.label[ind.val]))
                        # imputation
                        pred <- stats::predict(bst, X2[indNA, ind.col, drop = FALSE])
                        #X2[indNA, snp] <- as.raw(round(2 * pred) + 4)
                        impute_chromesome[indNA, snp - offset] <- as.double(round(2 * pred) + 4)
                    }
                    # this variant is done
#                     infos.imp[1, snp] <- nbna / n
                    impute_info[1, snp - offset] <- nbna / n
                }
            return(list("impute_chromesome" = impute_chromesome[, snp-offset], "impute_info" = impute_info[, snp-offset]))    
            }
            stopCluster(incl)
            print(paste('End imputation...'))
            print(Sys.time() - imputation_start)
        }
    #invisible()
    print(dim(impute_inner_chr_info$impute_chromesome))
    print(dim(impute_inner_chr_info$impute_info))
    print(impute_inner_chr_info$impute_chromesome[1:10, 1:10])
    print(impute_inner_chr_info$impute_info[,1:10])
    return(impute_inner_chr_info)
    }
    
    if (!requireNamespace("xgboost", quietly = TRUE))
    stop2("Please install package 'xgboost'.")

    print(paste('Preprocessing...'))
    bigsnpr:::check_args(infos.chr = "assert_lengths(infos.chr, cols_along(Gna))")
    
    CODE_IMPUTE_LABEL <- c(0, 0.5, 1, rep(NA, 253))
    CODE_IMPUTE_PRED  <- c(0, 1, 2, NA, 0, 1, 2, rep(NA, 249))

    X  <- Gna$copy(code = CODE_IMPUTE_LABEL)
    X2 <- Gna$copy(code = CODE_IMPUTE_PRED)

    infos.imp <- FBM_infos_cur(Gna)

    #   ind.chrs <- split(seq_along(infos.chr), infos.chr)
    # split CHR to 23 chromosomes
    # +1 for start the idx from 1
    idx <- seq_along(CHR)
    ind.chrs <- split(idx, cut_number(idx, 22))
    total_len <- as.numeric(length(ind.chrs))
    
    print('able slaves')
    able_slaves <- mpi.universe.size() - 1
    print(able_slaves)
    print('=========')
    
    MPIcl <- startMPIcluster(maxcores = 4, verbose = TRUE, logdir = "/home/u8294235/bigsnpr_medicine/snp_fastImpute_modify/MPI_slave_log/")
    print('=============')
    print(MPIcl)
    print('=============')
    registerDoMPI(MPIcl)
    clusterSize(MPIcl)

    all_packages <- c('bigsnpr','doParallel')

    imputed_chr_info <- foreach(i = 1:total_len, .packages=all_packages, .combine=comb, .multicombine=TRUE) %dopar%{ 
        imputeChr_GPU(X, X2, infos.imp, ind.chrs[[i]], alpha, size, p.train, n.cor, seed, ncores, cal_mode)
    }
    closeCluster(MPIcl)
    # replace X2
    print(dim(imputed_chr_info$impute_chromesome))
    print(dim(imputed_chr_info$impute_info))
    
    X2[] <- imputed_chr_info$impute_chromesome
    infos.imp[] <- imputed_chr_info$impute_info
    infos.imp
}


# === main ===

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
print(rds)
if(file.exists(rds)==FALSE) rds <- snp_readBed(file_path, backingfile = temp_path)
print(paste('rds path: ', rds))

# get bed file infos
obj<- snp_attach(rds)
G <- obj$genotypes
CHR <- obj$map$chromosome # ind.chr

NROWS <- dim(G)[1]
NCOLS <- dim(G)[2]
CHR_LEN <- length(CHR)
print(paste('CHR length: ', CHR_LEN))
print(paste('Matrix size: ', NROWS, ' x ', NCOLS))
print(paste('Core: ', NCORES))

# prevent second level parallelization occurs
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
}else{
    a = Sys.time()
    print(paste('Now using CPU version snpfastImpute'))
    G2 <- snp_fastImpute_GPU(G, CHR, ncores = NCORES, seed = SEED, cal_mode = 'hist')
    print(Sys.time()-a)
}
print(paste(' = = = = check imputed infos.imp= = = = '))
print(G2[,1:10])

print(paste('done'))
mpi.quit()








