FBM_infos <- function(Gna) {

  base <- sub_bk(Gna$backingfile, "-infos-impute", stop_if_not_ext = FALSE)
  rds <- paste0(base, ".rds")

  if (file.exists(rds)) {
    # cat("Attach\n")
    readRDS(rds)
  } else {
    # cat("Create\n")
    FBM(2, ncol(Gna), backingfile = base, init = NA_real_)$save()
  }
}

#=========

library(bigsnpr)
bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
rds <- snp_readBed(bedfile, backingfile = "a")

obj<- snp_attach(rds)
G <- obj$genotypes
CHR <- obj$map$chromosome # ind.chr
NCORES <- 1
G2 <- snp_fastImpute(G, CHR,NCORES)

# testing snp_corr running time
CODE_IMPUTE_LABEL <- c(0, 0.5, 1, rep(NA, 253))
CODE_IMPUTE_PRED  <- c(0, 1, 2, NA, 0, 1, 2, rep(NA, 249))

G_X <- G$copy(code = CODE_IMPUTE_LABEL)
G_X2 <- G$copy(code = CODE_IMPUTE_PRED)
infos.imp <- FBM_infos(G)

n <- nrow(G_X)
# use the following code if Error: Two levels of parallelism are used
options(bigstatsr.check.parallel.blas = FALSE)
corr <- snp_cor(
  Gna       = G_X,
  ind.row   = sort(sample(n, size = nrow(G))),
  ind.col   = CHR,
  size      = 200,
  alpha     = 1e-04,
  fill.diag = FALSE,
  ncores    = 2
)

library(rbenchmark)
benchmark(
  'snp_cor_1' = {
    snp_cor(
      Gna       = G_X,
      ind.row   = sort(sample(n, size = nrow(G))),
      ind.col   = CHR,
      size      = 200,
      alpha     = 1e-04,
      fill.diag = FALSE,
      ncores    = 1
    )
  },
  'snp_cor_4' = {
    snp_cor(
      Gna       = G_X,
      ind.row   = sort(sample(n, size = nrow(G))),
      ind.col   = CHR,
      size      = 200,
      alpha     = 1e-04,
      fill.diag = FALSE,
      ncores    = 4
    )
  },
  'snp_cor_8' = {
    snp_cor(
      Gna       = G_X,
      ind.row   = sort(sample(n, size = nrow(G))),
      ind.col   = CHR,
      size      = 200,
      alpha     = 1e-04,
      fill.diag = FALSE,
      ncores    = 8
    )
  },
  'snp_cor_16' = {
    snp_cor(
      Gna       = G_X,
      ind.row   = sort(sample(n, size = nrow(G))),
      ind.col   = CHR,
      size      = 200,
      alpha     = 1e-04,
      fill.diag = FALSE,
      ncores    = 16
    )
  },
  replications = 20,
  columns = c("test", "replications", "elapsed",
              "relative", "user.self", "sys.self")
)

# create new data for testing large dataset (for snp_cor)
NROWS <- 517 * 50 
NCOLS <- 4542
code <- rep(NA_real_, 256)
code[1:3] <- c(0,1,2)

G_toyBig <- FBM.code256(NROWS, NCOLS, code, init = sample(as.raw(0:2), NROWS*NCOLS, replace = TRUE)) 
G_toyBig[] <- sample(c(0,1,2), size = length(G_toyBig), replace = TRUE)
G_toyBig_CHR <- rep(1,NCOLS)

G_toyBig_X <- G_toyBig$copy(code = CODE_IMPUTE_LABEL)

snp_cor(
  Gna       = G_toyBig_X,
  ind.row   = sort(sample(nrow(G_toyBig_X), size = nrow(G_toyBig))),
  ind.col   = G_toyBig_CHR,
  size      = 200,
  alpha     = 1e-04,
  fill.diag = FALSE,
  ncores    = 1
)

benchmark(
  'snp_cor_1_toy' = {
    snp_cor(
      Gna       = G_toyBig_X,
      ind.row   = sort(sample(nrow(G_toyBig_X), size = nrow(G_toyBig))),
      ind.col   = G_toyBig_CHR,
      size      = 200,
      alpha     = 1e-04,
      fill.diag = FALSE,
      ncores    = 1
    )
  },
  'snp_cor_4_toy' = {
    snp_cor(
      Gna       = G_toyBig_X,
      ind.row   = sort(sample(nrow(G_toyBig_X), size = nrow(G_toyBig))),
      ind.col   = G_toyBig_CHR,
      size      = 200,
      alpha     = 1e-04,
      fill.diag = FALSE,
      ncores    = 4
    )
  },
  'snp_cor_8_toy' = {
    snp_cor(
      Gna       = G_toyBig_X,
      ind.row   = sort(sample(nrow(G_toyBig_X), size = nrow(G_toyBig))),
      ind.col   = G_toyBig_CHR,
      size      = 200,
      alpha     = 1e-04,
      fill.diag = FALSE,
      ncores    = 8
    )
  },
  'snp_cor_16_toy' = {
    snp_cor(
      Gna       = G_toyBig_X,
      ind.row   = sort(sample(nrow(G_toyBig_X), size = nrow(G_toyBig))),
      ind.col   = G_toyBig_CHR,
      size      = 200,
      alpha     = 1e-04,
      fill.diag = FALSE,
      ncores    = 16
    )
  },
  replications = 5,
  columns = c("test", "replications", "elapsed",
              "relative", "user.self", "sys.self")
)

