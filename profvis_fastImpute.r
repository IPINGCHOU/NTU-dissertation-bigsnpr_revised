CODE_IMPUTE_LABEL <- c(0, 0.5, 1, rep(NA, 253))
CODE_IMPUTE_PRED  <- c(0, 1, 2, NA, 0, 1, 2, rep(NA, 249))

NROWS <- 10000
NCOLS <- 10000
code <- rep(NA_real_, 256)
code[1:3] <- c(0,1,2)

# Testing fastImpute elaspe time (1st:xgboost 2:snp_cor)
# Testing fastImpute elaspe time without NAs (snp_cor will still be activated?)

profvis({
  G_toy <- FBM.code256(NROWS, NCOLS, code, init = sample(as.raw(0:2), NROWS*NCOLS, replace = TRUE)) 
  prob_of_NAs = 0.05
  res_of_probs = (1-prob_of_NAs) / 3
#  G_toy[] <- sample(c(0,1,2,3), size = length(G_toy), replace = TRUE, prob = c(res_of_probs,res_of_probs,res_of_probs,prob_of_NAs))
  G_toy[] <- sample(c(0,1,2), size = length(G_toy), replace = TRUE)
  G_toy_CHR <- rep(1,NCOLS)
  
  G_toy_X <- G_toy$copy(code = CODE_IMPUTE_LABEL)
  
  print(paste('The total NAs: ', length(G_toy)*prob_of_NAs))
  NCORES <- 8
  G2_toy <- snp_fastImpute(G_toy_X, G_toy_CHR, ncores = NCORES)
}
)

# Testing snp_fastImpute with different cores
# set a function for remove the info_impute bk and rds file
delete_info_impute <- function(bedfile_name){
  print(paste('Deleting bk rds files from ', bedfile_name))
  Loc_info <- sub_bk(bedfile_name, "-infos-impute", stop_if_not_ext = FALSE)
  bk_info_pos <- paste0(Loc_info, ".bk")
  rds_info_pos <- paste0(Loc_info, ".rds")
  system(paste0("rm -r ", bk_info_pos))
  system(paste0("rm -r ", rds_info_pos))
  print(paste('All bk rds files were deleted'))
}

# Init
CODE_IMPUTE_LABEL <- c(0, 0.5, 1, rep(NA, 253))
CODE_IMPUTE_PRED  <- c(0, 1, 2, NA, 0, 1, 2, rep(NA, 249))
# Matrix size
NROWS <- 1000
NCOLS <- 2000
code <- rep(NA_real_, 256)
code[1:3] <- c(0,1,2)
# Generating testing matrix
G_toy <- FBM.code256(NROWS, NCOLS, code, init = sample(as.raw(0:2), NROWS*NCOLS, replace = TRUE)) 
prob_of_NAs = 0.05
res_of_probs = (1-prob_of_NAs) / 3
G_toy[] <- sample(c(0,1,2,3), size = length(G_toy), replace = TRUE, prob = c(res_of_probs,res_of_probs,res_of_probs,prob_of_NAs))
G_toy_CHR <- rep(1,NCOLS)
G_toy_X <- G_toy$copy(code = CODE_IMPUTE_LABEL)

print(paste('The total NAs: ', length(G_toy)*prob_of_NAs))
# profiling
profvis({
  # Init
  # print(paste('Core: ', 1))
  # snp_fastImpute(G_toy_X, G_toy_CHR, ncores = 1)
  # # delete_info_impute(G_toy_X$backingfile)
  # print(paste('Deleting bk rds files from ', G_toy_X$backingfile))
  # Loc_info <- sub_bk(G_toy_X$backingfile, "-infos-impute", stop_if_not_ext = FALSE)
  # bk_info_pos <- paste0(Loc_info, ".bk")
  # rds_info_pos <- paste0(Loc_info, ".rds")
  # system(paste0("rm -r ", bk_info_pos))
  # system(paste0("rm -r ", rds_info_pos))
  # print(paste('All bk rds files were deleted'))
  # 
  # print(paste('Core: ', 2))
  # snp_fastImpute(G_toy_X, G_toy_CHR, ncores = 2)
  # # delete_info_impute(G_toy_X$backingfile)
  # print(paste('Deleting bk rds files from ', G_toy_X$backingfile))
  # Loc_info <- sub_bk(G_toy_X$backingfile, "-infos-impute", stop_if_not_ext = FALSE)
  # bk_info_pos <- paste0(Loc_info, ".bk")
  # rds_info_pos <- paste0(Loc_info, ".rds")
  # system(paste0("rm -r ", bk_info_pos))
  # system(paste0("rm -r ", rds_info_pos))
  # print(paste('All bk rds files were deleted'))
  # 
  # print(paste('Core: ', 4))
  # snp_fastImpute(G_toy_X, G_toy_CHR, ncores = 4)
  # # delete_info_impute(G_toy_X$backingfile)
  # print(paste('Deleting bk rds files from ', G_toy_X$backingfile))
  # Loc_info <- sub_bk(G_toy_X$backingfile, "-infos-impute", stop_if_not_ext = FALSE)
  # bk_info_pos <- paste0(Loc_info, ".bk")
  # rds_info_pos <- paste0(Loc_info, ".rds")
  # system(paste0("rm -r ", bk_info_pos))
  # system(paste0("rm -r ", rds_info_pos))
  # print(paste('All bk rds files were deleted'))
  # 
  # print(paste('Core: ', 8))
  # snp_fastImpute(G_toy_X, G_toy_CHR, ncores = 8)
  # # delete_info_impute(G_toy_X$backingfile)
  # print(paste('Deleting bk rds files from ', G_toy_X$backingfile))
  # Loc_info <- sub_bk(G_toy_X$backingfile, "-infos-impute", stop_if_not_ext = FALSE)
  # bk_info_pos <- paste0(Loc_info, ".bk")
  # rds_info_pos <- paste0(Loc_info, ".rds")
  # system(paste0("rm -r ", bk_info_pos))
  # system(paste0("rm -r ", rds_info_pos))
  # print(paste('All bk rds files were deleted'))
  
  print(paste('Core: ', 16))
  snp_fastImpute(G_toy_X, G_toy_CHR, ncores = 16)
  # delete_info_impute(G_toy_X$backingfile)
  # print(paste('Deleting bk rds files from ', G_toy_X$backingfile))
  # Loc_info <- sub_bk(G_toy_X$backingfile, "-infos-impute", stop_if_not_ext = FALSE)
  # bk_info_pos <- paste0(Loc_info, ".bk")
  # rds_info_pos <- paste0(Loc_info, ".rds")
  # system(paste0("rm -r ", bk_info_pos))
  # system(paste0("rm -r ", rds_info_pos))
  # print(paste('All bk rds files were deleted'))
  
  }
)


.rs.restartR()


