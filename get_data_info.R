# read data
library(bigsnpr)
file_path <- './med_data/sample4.bed'
rds <- './med_temp/temp_file_1.rds'
# temp_path <- './med_temp/tempfile'

if(file.exists(rds)==FALSE){
  rds <- snp_readBed(file_path, backingfile = temp_path)
}
print(paste('rds path: ', rds))

obj<- snp_attach(rds)
G <- obj$genotypes
CHR <- obj$map$chromosome # ind.chr

# get chromesome length
chr_len <- table(CHR)

# get NAs on each chromesome
NAs_chromosome <- c()

head <- 1
end <- 0
for(i in c(1:22)){
  end <- end + chr_len[i]
  get_counts <- big_counts(G, ind.row = 1:dim(G)[1], ind.col = head:end, byrow = FALSE)
  head <- head + chr_len[i]
  now_counts <- sum(get_counts[4,])
  print(now_counts)
  NAs_chromosome[i] <- now_counts
}


# get NA% of each chromesome
NA_ratio <- c()
rol <- dim(G)[1]
for(i in c(1:22)){
  total <- rol * chr_len[i]
  NA_ratio[i] <- NAs_chromosome[i]/total
}

NA_df <- data.frame(Chr_len = as.vector(chr_len), NAs = NAs_chromosome, NA_ratio = NA_ratio)
NA_df


