
########################## usual functions #########################
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))
movetolast <- function(data, move) {  data[c(setdiff(names(data), move), move)] }
normalization_01 <- function(x) { x <- (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) ; return(x) }
percent <- function(x, digits = 0, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}
## vapply(strsplit( vector , split="[_]"), "[", "", 1) take first element of a split 
## read_excel:  library(readxl) library(openxlsx)
## Reduce(intersect, list(a,b,c))
## gsub("\\ ",".",names(dataset))  # substitute text/string
## aggregate( .~ gene, data = GEX, mean)
## matches <- unique(grep(paste(toMatch,collapse="|"),myfile$Letter,value=TRUE))
## t(do.call("rbind",HiRes_list))
## features <- data.frame(lapply(features, as.character), stringsAsFactors=FALSE)
## features <- data.frame(lapply(features, as.numeric), stringsAsFactors=FALSE)

## library(stringr)
## str_sub(x, end= 2)             # take first 2 strings
## str_sub(x, start= 3)           # remove first 2 strings
## str_sub(x, start= -2)          # take last 2 strings
## str_sub(x, end= -3)            # remove last 2 strings
 

## yyz[] <- lapply(yyz, function(x) as.numeric(as.character(x)))

## from matrix of character to matrix of numeric
## df <- matrix(as.numeric(df),ncol = ncol(df))

###### from table to matrix
# makes sure there is no duplicates first.
# mat_table_unique <- mat_table[!duplicated(mat_table[c(1,2,3)]),] # columns 1,2,3 correspond to value, ind, Var1
# df <- xtabs(value ~ ind + Var1, mat_table)
# df <-  as.data.frame.matrix(df)

###### ID conversion
# from Uniprot to HGNC 
# library(AbHAC) ; uniprot.to.hgnc(uniprot, id.con.set = id.conversion.set)

###### list of dataframes
# lapply( df_list, "[", , names )  # extract specific columns from a list of data frames

###### count Missing/Specific values in dataframe
# NA_per_row <- rowSums(is.na(df) | df == "")  # NA per row in df 

###### Directories
# dir.create(file.path(result_folder), recursive = T, showWarnings = F)

###### Divers
# split vector into chunks of 5
# split(my_vec, ceiling(seq_along(my_vec) / 5))

# library(impute) : GEX_knn <- impute.knn(data.matrix(GEX),k = 10)$data

## add specific numbers to each row
#mat <- matrix(1:15, nrow=5)
#sweep(mat, 1, c(5, 10, 15, 20, 25), "+")

convert_ENS_HGNC_df_row <- function(df) {
  conversion_table <- read.csv("/Users/i0535027/Documents/Sanofi/ID_CONVERSION/table_conversion_ENS_HGNC_ENTREZ.txt",sep="\t")
  conversion_table <- conversion_table[conversion_table$Ensembl.Gene.ID %in% rownames(df), ]
  df <- df[conversion_table$Ensembl.Gene.ID, ]
  rownames(df) <- conversion_table$Approved.Symbol
  return(df)
}

convert_ENTREZ_HGNC_df_row <- function(df) {
  conversion_table <- read.csv("/Users/i0535027/Documents/Sanofi/ID_CONVERSION/table_conversion_ENS_HGNC_ENTREZ.txt",sep="\t")
  conversion_table <- conversion_table[conversion_table$Entrez.Gene.ID %in% rownames(df), ]
  df <- df[conversion_table$Entrez.Gene.ID, ]
  rownames(df) <- conversion_table$Approved.Symbol
  return(df)
}

convert_protein_product_df_row <- function(df) {
  HGNC_gene_with_protein_product <- read.delim("/Users/i0535027/Documents/Sanofi/ID_CONVERSION/HGNC_gene_with_protein_product.txt")
  common <- intersect(rownames(df),HGNC_gene_with_protein_product$symbol)
  HGNC_gene_with_protein_product <- HGNC_gene_with_protein_product[HGNC_gene_with_protein_product$symbol %in% common , ]
  df <- df[HGNC_gene_with_protein_product$symbol , ]
  return(df)
}

progeny <- function(GEX) {
  model <- data.matrix( read.csv("/Users/miyang/Documents/PrognomIQ/GENERAL_DATA/model_14PW.csv", row.names=1) )
  common <- intersect( rownames(GEX) , rownames(model) ) 
  GEX <- GEX[common, ]
  model <- model [common, ]
  PROGENy <- t(GEX) %*% model
  return(t(PROGENy))
}

ssGSEA <-function(GEX, genesets) {
  #################################### ssGSEA ####################################
  library(GSVA) ; library(GSA)
  # canonical_pathway_1329  immunologic_signatures_1888
  MSigDB <- GSA.read.gmt(paste0(path,"/GENERAL_DATA/MSigDB/",genesets,".gmt"))
  geneSets <- MSigDB$genesets ; names(geneSets) <- MSigDB$geneset.names
  GEX_ssGSEA <- t(gsva(data.matrix(t(GEX)), geneSets, method="ssgsea", verbose=FALSE, parallel.sz=4))
  return(GEX_ssGSEA)
}

remove_string <- function(string) {
  library(stringr)
  string <- str_remove(string=string, pattern="REACTOME_")
  string <- str_remove(string=string, pattern="BIOCARTA_")
  string <- str_remove(string=string, pattern="KEGG_")
  string <- str_remove(string=string, pattern="PID_")
  string <- str_remove(string=string, pattern="SIGNALING_BY_")
  string <- str_remove(string=string, pattern="REGULATION_OF_")
  string <- str_remove(string=string, pattern="CANONICAL_")
  string <- str_remove(string=string, pattern="_MEDIATED_ANTIGEN_PROCESSING")
  string <- str_remove(string=string, pattern="INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_")
  
  string <- str_remove(string=string, pattern="ST_")
  return(string)
}


load_folder <- function (workdir ) {
  setwd(workdir)
  files <- list.files(path=getwd() , pattern="") # 
  n <- 0
  imp <- list()
  for(file in files)
  {
    n <- n + 1
    imp[[n]] <- read.csv(paste(workdir, "/",files[n], sep = ""), row.names = 1 , check.names = T)
  }
  names(imp) <- files
  return(imp)
}

stopQuietly <- function(...) {
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" ")); stop(simpleError(blankMsg));
} # stopQuietly()

ClosestMatch2 = function(string, stringVector, maxDist=1){ 
  library(dplyr)
  library(stringdist)
  stringVector[amatch(string, stringVector, maxDist=maxDist)] 
}

ID_to_celllines <- function(x) {
  
  MASTER_LIST <- as.matrix(MASTER_LIST)
  MASTER_LIST[753,2] <- "PC-3-[JPC-3]"
  M <- MASTER_LIST[ ,1:2]
  names <- c()
  for (i in 1:length(rownames(x))) {
    if(length(grep(rownames(x)[i], M[,1])) == 1){
      names <- c(names, as.character(M[,2][grep(rownames(x)[i], M[,1])]) ) } else {names <- c(names, rownames(x)[i])}
  }
  rownames(x) <- names
  return (x)
}

celllines_to_ID <- function(x) {
  
  MASTER_LIST <- as.matrix(MASTER_LIST)
  MASTER_LIST[753,2] <- as.character(MASTER_LIST[753,2])
  MASTER_LIST[753,2] <- "PC-3-[JPC-3]"
  M <- MASTER_LIST[ ,1:2]
  names <- c()
  for (i in 1:length(rownames(x))) {
    if(length(grep(rownames(x)[i], M[,2])) == 1){
      names <- c(names, M[,1][grep(rownames(x)[i], M[,2])] ) } else {names <- c(names, rownames(x)[i])
      }
  }
  names <- trimws(names)
  rownames(x) <- names
  return (x)
}

ID_to_celllines_RNAseq <- function(x) {
  CLI <- CLI[ ,c(16, 3, 11)]
  M <- CLI[ ,1:2]
  names <- c()
  for (i in 1:length(rownames(x))) {
    if(length(grep(rownames(x)[i], M[,1])) == 1){
      names <- c(names, as.character(M[,2][grep(rownames(x)[i], M[,1])]) ) } else {names <- c(names, rownames(x)[i])}
  }
  rownames(x) <- names
  return (x)
}

Hugo_to_EntrezID <- function(gene) {
  library(biomaRt)
  ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  mapTab <- getBM(attributes = c("hgnc_symbol", "entrezgene"), filters = "hgnc_symbol", values = gene, mart = ensembl, uniqueRows=FALSE)
  dupRows <- union(which(duplicated(mapTab[,1])), which(duplicated(mapTab[,2])))
  entrezIds <- mapTab[-dupRows, 2]
  names(entrezIds) <- mapTab[-dupRows, 1]
  entrezIds <- entrezIds[!is.na(entrezIds)]
  return(entrezIds)
}

EntrezID_to_Hugo <- function(gene) {
  library(biomaRt)
  ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  mapTab <- getBM(attributes = c("entrezgene", "hgnc_symbol"), filters = "entrezgene", values = gene, mart = ensembl, uniqueRows=FALSE)
  dupRows <- union(which(duplicated(mapTab[,1])), which(duplicated(mapTab[,2])))
  hgnc <- mapTab[-dupRows, 2]
  names(hgnc) <- mapTab[-dupRows, 1]
  hgnc <- hgnc[!is.na(hgnc)]
  return(hgnc)
}

ENS_to_Hugo <- function(gene) {
  library(biomaRt) ## mart = useMart('ensembl') ; listDatasets(mart) ; BiocInstaller::biocLite('grimbough/biomaRt')
  ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") ## "hsapiens_gene_ensembl"
  mapTab <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = gene, mart = ensembl, uniqueRows=FALSE)
  dupRows <- union(which(duplicated(mapTab[,1])), which(duplicated(mapTab[,2])))
  hgnc <- mapTab[-dupRows, 2]
  names(hgnc) <- mapTab[-dupRows, 1]
  hgnc <- hgnc[!is.na(hgnc)]
  return(hgnc)
}

get_ensembl_table <- function() {
  library(biomaRt)
  ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  mapTab <- getBM(attributes = c("ensembl_gene_id","description","hgnc_symbol"),mart = ensembl, uniqueRows=FALSE) # listAttributes(ensembl)
  return(mapTab)
}

corr_by_row <- function(mat1, mat2, method="pearson")  {
  mat1 <- as.matrix(mat1)
  mat2 <- as.matrix(mat2) 

  corr_vec <- c()
  for(i in 1:length(mat1[,1]) ) {
    
    c <- rbind(mat1[i, ], mat2[i, ]) ; c <- c[ ,complete.cases(c)]
    if(length(which(apply(c, 1, var) == 0)) > 0) { corr_vec <- c(corr_vec , 0 ) } else 
      {
      temp <- cor.test(mat1[i, ], mat2[i, ],method=method)
      pcorr <- temp$estimate # pearson correlation
      corr_vec <- c(corr_vec , pcorr)
      }
  }
  boxplot(corr_vec ) ; print(paste("mean: ", mean(corr_vec, na.rm = T), sep = "") )   ; print(paste("median: ", median(corr_vec, na.rm = T), sep = "")  ) 
  names(corr_vec) <- rownames(mat1)
  return(corr_vec)
}


rmse_by_row <- function(mat1, mat2)  {
  mat1 <- as.matrix(mat1)
  mat2 <- as.matrix(mat2) 
  
  rmse_vec <- c()
  for(i in 1:length(mat1[,1]) ) {
    rmse_vec <- c(rmse_vec , rmse(mat1[i, ], mat2[i, ], na.rm=TRUE))
  }
  boxplot(rmse_vec ) ; print(paste("mean: ", mean(rmse_vec, na.rm = T), sep = "") )   ; print(paste("median: ", median(rmse_vec, na.rm = T), sep = "")  ) 
  names(rmse_vec) <- rownames(mat1)
  return(rmse_vec)
}

common_full <- function(MAT_LIST)  {
  for (iteration in 1:2) {
  row_list <- list()
  col_list <- list()
  for(i in 1:length(MAT_LIST)) {
    row_list[[i]] <- rownames(MAT_LIST[[i]])
    col_list[[i]] <- colnames(MAT_LIST[[i]])
  }
  common_row <- Reduce(intersect, row_list )
  common_col <- Reduce(intersect, col_list )
  for(i in 1:length(MAT_LIST)) {
    MAT_LIST[[i]] <-MAT_LIST[[i]][common_row, common_col]
    MAT_LIST[[i]] <- MAT_LIST[[i]][rowSums(is.na(MAT_LIST[[i]])) < length(MAT_LIST[[i]][1,]) - 2 , ] # remove rows with too many NAs
   }
  }
  return(MAT_LIST)
}

common_row <- function(MAT_LIST)  {
  MAT_LIST_output <- list()
  for (iteration in 1:2) { # iteration=1
    row_list <- list()
    for(i in 1:length(MAT_LIST)) {
      row_list[[i]] <- rownames(MAT_LIST[[i]])
    }
    common_row <- Reduce(intersect, row_list )
    for(i in 1:length(MAT_LIST)) {
      MAT_LIST_output[[i]] <-MAT_LIST[[i]][common_row, ]
      # MAT_LIST_output[[i]] <- MAT_LIST_output[[i]][rowSums(is.na(MAT_LIST_output[[i]])) < length(MAT_LIST_output[[i]][1,]) - 2 , ] # remove rows with too many NAs
    }
  }
  return(MAT_LIST_output)
}

# MAT_LIST <- list(LM22, clinical_features)
common_row_bind <- function(MAT_LIST)  {
  for (iteration in 1:1) {
    row_list <- list()
    for(i in 1:length(MAT_LIST)) {
      row_list[[i]] <- rownames(MAT_LIST[[i]])
    }
    common_row <- Reduce(intersect, row_list )
    for(i in 1:length(MAT_LIST)) {
      MAT_LIST[[i]] <- MAT_LIST[[i]][common_row, ]
      # MAT_LIST[[i]] <- MAT_LIST[[i]][rowSums(is.na(MAT_LIST[[i]])) < length(MAT_LIST[[i]][1,]) - 2 , ] # remove rows with too many NAs
    }
  }
  
  MAT_LIST_bind <- c()
  for(i in 1:length(MAT_LIST)) {
    if(i == 1) { MAT_LIST_bind <- MAT_LIST[[i]] } else {
      MAT_LIST_bind <- cbind(MAT_LIST_bind, MAT_LIST[[i]])
    }
  }
  return(MAT_LIST_bind)
}

common_col <- function(MAT_LIST)  {
  for (iteration in 1:1) {
    col_list <- list()
    for(i in 1:length(MAT_LIST)) {
      col_list[[i]] <- colnames(MAT_LIST[[i]])
    }
    common_col <- Reduce(intersect, col_list )
    for(i in 1:length(MAT_LIST)) {
      MAT_LIST[[i]] <-MAT_LIST[[i]][ , common_col ]
      # MAT_LIST[[i]] <- MAT_LIST[[i]][rowSums(is.na(MAT_LIST[[i]])) < length(MAT_LIST[[i]][1,]) - 2 , ] # remove rows with too many NAs
    }
  }
  return(MAT_LIST)
}

common_col_bind <- function(MAT_LIST)  {
  for (iteration in 1:2) {
    col_list <- list()
    for(i in 1:length(MAT_LIST)) {
      col_list[[i]] <- colnames(MAT_LIST[[i]])
    }
    common_col <- Reduce(intersect, col_list )
    for(i in 1:length(MAT_LIST)) {
      MAT_LIST[[i]] <-MAT_LIST[[i]][ , common_col ]
      MAT_LIST[[i]] <- MAT_LIST[[i]][rowSums(is.na(MAT_LIST[[i]])) < length(MAT_LIST[[i]][1,]) - 2 , ] # remove rows with too many NAs
    }
  }
  
  MAT_LIST_bind <- c()
  for(i in 1:length(MAT_LIST)) {
    if(i == 1) { MAT_LIST_bind <- MAT_LIST[[i]] } else {
      MAT_LIST_bind <- rbind(MAT_LIST_bind, MAT_LIST[[i]])
    }
  }
  return(MAT_LIST_bind)
}




find_pubmed_association <- function (vector_row, vector_column, constant_term="") {
library(RISmed)
m <- matrix(nrow = length(vector_row), ncol = length(vector_column))
rownames(m) <- vector_row ; colnames(m) <- vector_column
for(i in 1:length(vector_row)) {
  for(j in 1:length(vector_column)) {
    tryCatch({
      query <- paste("(",vector_row[i]," AND ",vector_column[j]," AND ",constant_term,")", sep = "" )
      mindate=2008 ; maxdate=2018
      ngs_search <- EUtilsSummary(query, type="esearch",db = "pubmed",mindate=mindate, maxdate=maxdate, retmax=500)
      m[i,j] <- QueryCount(ngs_search)
    }, error=function(e){})
  }
}
return(m)
}


binarize_by_column <- function(mat, th, binary=F) {
  for(i in 1:length(colnames(mat)) ){ # i=2
    limit <- quantile(mat[ ,i], th)
    mat[ ,i][mat[ ,i]<=limit]  <- "low" 
    mat[ ,i][mat[ ,i] %not in% "low"] <- "high" 
  }
  if(binary==T) {
    mat[mat=="high"] <- 1
    mat[mat=="low" ] <- 0
  }
  return(mat)
}

binarize_columns <- function(mat,column, th, binary=F) {
  for(i in column ){ # i=2
    limit <- quantile(mat[ ,i], th)
    mat[ ,i][mat[ ,i]<=limit]  <- "low" 
    mat[ ,i][mat[ ,i] %not in% "low"] <- "high" 
  }
  if(binary==T) {
    mat[mat=="high"] <- 1
    mat[mat=="low" ] <- 0
  }
  return(mat)
}

# mat=clinical ; column=1 ; th=0.3 ; binary=F
binarize_columns_reduce <- function(mat,column, th, binary=F) {
  # here th mean take top X% or down X%
  for(i in column ){ # i=1
    limit_down <- quantile(mat[ ,i], th)
    limit_up <- quantile(mat[ ,i], 1-th )
    
    mat[ ,i][mat[ ,i]<=limit_down]  <- "low" 
    mat_down <- mat[ mat[ ,i] == "low" , ]
    
    mat <- mat[ mat[ ,i] != "low" , ]
    
    mat[ ,i][mat[ ,i]>=limit_up]  <- "high" 
    mat_high <- mat[ mat[ ,i] == "high" , ]
    
    mat = rbind(mat_high,mat_down)
    
  }
  if(binary==T) {
    mat[mat=="high"] <- 1
    mat[mat=="low" ] <- 0
  }
  return(mat)
}

convert_pathway <- function(GEX,geneset_name,N_core=1) {
  library(GSVA) ; library(GSA)
  MSigDB <- GSA.read.gmt(paste0(path,"/GENERAL_DATA/MSigDB/",geneset_name,".gmt"))
  geneSets <- MSigDB$genesets ; names(geneSets) <- MSigDB$geneset.names
  GEX <- gsva(data.matrix(GEX), geneSets, method="ssgsea", verbose=FALSE, parallel.sz=N_core)
  return(GEX)
}

