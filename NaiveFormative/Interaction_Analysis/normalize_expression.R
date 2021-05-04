library(edgeR)


read.rna.seq.data <- function(rna.seq.list, field) {
  
  # Read in file containing RNA-seq results list.
  rna.seq.files <- read.table(rna.seq.list, sep="\t", stringsAsFactors=F, header=F)
  
  # Read in RNA-seq results by cell type.
  cell.types <- unique(rna.seq.files[, 2])
  rna.seq.summary <- list()
  for (cell.type in cell.types) {
    
    # Read in raw data for all replicates corresponding to the current cell type.
    subset.files <- rna.seq.files[rna.seq.files[, 2] == cell.type, ]
    raw.data <- list()
    num.replicates <- length(subset.files[, 1])
    for (i in 1:num.replicates) {
      raw.data[[i]] <- read.table(paste0( "RSEM/", subset.files[i, 1]), sep="\t", stringsAsFactors=F, header=T, 
                                  colClasses=c("character", "character", rep("numeric", 5)))
    }
    
    # Get list of gene IDs. Rows of the RSEM results files should follow the same order of gene IDs.
    gene.ids <- unlist(strsplit(raw.data[[1]]$gene_id, split="\\."))[(1:length(raw.data[[1]]$gene_id))*2-1]
    num.genes <- length(gene.ids)
    
    # Create and populate results table for each replicate. Rows are gene IDs, columns are replicates.
    rna.seq.summary[[cell.type]] <- data.frame(matrix(NA, num.genes, num.replicates))
    rownames(rna.seq.summary[[cell.type]]) <- gene.ids
    colnames(rna.seq.summary[[cell.type]]) <- c(paste0(cell.type, "_", 1:num.replicates))
    for (i in 1:num.replicates) {
      
      if (all(unlist(strsplit(raw.data[[i]]$gene_id, split="\\."))[(1:length(raw.data[[i]]$gene_id))*2-1] == gene.ids)) {
        rna.seq.summary[[cell.type]][, i] <- raw.data[[i]][, field]
      } else {
        print("ERROR: Rows of the RSEM results files should follow the same order of gene IDs.")
      }
      
    }
    
  }
  
  # Return RNA-seq results.
  return(rna.seq.summary)
  
}


# Calculate mean TPM or TMM-normalized FPKM expression values for each gene and cell type.
# Notes: RNA-seq results for each cell type should have the same # and ordering of genes.
#
summarizeExpressionResults <- function(expression.results, gene.lengths=NULL, type) {
  
  # Create results matrix where rows are genes and columns are cell types.
  cell.types <- names(expression.results)
  num.genes <- length(expression.results[[cell.types[1]]][, 1])
  expression.summary <- data.frame(matrix(NA, num.genes, length(cell.types)), stringsAsFactors=F)
  colnames(expression.summary) <- cell.types
  rownames(expression.summary) <- rownames(expression.results[[cell.types[1]]])
  
  # Calculate the specified metric for gene expression.
  if (type == "TPM") {
    
    # Calculate the mean TPM across all replicates.
    for (i in 1:length(cell.types))
      expression.summary[, i] <- rowMeans(expression.results[[cell.types[i]]])
    
  } else if (type == "RPKM") {
    
    # Count number of replicates across all cell types.
    num.replicates <- 0
    for (cell.type in cell.types)
      num.replicates <- num.replicates + length(names(expression.results[[cell.type]]))
    
    # Compile matrices containing gene lengths and expected counts for each replicate across all cell types.
    all.results <- data.frame(matrix(NA, num.genes, num.replicates))
    all.lengths <- data.frame(matrix(NA, num.genes, num.replicates))
    rownames(all.results) <- rownames(expression.results[[cell.types[1]]])
    rownames(all.lengths) <- rownames(expression.results[[cell.types[1]]])
    groups <- c()
    col <- 1
    for (cell.type in cell.types) {
      
      all.results[, col:(col + length(names(expression.results[[cell.type]])) - 1)] <- expression.results[[cell.type]]
      all.lengths[, col:(col + length(names(expression.results[[cell.type]])) - 1)] <- gene.lengths[[cell.type]]
      colnames(all.results)[col:(col + length(names(expression.results[[cell.type]])) - 1)] <- names(expression.results[[cell.type]])
      colnames(all.lengths)[col:(col + length(names(expression.results[[cell.type]])) - 1)] <- names(expression.results[[cell.type]])
      groups <- c(groups, rep(cell.type, length(expression.results[[cell.type]])))
      col <- col + length(names(expression.results[[cell.type]]))
      
    }
    
    # Calculate TMM normalization factors for each replicate using edgeR.
    y <- DGEList(counts=all.results, group=groups)
    y <- calcNormFactors(y, method="TMM")
    
    # Calculate the TMM-normalized RPKM for each replicate (accounting for its reported gene lengths) using edgeR.
    all.rpkms <- data.frame(matrix(NA, num.genes, num.replicates))
    for (i in 1:length(all.results)) {
      
      RPKM <- rpkm(y, gene.length=all.lengths[, i])
      all.rpkms[, i] <- RPKM[, i]
      
    }

    # Calculate the mean TMM-normalized RPKM across all replicates.
    for (i in 1:length(cell.types))
      expression.summary[, i] <- rowMeans(as.matrix(all.rpkms[, which(groups==cell.types[i])]))
    
  }
  
  # Return summarized expression results.
  return(expression.summary)
  
}


transcript.tpm.individual<-read.rna.seq.data('RSEM/rna_seq_output_replicates.txt',"TPM")
transcript.counts.individual<-read.rna.seq.data('RSEM/rna_seq_output_replicates.txt',"expected_count")
transcript.lengths.individual<-read.rna.seq.data('RSEM/rna_seq_output_replicates.txt',"length")

save(transcript.lengths.individual, file="files/transcript.lengths.individual.Rdata")
save(transcript.counts.individual, file="files/transcript.counts.individual.Rdata")
save(transcript.tpm.individual, file="files/transcript.tpm.individual.Rdata")
load(file="files/transcript.lengths.individual.Rdata")
load(file="files/transcript.counts.individual.Rdata")
load(file="files/transcript.tpm.individual.Rdata")


# Calculate mean TPM or TMM-normalized FPKM expression values for each gene and cell type.
tpm.data.individual <- summarizeExpressionResults(transcript.tpm.individual, type="TPM")
rpkm.data.individual <- summarizeExpressionResults(transcript.counts.individual, transcript.lengths.individual, type="RPKM")


# Write results to file.
write.table(tpm.data.individual, file="RSEM/tpm.data.individual.txt", sep="\t", row.names=T, col.names=T, quote=F)
write.table(rpkm.data.individual, file="RSEM/rpkm.data.individual.txt", sep="\t", row.names=T, col.names=T, quote=F)
