#
# gwas_mining.R
# Author: Michael Song
# Last modified: 2038-01-21
# This script mines SNPs from the GWAS Catalog, then imputes them and formats them for downstream analysis.
#


# Load dependencies.
require(haploR)
require(GenomicRanges)
require(rtracklayer)
require(rsnps)
require(dplyr)

# Set up environment and load settings ------------------------------------


# Clear workspace before running script.
rm(list=ls())

# Turn off scientific notation for writing output.
options(scipen=999)

# Convert warnings into errors.
# options(warn=2)

# Set home directory.
home.dir <- '/shen/shenlabstore3/neng/20200509_scATAC_sneddon/SNP_analysis'


# Load auxiliary functions.
scripts.dir <- paste0(home.dir, "/scripts")

exons<-read.table(gzfile(paste0(home.dir,'/resources/hg38_gencode_v32_exons.bed.gz'),'rt'),header=F,stringsAsFactors=F)

source(paste0(scripts.dir,"/utilities.R"))

# Define the significance cutoff for filtering SNPs from the GWAS Catalog.
sig.cutoff <- 10^-6

# Define the LD threshold for performing imputation.
ld.threshold <- 0.8

# Define the populations to be used for imputation.
populations <- c("AFR", "AMR", "ASN", "EUR")

# List of diseases to mine SNPs for, as well ass the dates where the associations and study metadata were downloaded from the GWAS Catalog.
diseases <- c("T1D", "T2D")
dates <- rep("2018-10-30", 2)


# Load resources ----------------------------------------------------------


# Read in the RsMergeArch table for updating rsids.
rs.merge.file <- paste0(home.dir, "/resources/RsMergeArch.bcp")
rs.merge <- read.table(rs.merge.file, sep="\t", header=F, stringsAsFactors=F)
rs.merge[, 1] <- paste0("rs", rs.merge[, 1])
rs.merge[, 2] <- paste0("rs", rs.merge[, 2])


# SNP imputation ----------------------------------------------------------


# Read in study and association data for each disease and perform imputation by population(s).
studies <- list()
associations <- list()
snps.by.population <- list()
imputation.results <- list()
for (i in 1:length(diseases)) {
  
  # Set the current disease.
  current.disease <- diseases[i]
  
  # Read in study and association data for each disease.
  print(paste0("Imputing SNPs for: ", current.disease))
  studies[[current.disease]] <- read.table(paste0(home.dir, "/", current.disease, "_", dates[i], "_studies.txt"), 
                                           sep="\t", header=T, stringsAsFactors=F, quote="",
                                           colClasses=c(rep("character", 9), "integer", "character", "integer"))
  associations[[current.disease]] <- read.table(paste0(home.dir, "/", current.disease, "_", dates[i], "_withChildTraits.tsv"), 
                                                sep="\t", header=T, stringsAsFactors=F, fill=T, quote="")
  
  # Print #s of studies and associations.
  print(paste0("Number of studies: ", length(studies[[current.disease]][, 1])))
  print(paste0("Number of associations: ", sum(studies[[current.disease]]$Association.count)))
  print("Studies with associations that also have metadata available:")
  print(table(unique(associations[[current.disease]]$STUDY.ACCESSION) %in% studies[[current.disease]]$Study.accession))
  # print("Studies with associations without metadata available:")
  # print(unique(associations[[current.disease]]$STUDY.ACCESSION[!(associations[[current.disease]]$STUDY.ACCESSION %in% studies[[current.disease]]$Study.accession)]))
  
  # Filter out studies unrelated to the current disease.
  # studies[[current.disease]] <- studies[[current.disease]][studies[[current.disease]]$Include == 1, ]
  # print(paste0("Number of studies used for imputation: ", length(unique(studies[[current.disease]]$Study.accession))))
  
  # Filter SNPs according to the significance cutoff.
  snps.filtered <- associations[[current.disease]][as.numeric(associations[[current.disease]]$P.VALUE) < sig.cutoff, ]
  print(paste0("Number of SNPs passing significance cutoff of ", sig.cutoff, ": ", length(snps.filtered[, 1])))
  
  # Count how many SNPs do not have a current rsid ('SNP_ID_CURRENT' field is NA).
  snps.current.ids <- snps.filtered$SNP_ID_CURRENT
  snps.na.ids <- which(is.na(snps.filtered$SNP_ID_CURRENT))
  print(paste0("Number of SNPs without a current rsid ('SNP_ID_CURRENT' field is NA): ", length(snps.na.ids)))
  
  # Detect instances where the rsid in the 'SNPS' field doesn't agree with the 'SNP_ID_CURRENT' field. Print out all such instances.
  # 'SNP_ID_CURRENT' takes precedence over 'SNPS', so use the former as long as it is available.
  snps.reported <- snps.filtered$SNPS
  num.disagreements <- length(which(paste0("rs", snps.current.ids) != snps.reported))
  print(paste0("Number of instances where the rsid in the 'SNPS' field doesn't agree with the 'SNP_ID_CURRENT' field: ", num.disagreements))
  if (num.disagreements > 0) {
    
    which.disagreements <- !(paste0("rs", snps.current.ids) == snps.reported)
    mismatches <- cbind(snps.current.ids[which.disagreements], snps.reported[which.disagreements], is.na(snps.filtered$SNP_ID_CURRENT)[which.disagreements])
    colnames(mismatches) <- c("SNP_ID_CURRENT", "SNPS", "IS.NA")
     print(mismatches)
    
  }

  # Try to use the information in the 'SNPS' field to resolve instances where the 'SNP_CURRENT_ID' is NA. Keep track of which instances were resolved.
  resolved <- rep(FALSE, length(snps.na.ids))
  if (length(snps.na.ids) > 0) {
    
    # Try to resolve each NA case one by one.
    for (j in 1:length(snps.na.ids)) {
      
      # Set current query to try to resolve.
      query <- snps.reported[snps.na.ids[j]]
      
      # First check if it is a delimited lists of multiple rsids by splitting via comma, then by semicolon.
      query.split <- unlist(strsplit(query, split=", "))
      if (length(query.split) == 1) {
        query.split <- unlist(strsplit(query, split="; "))
      }
      
      # If splitting resulted in multiple entries, process them--otherwise, turn to the translation table.
      if (length(query.split) > 1) {
        
        # Remove "rs" prefix for conformity, then update the 'SNP_ID_CURRENT' field (the "#" character is used as a delimiter in this script).
        query.split <- gsub("rs", "", query.split)
        snps.current.ids[snps.na.ids[j]] <- paste0(query.split, collapse="#")
        resolved[j] <- TRUE
        
      }
      
    }

  }
  
  # Print unresolved and resolved cases.
  print(paste0("Unresolved cases (", length(which(!resolved)), "):"))
 # print(data.frame(mismatches[mismatches[, 3] == "TRUE", ], stringsAsFactors=F)[!resolved, ])
  print(paste0("Resolved cases (", length(which(resolved)), "):"))
  # print(cbind(as.data.frame(mismatches[mismatches[, 3] == "TRUE", ], stringsAsFactors=F)[resolved, ], NEW_SNP_ID_CURRENT=snps.current.ids[snps.na.ids][resolved]))

  # For each population, get a list of SNPs associated with that population.
  snps.by.population[[current.disease]] <- list()
  for (j in 1:length(populations)) {
    
    # For each study, look up what populations are associated with it. If that includes the current population, add all the SNPs from that study to the master list.
    snps.by.population[[current.disease]][[populations[j]]] <- c()
    for (k in 1:length(studies[[current.disease]][, 1])) {
      
      # See what populations are associated with the current study.
      current.populations <- unlist(strsplit(studies[[current.disease]]$Population[k], split="\\|"))
      
      # If those populationss include the current population, add all the SNPs from that study to the master list.
      if (populations[j] %in% current.populations) {
        
        # Retrieve the rsids for the current study.
        snps.study <- as.character(snps.current.ids[snps.filtered$STUDY.ACCESSION == studies[[current.disease]]$Study.accession[k]])
        for (snp in snps.study) {
          
          # Explode each SNP by the "#" delimiter and add them to the master list.
          exploded <- unlist(strsplit(snp, split="#"))
          if (length(exploded) >= 1)
            snps.by.population[[current.disease]][[populations[j]]] <- c(snps.by.population[[current.disease]][[populations[j]]], exploded)
          
        }
        
      }
      
    }
    
    # Remove NA values, add "rs" prefix to rsids, and retain only unique entries.
    if (!is.null(snps.by.population[[current.disease]][[populations[j]]])) {
      
      snps.by.population[[current.disease]][[populations[j]]] <- snps.by.population[[current.disease]][[populations[j]]][!is.na(snps.by.population[[current.disease]][[populations[j]]])]
      snps.by.population[[current.disease]][[populations[j]]] <- unique(snps.by.population[[current.disease]][[populations[j]]])
      snps.by.population[[current.disease]][[populations[j]]] <- paste0("rs", snps.by.population[[current.disease]][[populations[j]]])
      
    }
    
  }
  
  # Perform imputation for each population using HaploReg.
  imputation.results[[current.disease]] <- c()
  chunk.size <- 500
  if (length(snps.by.population[[current.disease]][["AFR"]]) > 0) {
    
    AFR <- snps.by.population[[current.disease]][["AFR"]]
    print(paste0("AFR unique SNPs: ", length(AFR)))
    write.table(AFR, file=paste0(home.dir, "/", current.disease, ".AFR.txt"), sep="\t", row.names=F, col.names=F, quote=F)
    AFR.results <- queryHaploreg(query=AFR, file=NULL, study=NULL, ldThresh=ld.threshold, ldPop="AFR", epi="vanilla", cons="siphy", genetypes="gencode",
                                 url="http://pubs.broadinstitute.org/mammals/haploreg/haploreg.php", timeout=500000, encoding="UTF-8", verbose=FALSE)
    AFR.results <- data.frame(AFR.results, stringsAsFactors=F)
    print(paste0("AFR imputation resulted in ", length(AFR.results[, 1]), " entries."))
    imputation.results[[current.disease]] <- rbind(imputation.results[[current.disease]], AFR.results)
    
  }
  if (length(snps.by.population[[current.disease]][["AMR"]]) > 0) {
    
    AMR <- snps.by.population[[current.disease]][["AMR"]]
    if (current.disease == "UD") {
      AMR <- AMR[-which(AMR == "rs748443812")]
      AMR <- AMR[-which(AMR == "rs782472239")]
    }
    print(paste0("AMR unique SNPs: ", length(AMR)))
    write.table(AMR, file=paste0(home.dir, "/", current.disease, ".AMR.txt"), sep="\t", row.names=F, col.names=F, quote=F)
    AMR.results <- queryHaploreg(query=AMR, file=NULL, study=NULL, ldThresh=ld.threshold, ldPop="AMR", epi="vanilla", cons="siphy", genetypes="gencode",
                                 url="http://pubs.broadinstitute.org/mammals/haploreg/haploreg.php", timeout=500000, encoding="UTF-8", verbose=FALSE)
    AMR.results <- data.frame(AMR.results, stringsAsFactors=F)
    print(paste0("AMR imputation resulted in ", length(AMR.results[, 1]), " entries."))
    imputation.results[[current.disease]] <- rbind(imputation.results[[current.disease]], AMR.results)
    
  }
  if (length(snps.by.population[[current.disease]][["ASN"]]) > 0) {
    
    ASN <- snps.by.population[[current.disease]][["ASN"]]
    print(paste0("ASN unique SNPs: ", length(ASN)))
    write.table(ASN, file=paste0(home.dir, "/", current.disease, ".ASN.txt"), sep="\t", row.names=F, col.names=F, quote=F)
    if (length(ASN) > chunk.size) {
      ASN.results <- c()
      for (i in 1:floor(length(ASN)/chunk.size)) {
        ASN.results.temp <- queryHaploreg(query=ASN[(1+(i-1)*chunk.size):(i*chunk.size)], file=NULL, study=NULL, ldThresh=ld.threshold, ldPop="ASN", epi="vanilla", cons="siphy", genetypes="gencode",
                                          url="http://pubs.broadinstitute.org/mammals/haploreg/haploreg.php", timeout=500000, encoding="UTF-8", verbose=FALSE)
        ASN.results <- rbind(ASN.results, data.frame(ASN.results.temp))
      }
      ASN.results.temp <- queryHaploreg(query=ASN[(1+floor(length(ASN)/chunk.size)*chunk.size):(length(ASN))], file=NULL, study=NULL, ldThresh=ld.threshold, ldPop="ASN", epi="vanilla", cons="siphy", genetypes="gencode",
                                        url="http://pubs.broadinstitute.org/mammals/haploreg/haploreg.php", timeout=500000, encoding="UTF-8", verbose=FALSE)
      ASN.results <- rbind(ASN.results, data.frame(ASN.results.temp))
    } else {
      ASN.results <- queryHaploreg(query=ASN, file=NULL, study=NULL, ldThresh=ld.threshold, ldPop="ASN", epi="vanilla", cons="siphy", genetypes="gencode",
                                   url="http://pubs.broadinstitute.org/mammals/haploreg/haploreg.php", timeout=500000, encoding="UTF-8", verbose=FALSE)
      ASN.results <- data.frame(ASN.results, stringsAsFactors=F)
    }
    print(paste0("ASN imputation resulted in ", length(ASN.results[, 1]), " entries."))
    imputation.results[[current.disease]] <- rbind(imputation.results[[current.disease]], ASN.results)
    
  }
  if (length(snps.by.population[[current.disease]][["EUR"]]) > 0) {
    
    EUR <- snps.by.population[[current.disease]][["EUR"]]
    print(paste0("EUR unique SNPs: ", length(EUR)))
    write.table(EUR, file=paste0(home.dir, "/", current.disease, ".EUR.txt"), sep="\t", row.names=F, col.names=F, quote=F)
    if (length(EUR) > chunk.size) {
      EUR.results <- c()
      for (i in 1:floor(length(EUR)/chunk.size)) {
        EUR.results.temp <- queryHaploreg(query=EUR[(1+(i-1)*chunk.size):(i*chunk.size)], file=NULL, study=NULL, ldThresh=ld.threshold, ldPop="EUR", epi="vanilla", cons="siphy", genetypes="gencode",
                                     url="http://pubs.broadinstitute.org/mammals/haploreg/haploreg.php", timeout=500000, encoding="UTF-8", verbose=FALSE)
        EUR.results <- rbind(EUR.results, data.frame(EUR.results.temp))
      }
      EUR.results.temp <- queryHaploreg(query=EUR[(1+floor(length(EUR)/chunk.size)*chunk.size):(length(EUR))], file=NULL, study=NULL, ldThresh=ld.threshold, ldPop="EUR", epi="vanilla", cons="siphy", genetypes="gencode",
                                        url="http://pubs.broadinstitute.org/mammals/haploreg/haploreg.php", timeout=500000, encoding="UTF-8", verbose=FALSE)
      EUR.results <- rbind(EUR.results, data.frame(EUR.results.temp))
    } else {
      EUR.results <- queryHaploreg(query=EUR, file=NULL, study=NULL, ldThresh=ld.threshold, ldPop="EUR", epi="vanilla", cons="siphy", genetypes="gencode",
                                   url="http://pubs.broadinstitute.org/mammals/haploreg/haploreg.php", timeout=500000, encoding="UTF-8", verbose=FALSE)
      EUR.results <- data.frame(EUR.results, stringsAsFactors=F)
    }
    print(paste0("EUR imputation resulted in ", length(EUR.results[, 1]), " entries."))
    imputation.results[[current.disease]] <- rbind(imputation.results[[current.disease]], EUR.results)
    
  }
  
  # Print total # of imputed entries.
  print(paste0("Total number of imputed entries: ", length(imputation.results[[current.disease]][, 1])))

}



# Save results.
save(studies, file=paste0(home.dir, "/studies.Rdata"))
save(associations, file=paste0(home.dir, "/associations.Rdata"))
save(snps.by.population, file=paste0(home.dir, "/snps.by.population.Rdata"))
save(imputation.results, file=paste0(home.dir, "/imputation.results.Rdata"))


# Process HaploReg results ------------------------------------------------


# Load results.
load(file=paste0(home.dir, "/imputation.results.Rdata"))

# Process imputed data.
processed.snps.noncoding <- list()
processed.snps.merged <- c()
processed.snps.all <- list()
chunk.size <- 200
for (i in 1:length(diseases)) {
  
  # Set the current disease and summarize imputation results.
  current.disease <- diseases[i]
  print(paste0("Processing SNPs for disease: ", current.disease))
  print(paste0("# of entries from imputation: ", length(imputation.results[[current.disease]][, 1])))
  print(paste0("# of unique rsids: ", length(unique(imputation.results[[current.disease]]$rsID))))
  
  # Count the # of unique SNPs w/ and w/o positional info.
  # print(table(which(imputation.results[[current.disease]][, 2] == "") %in% which(is.na(imputation.results[[current.disease]][, 1])))) # Should be TRUE.
  has.pos <- imputation.results[[current.disease]][!is.na(imputation.results[[current.disease]][, 1]), ]
  missing.pos <- imputation.results[[current.disease]][is.na(imputation.results[[current.disease]][, 1]), ]
  print(paste0("# of unique rsids with positional info: ", length(unique(has.pos$rsID))))
  print(paste0("# of unique rsids without positional info: ", length(unique(missing.pos$rsID))))
  # print(table(unique(missing.pos$rsID) %in% unique(has.pos$rsID))) # Should be FALSE.
  
  # Update rsids w/o positional info. First update any merged rsids then query NCBI's dbSNP for information.
  num.updated <- 0
  num.recovered <- 0
  original.rsids <- missing.pos$rsID
  for (j in 1:length(missing.pos$rsID)) {
    
    missing.pos$rsID[j] <- update.rsid(missing.pos$rsID[j], rs.merge)
    if (missing.pos$rsID[j] != original.rsids[j])
      num.updated <- num.updated + 1
    
  }
  if (length(missing.pos$rsID) > chunk.size) {
    
    results <- c()
    for (j in 1:floor(length(missing.pos$rsID)/chunk.size)) {
      results <- rbind(results, ncbi_snp_query(snps=missing.pos$rsID[(1 + (j - 1)*chunk.size):(j*chunk.size)]))
    }
    results <- rbind(results, ncbi_snp_query(snps=missing.pos$rsID[(1 + (ceiling(length(missing.pos$rsID)/chunk.size) - 1)*chunk.size):length(missing.pos$rsID)]))

  } else {
    
    results <- ncbi_snp_query(snps=missing.pos$rsID)
    
  }
  
  # Update the information for SNPs w/o positional info using the query results.
  unresolved.rsids <- c()
  for (j in 1:length(missing.pos[, 1])) {
    
    # There should be one unique match per entry.
    match <- unique(results[results$Marker == missing.pos$rsID[j], ])
    if (length(match[, 1]) == 1) {
    
      if (!is.na(match$Chromosome) & !is.na(match$BP)) {
        
        missing.pos$chr[j] <- match$Chromosome
        missing.pos$pos_hg38[j] <- match$BP
        num.recovered <- num.recovered + 1
        
      } else {
        
        # print(paste0(missing.pos$rsID[j], " not recovered"))
        unresolved.rsids <- c(unresolved.rsids, missing.pos$rsID[j])
        
      }
        
    } else if (length(match[, 1]) == 0) {
      
      print(paste0(missing.pos$rsID[j], " did not return any results"))
      unresolved.rsids <- c(unresolved.rsids, missing.pos$rsID[j])
      
    }
    
  }
  
  
  # Print how many rsids were updated and how many had their positional info recovered.
  print(paste0(num.updated, " of ", length(missing.pos[, 1]), " rsids updated."))
  print(paste0(num.recovered, " of ", length(missing.pos[, 1]), " rsids had their positional info recovered."))
  
  # Use the old rsids for now (we can add the new rsids back later anytime).
  mapping <- unique(cbind(original.rsids, missing.pos$rsID))
  if (!(length(original.rsids) == dim(missing.pos)[1]))
    print("Inconsistency between original and updated rsids")
  missing.pos$rsID <- original.rsids
  
  # Recombine the recovered rsids with the original rsids with positional info. Also print any rsids that were not recovered.
  combined.pos <- rbind(has.pos, missing.pos[!is.na(missing.pos$chr), ])
  print(paste0("# of unique rsids recovered: ", length(unique(missing.pos[!is.na(missing.pos$chr), ]))))
  print(paste0("# of unique rsids which could not be recovered: ", length(unique(missing.pos$rsID[is.na(missing.pos$chr)]))))
  # print(table(is.na(combined.pos[, 1]))) # Should be FALSE.
  print(paste0("Check: ", length(unique(unresolved.rsids))))
  print(unique(unresolved.rsids))
  
  # Show table of converted rsids in case one of the rsids to be recovered was updated.
  unresolved.mapping <- mapping[(mapping[, 1] %in% unresolved.rsids) | (mapping[, 2] %in% unresolved.rsids), ]
  print(table(unresolved.mapping[, 1] == unresolved.mapping[, 2]))
  print(unresolved.mapping[unresolved.mapping[, 1] != unresolved.mapping[, 2], ])

  # Create GenomicRanges object with rsids to be lifted over.
  gr <- GRanges(seqnames=as.character(combined.pos$chr),
          ranges=IRanges(as.numeric(combined.pos$pos_hg38) - 1, end=as.numeric(combined.pos$pos_hg38)),
          strand=rep("*", length(combined.pos[, 1])),
          rsid=as.character(combined.pos$rsID),
          ref=as.character(combined.pos$ref),
          alt=as.character(combined.pos$alt),
          query_snp=as.character(combined.pos$query_snp_rsid),
          is_query_snp=as.character(combined.pos$is_query_snp))


  hg38.coords <- data.frame(gr, stringsAsFactors=F)
  hg38.coords <- hg38.coords[, c(1, 2, 3, 6, 7, 8, 9, 10)]

  unique.coords <- hg38.coords %>% group_by(.dots=c("seqnames", "start", "end")) %>% 
    summarize(rsid=paste0(unique(rsid), collapse=","), ref=paste0(ref, collapse=","), alt=paste0(alt, collapse=","), query_snp=paste0(unique(query_snp), collapse=","))
  unique.coords <- data.frame(unique.coords, stringsAsFactors=F)

  num.pos <- length(unique.coords[, 1])
  for (j in 1:num.pos) {
    
    unique.coords$ref[j] <- paste0(unique(unlist(strsplit(unique.coords$ref[j], split=","))), collapse=",")
    unique.coords$alt[j] <- paste0(unique(unlist(strsplit(unique.coords$alt[j], split=","))), collapse=",")
    
  }
  

  # Print entries where two rsids have the same coordinates, and the reconciled entry under the new approach.
  unique.coords.old <- unique(hg38.coords[, 1:6])
  duplicate.coords <- unique.coords.old[duplicated(unique.coords.old[, 1:3]), 1:3]
  if (length(duplicate.coords[, 1]) > 0) {

    for (j in 1:length(duplicate.coords[, 1])) {

      duplicate.indices <- which((unique.coords.old[, 1] == duplicate.coords[j, 1]) & (unique.coords.old[, 2] == duplicate.coords[j, 2]) & (unique.coords.old[, 3] == duplicate.coords[j, 3]))
      # print(unique.coords.old[duplicate.indices, ])

    }
    
  }
  unique.coords[grep(",", unique.coords$rsid), ]
  
  # Confirm that the #s of unique rsids and unique positions agree. See above comment.
  print(length(unique(unique.coords[, 1:3])[, 1]) == length(unique(unique.coords$rsid))) # Should be TRUE.

# Add logical field determining whether entry is a query SNP.
  processed.snps <- cbind(unique.coords, is_query_snp=rep(FALSE, num.pos))
  colnames(processed.snps)[7] <- "query_snp_rsid"
  processed.snps$seqnames <- as.character(processed.snps$seqnames)
  processed.snps$seqnames <- as.character(processed.snps$seqnames)
  processed.snps$start <- as.integer(processed.snps$start)
  processed.snps$end <- as.integer(processed.snps$end)
  processed.snps$rsid <- as.character(processed.snps$rsid)
  processed.snps$ref <- as.character(processed.snps$ref)
  processed.snps$alt <- as.character(processed.snps$alt)
  processed.snps$query_snp_rsid <- as.character(processed.snps$query_snp_rsid)
  processed.snps$is_query_snp <- as.logical(processed.snps$is_query_snp)
  for (j in 1:num.pos) {
    
    rsids <- unlist(strsplit(processed.snps$rsid[j], split=","))
    query.rsids <- unlist(strsplit(processed.snps$query_snp_rsid[j], split=","))
    if (length(rsids) == 1) {
      
      if (rsids %in% query.rsids)
        processed.snps$is_query_snp[j] <- TRUE
      
    } else {
      
      print(processed.snps[j, ])
      if (any(rsids %in% query.rsids)) {
        
        print("Query SNP detected.")
        processed.snps$is_query_snp[j] <- TRUE
        
      }
      
    }
    
  }
  
  # Set each SNP to have a width of at least 1.
  interval.width <- processed.snps$end - processed.snps$start
  processed.snps$start[interval.width == 0] <- processed.snps$start[interval.width == 0] - 1


  # Determine which SNPs intersect exonic intervals and which are noncoding.
  intersect <- bedTools.2in(bed1=processed.snps, bed2=exons, opt.string="-c")
  processed.snps.noncoding[[current.disease]] <- processed.snps[intersect[, 9] == 0, ]

  # Store all processed SNPs for reference later.
  processed.snps.all[[current.disease]] <- processed.snps

  # Print final #s of processed SNPs.
  print(paste0("# of unique tag SNPs: ", length(unique(processed.snps[processed.snps[, 8] == TRUE, 4]))))
  print(paste0("# of unique imputed SNPs: ", length(unique(processed.snps[, 4]))))
  print(paste0("# of unique noncoding imputed SNPs: ", length(processed.snps.noncoding[[current.disease]][, 1])))

  # Add "exon_overlap" column as well as any redundant SNPs (rsids at the same position).
  intersect <- bedTools.2in(bed1=processed.snps.all[[current.disease]][, 1:8], bed2=exons, opt.string="-c")
  processed.snps.all[[current.disease]] <- cbind(processed.snps.all[[current.disease]][, 1:8], intersect[, 9] > 0)
  
  # Write processed SNPs to files for downstream analysis.
  write.table(processed.snps.all[[current.disease]], file=paste0(home.dir, "/features/", current.disease, ".complete.snps.features.bed"),
              sep="\t", row.names=F, col.names=F, quote=F)
  # write.table(processed.snps.noncoding[[current.disease]], file=paste0(home.dir, "features/", current.disease, ".noncoding.snps.features.bed"),
  #             sep="\t", row.names=F, col.names=F, quote=F)
  write.table(sortBed(unique(processed.snps.all[[current.disease]][, 1:4])), file=paste0(home.dir, "/", current.disease, ".complete.snps.display.bed"),
              sep="\t", row.names=F, col.names=F, quote=F)
  # write.table(sortBed(unique(processed.snps.noncoding[[current.disease]][, 1:4])), file=paste0(home.dir, "/", current.disease, ".noncoding.snps.display.bed"),
  #             sep="\t", row.names=F, col.names=F, quote=F)
  
  # Final check to make sure all query SNPs are also entries.
  print(table(hg38.coords$query_snp %in% hg38.coords$rsid))
  
}

# Save results.
save(processed.snps.all, file=paste0(home.dir, "/processed.snps.all.Rdata"))


# End ---------------------------------------------------------------------

