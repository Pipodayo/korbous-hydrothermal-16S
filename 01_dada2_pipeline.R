############################################################
## Korbous (Aïn Kanassira) - 01_dada2_pipeline.R
## From raw FASTQ files to the ASV tables, taxonomy, sequences, tree and phyloseq object
## used in the article (Gannoun et al., submitted to Microbial Ecology)
##
## Part A (sections 1-4): FASTQ -> chimera-free ASV table (DADA2)
## Part B (sections 5-10): ASV table -> taxonomy, tree, phyloseq object, cleaned tables
##
## The chimera-free table used in the article is distributed as data/seqtab_nochim.rds.
## With start_from_seqtab = TRUE (default) Part A is skipped and Part B reproduces
## data/asv_counts_raw_corrected.csv exactly (3,559 ASVs, identical counts).
## Re-running Part A from the FASTQ files gives counts that may differ by a few reads,
## depending on the DADA2 version.
##
## Run from the folder containing this script. Expected layout:
##   Rawdata/            30 FASTQ files (15 samples x R1/R2), as deposited in the ENA (Part A only)
##   ref/                silva_nr99_v138.2_toSpecies_trainset.fa.gz (https://zenodo.org/records/14169026)
##   data/mapfile.txt    sample metadata (one row per sample, column SampleID)
## Outputs go to data_rebuilt/; the files distributed with the article (data/) are never overwritten.
## The taxonomy step needs about 16 GB of RAM.
##
## R 4.3; dada2, phyloseq, Biostrings, DECIPHER, phangorn
############################################################
suppressPackageStartupMessages({
  library(dada2); library(phyloseq); library(Biostrings)
  library(DECIPHER); library(phangorn)
})

## ---------------- parameters ----------------
start_from_seqtab <- TRUE                        # FALSE: re-run DADA2 from the FASTQ files (Part A)
raw_dir    <- "Rawdata"
ref_train  <- "ref/silva_nr99_v138.2_toSpecies_trainset.fa.gz"
map_file   <- "data/mapfile.txt"
out_dir    <- "data_rebuilt"
threads    <- .Platform$OS.type != "windows"   # dada2 multithreading is not available on Windows
min_rel    <- 0.01                              # % of reads: an ASV is kept if it reaches this in >= 1 sample
human_asv  <- c("ASV29", "ASV30")               # Staphylococcus and Lactobacillus ASVs removed as human-associated
dir.create(out_dir, showWarnings = FALSE)
set.seed(123)

meta <- read.delim(map_file, check.names = FALSE, stringsAsFactors = FALSE)
rownames(meta) <- meta$SampleID

if (start_from_seqtab) {
  seqtab.nochim <- readRDS("data/seqtab_nochim.rds")
} else {
  ## ---------------- 1. FASTQ files and sample names ----------------
  dir.create(file.path(raw_dir, "filtered"), showWarnings = FALSE)
  fnFs <- sort(list.files(raw_dir, pattern = "_R1\\.fastq\\.gz$", full.names = TRUE))
  fnRs <- sort(list.files(raw_dir, pattern = "_R2\\.fastq\\.gz$", full.names = TRUE))
  smp  <- sub("_.*$", "", basename(fnFs))                     # "MWC_GGTCAC-JM4Y5_L001_R1.fastq.gz" -> "MWC"
  stopifnot(length(fnFs) == 15, identical(smp, sub("_.*$", "", basename(fnRs))), setequal(smp, meta$SampleID))

  ## ---------------- 2. Quality filtering and primer removal ----------------
  ## Reads are 2 x 250 bp. The first 20 nt (primer region: 515F-Y is 19 nt, 806R-B is 20 nt)
  ## and the last 20 nt of each read are removed; no cutadapt step.
  filtFs <- file.path(raw_dir, "filtered", paste0(smp, "_F_filt.fastq.gz")); names(filtFs) <- smp
  filtRs <- file.path(raw_dir, "filtered", paste0(smp, "_R_filt.fastq.gz")); names(filtRs) <- smp
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = 20, trimRight = 20, minLen = 150,
                       maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = threads)

  ## ---------------- 3. Error models and ASV inference (sample by sample) ----------------
  errF <- learnErrors(filtFs, multithread = threads)
  errR <- learnErrors(filtRs, multithread = threads)
  dadaFs <- dada(filtFs, err = errF, multithread = threads)
  dadaRs <- dada(filtRs, err = errR, multithread = threads)

  ## ---------------- 4. Read merging, sequence table, chimera removal ----------------
  ## No length filter is applied (merged sequences of 214-377 bp are retained).
  mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, maxMismatch = 1)
  seqtab  <- makeSequenceTable(mergers)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = threads)
  rownames(seqtab.nochim) <- smp
  saveRDS(seqtab.nochim, file.path(out_dir, "seqtab_nochim.rds"))

  getN <- function(x) sum(getUniques(x))
  track <- data.frame(Sample = smp, input = out[, 1], filtered = out[, 2], denoisedF = sapply(dadaFs, getN),
                      denoisedR = sapply(dadaRs, getN), merged = sapply(mergers, getN), nonchim = rowSums(seqtab.nochim))
  write.csv(track, file.path(out_dir, "read_tracking.csv"), row.names = FALSE)
}
smp <- rownames(seqtab.nochim)
stopifnot(setequal(smp, meta$SampleID))

## ---------------- 5. ASV identifiers ----------------
## ASVn is the n-th column of the chimera-free table (decreasing total abundance).
## These identifiers are kept unchanged through all later filtering steps.
asv_seq <- colnames(seqtab.nochim)
asv_id  <- paste0("ASV", seq_along(asv_seq))
names(asv_seq) <- asv_id

## ---------------- 6. Taxonomy (SILVA 138.2) ----------------
## Naive Bayesian classifier down to species with the "toSpecies" training set, minimum bootstrap 60
taxa <- assignTaxonomy(seqtab.nochim, ref_train, minBoot = 60, multithread = threads,
                       taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
stopifnot(identical(rownames(taxa), asv_seq))
rownames(taxa) <- asv_id

## ---------------- 7. Phylogenetic tree ----------------
## DECIPHER alignment, Jukes-Cantor distances, neighbour-joining tree, midpoint rooting.
## The tree used in the article is distributed as data/asv_tree.nwk.
dna  <- DNAStringSet(asv_seq)
aln  <- AlignSeqs(dna, anchor = NA, verbose = FALSE)
dm   <- dist.ml(phyDat(as.matrix(aln), type = "DNA"), model = "JC69")
tree <- midpoint(NJ(dm))

## ---------------- 8. Phyloseq object ----------------
otu <- seqtab.nochim; colnames(otu) <- asv_id
ps_all <- phyloseq(otu_table(otu, taxa_are_rows = FALSE), tax_table(taxa),
                   sample_data(meta[smp, ]), phy_tree(tree), dna)
saveRDS(ps_all, file.path(out_dir, "phyloseq_all_ASVs.rds"))

## ---------------- 9. Cleaning ----------------
## (a) Cutibacterium (human skin) is removed (4 ASVs: ASV25, ASV303, ASV468, ASV2276)
ps <- subset_taxa(ps_all, Genus != "Cutibacterium" | is.na(Genus))
## (b) ASVs never reaching min_rel % of the reads of a sample (relative abundances computed
##     after step a) are removed
rel  <- transform_sample_counts(ps, function(x) 100 * x / sum(x))
keep <- taxa_names(ps)[apply(as(otu_table(rel), "matrix"), 2, max) >= min_rel]
## (c) one Staphylococcus and one Lactobacillus ASV, likely of human origin, are removed
stopifnot(all(tax_table(ps)[human_asv, "Genus"] == c("Staphylococcus", "Lactobacillus")))
keep <- setdiff(keep, human_asv)
ps <- prune_taxa(keep, ps)
## chloroplast, mitochondrial and unassigned sequences are retained (< 1% of reads in every sample)
saveRDS(ps, file.path(out_dir, "phyloseq_article.rds"))

## (d) rarefaction to the smallest library. The random draw differs between R/phyloseq versions;
##     the rarefied table used in the article is data/asv_counts_rarefied_corrected.csv.
ps_rar <- rarefy_even_depth(ps, rngseed = TRUE, replace = TRUE, trimOTUs = FALSE, verbose = FALSE)

## ---------------- 10. Exports (same format as data/) ----------------
smp_order <- c("MWC","MSed","TWS","MAWS","MBWS","MASed","MA10","MA40","MBSed","MB10","MB40","TSed1","TSed2","T10","T40")
tab <- function(x) t(as(otu_table(x), "matrix"))[, smp_order]
write.csv(tab(ps), file.path(out_dir, "asv_counts_raw_corrected.csv"))
write.csv(tab(ps_rar), file.path(out_dir, "asv_counts_rarefied_rerun.csv"))
write.csv(as(tax_table(ps), "matrix"), file.path(out_dir, "taxonomy.csv"))
writeXStringSet(refseq(ps), file.path(out_dir, "asv_sequences.fasta"), width = 1000)
write.tree(phy_tree(ps), file.path(out_dir, "asv_tree.nwk"))
cat("Retained:", ntaxa(ps), "ASVs x", nsamples(ps), "samples\n")

## comparison with the files distributed with the article
if (file.exists("data/asv_counts_raw_corrected.csv")) {
  ref <- as.matrix(read.csv("data/asv_counts_raw_corrected.csv", row.names = 1, check.names = FALSE))
  new <- tab(ps)
  cat("Same ASV set as data/:", setequal(rownames(ref), rownames(new)), "\n")
  common <- intersect(rownames(ref), rownames(new))
  cat("Identical counts for", sum(rowSums(ref[common, colnames(new)] != new[common, ]) == 0), "of", length(common), "shared ASVs\n")
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_dada2.txt"))
