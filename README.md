# Korbous (Aïn Kanassira) coastal hydrothermal system: 16S rRNA gene amplicon data, metadata and analysis scripts

Data and code supporting:

Gannoun H, Ben Gaida L, Chebbi H, Dhifalli K, Fersi C, Valette C, Casalot L, Garel M, Pillot G, Armougom F, Liebgott P-P. *Chemical and microbial signatures of seawater input and potential subsurface connectivity across a coastal hydrothermal land–sea interface (Korbous, Tunisia).* Submitted to *Microbial Ecology*.

Corresponding authors: Pierre-Pol Liebgott (pierre-pol.liebgott@ird.fr), Fabrice Armougom (fabrice.armougom@mio.osupytheas.fr)

Raw sequencing reads: European Nucleotide Archive, BioProject PRJEB129037

## Study overview

Fifteen composite samples (water and sediment) were collected on 18 June 2021 at the Aïn Kanassira hydrothermal outflow system (Korbous, Cap Bon, north-eastern Tunisia): one terrestrial spring (T), two submarine venting zones (MA, MB) and a marine reference site. The V4 region of the 16S rRNA gene was amplified with primers 515F-Y/806R-B and sequenced on an Illumina MiSeq (2 × 250 bp). Amplicon sequence variants (ASVs) were inferred with DADA2 and classified against SILVA 138.2.

| Sample | Compartment | Site | Material |
|---|---|---|---|
| MWC | Marine reference | Reference | Seawater (chemical reference only, see below) |
| MSed | Marine reference | Reference | Surface sediment |
| TWS, MAWS, MBWS | Hydrothermal water | T, MA, MB | Discharge water |
| MASed, MA10, MA40, MBSed, MB10, MB40 | Submarine hydrothermal sediment | MA, MB | Surface, ~10 cm and ~40 cm |
| TSed1, TSed2, T10, T40 | Terrestrial hydrothermal sediment | T | Outflow, spring-fed basin, ~10 cm and ~40 cm |

**MWC.** The seawater sample could only be amplified with a two-step PCR; 73% of its reads belong to genera commonly reported as reagent contaminants and it contains no SAR11 read. It is therefore excluded from all community analyses (alpha diversity, PCA, RDA, heatmap, UpSet, Bray–Curtis, Fig. 5A and C) and used only as the chemical seawater reference (ion ratios of Fig. 5B). Its counts are kept in the data files, and `TableS7_robustness.R` repeats the PCA and RDA with MWC included. Community analyses therefore use 14 samples.

## Contents

```
README.md
run_all.R                      runs every analysis script in order (not 01_dada2_pipeline.R)
00_common.R                    shared objects: data import (MWC excluded), palettes, CLR transformation, contaminant genera
01_dada2_pipeline.R            sequence processing: DADA2, SILVA classification, filtering, tree, phyloseq object
Fig2_alpha_diversity.R         Fig. 2, Fig. S1, Table S1
Fig3_PCA_heatmap.R             Fig. 3
Fig4_RDA.R                     Fig. 4, Tables S3 to S6 (Table S5 includes the global test of the candidate model)
Fig5_marine_influence.R        Fig. 5, Tables S8 (Bray–Curtis similarity) and S9
FigS2_UpSet.R                  Fig. S2, Table S2
TableS7_robustness.R           Tables S7 and S7b
sessionInfo.txt                R and package versions used for the analyses
data/
  seqtab_nochim.rds            chimera-free DADA2 sequence table (15 libraries), input of 01_dada2_pipeline.R
  mapfile.txt                  sample sheet: FASTQ prefix, barcode, compartment, site, material, depth, chemistry
  asv_counts_raw_corrected.csv
  asv_counts_rarefied_corrected.csv
  taxonomy.csv
  asv_sequences.fasta
  asv_tree.nwk
  metadata.csv
  TableS1_alpha_diversity_Fabrice.csv
  Korbous_ENA_metadata.xlsx    sample and run metadata prepared for the ENA submission
```

## Sequence processing (01_dada2_pipeline.R)

Reads (2 × 250 bp) were processed with DADA2 following the CNRS ANF MEtaBioDiv workflow. `filterAndTrim` removed the first 20 nucleotides of each read (primer region) and the last 20 nucleotides, and discarded reads with ambiguous bases, more than two expected errors or fewer than 150 nucleotides (`trimLeft = 20, trimRight = 20, minLen = 150, maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE`). No cutadapt step was used. Error models were learned separately for forward and reverse reads, reads were denoised sample by sample, pairs were merged with `maxMismatch = 1`, no length filter was applied, and chimeras were removed with the consensus method.

Taxonomy: SILVA SSU 138.2 species-level training set (`silva_nr99_v138.2_toSpecies_trainset.fa.gz`, https://doi.org/10.5281/zenodo.14169026), `assignTaxonomy(minBoot = 60)`, seven ranks. The Species column comes from this training set; `addSpecies` was not used.

ASV identifiers: ASVn is the n-th column of `seqtab_nochim.rds` (decreasing total abundance); identifiers are kept through all filtering steps.

Filtering, in this order:
1. ASVs affiliated with *Cutibacterium* removed (ASV25, ASV303, ASV468, ASV2276);
2. ASVs never reaching 0.01% of the reads of a sample removed (relative abundances computed after step 1);
3. ASV29 (*Staphylococcus*) and ASV30 (*Lactobacillus*), likely of human origin, removed.

Unassigned, chloroplast and mitochondrial sequences are retained (< 1% of reads in every sample). Result: 3,559 ASVs in the 15 libraries (3,436 in the 14 samples used for community analyses).

Tree: DECIPHER alignment, Jukes–Cantor distances (phangorn), neighbour-joining tree, midpoint rooting.

By default (`start_from_seqtab = TRUE`) the script starts from `data/seqtab_nochim.rds` and rebuilds the article tables exactly (same 3,559 ASVs and counts); with `FALSE` it re-runs DADA2 from the FASTQ files, and counts may differ by a few reads with other DADA2 versions. Outputs go to `data_rebuilt/` and are compared automatically with `data/`. The taxonomy step needs about 16 GB of RAM. The SILVA training set must be placed in `ref/`.

Rarefaction (`rarefy_even_depth`) cannot be reproduced exactly across R and phyloseq versions: the rarefied table used in the article is `data/asv_counts_rarefied_corrected.csv`.

## Data files

**asv_counts_raw_corrected.csv.** ASV count table after the filtering described above (3,559 ASVs × 15 libraries; rows = ASVs, columns = samples). Not rarefied. Used for the PCA, RDA and rarefaction curves.

**asv_counts_rarefied_corrected.csv.** Same table rarefied to 91,631 reads per sample, the lowest sequencing depth (T10). Used for alpha diversity, the heatmap, the UpSet analysis and relative abundances.

**taxonomy.csv.** SILVA SSU 138.2 classification of each ASV (naïve Bayesian classifier, species-level training set, minBoot = 60): Kingdom to Species. Species assignments are not used as identifications in the article.

**asv_sequences.fasta.** Representative sequence of each ASV, with identifiers matching the count and taxonomy tables.

**asv_tree.nwk.** Neighbour-joining tree of the ASV sequences, used to calculate Faith's PD.

**metadata.csv.** One row per sample.

| Column | Description | Unit |
|---|---|---|
| Sample | Sample code | |
| Zone | Compartment | |
| Site | Sampling site (Reference, T, MA, MB) | |
| Temperature | In situ temperature | °C |
| pH | pH measured immediately after collection | |
| Calcium, Magnesium, Potassium, Sodium | Atomic absorption spectrometry | g L⁻¹ |
| Chloride | Ion chromatography | g L⁻¹ |
| Sulfate | Spectrophotometry (Hach Lange DR 1900) | g L⁻¹ |
| Sulfur | Dissolved sulfide (colloidal copper sulfide method) | mM |

For sediments, ion and sulfide concentrations refer to interstitial water recovered by centrifugation. Dissolved sulfide was not detected in MWC and MSed; these values were replaced by 0.012 mM (MSed) and 0.001 mM (MWC, used only in the robustness analysis that includes it).

**TableS1_alpha_diversity_Fabrice.csv.** Alpha-diversity indices of the 15 libraries computed on the rarefied table (observed ASVs, Shannon, Gini–Simpson, Faith's PD). PD values are read from this file by `Fig2_alpha_diversity.R`, which checks the other indices against those recomputed from the rarefied table. The PD_ses columns of this file are not used in the article.

## Reproducing the analyses

Set the working directory to the folder containing `00_common.R` and run:

```r
source("run_all.R")
```

Figures are written to `figures/` and tables to `tables/`. Each script can also be run on its own. Required R packages: vegan, zCompositions, compositions, ggplot2, patchwork, ggrepel, dplyr, tidyr, ComplexHeatmap, circlize, scales. Exact versions are listed in `sessionInfo.txt`; `run_all.R` writes the session information of each run to `tables/sessionInfo.txt`.

Checks reproduced by the scripts (14 samples): CLR matrix of 14 samples × 1,118 ASVs; PC1 47.47% and PC2 11.52%; global test of the model with the six predictors retained after VIF screening, adjusted R² 43.7%, p = 0.0004; forward selection retains dissolved sulfide, pH and temperature (adjusted R² 26.1%, 36.0% and 41.0%); RDA axes 40.44% and 9.32% (RDA1 p < 0.001, RDA2 p = 0.04); UpSet (≥ 10 reads): 1,855 ASVs found only in terrestrial sediments, 30 ASVs shared by the four compartments (34% of reads), 118 ASVs below 10 reads in every sample. With MWC included (Table S7): 15 × 1,171 ASVs, PC1 44.34% and PC2 10.99%, sulfide + pH retained (adjusted R² 31.7%).

## License

Data: Creative Commons Attribution 4.0 International (CC BY 4.0). Code: MIT License.

## Funding

IRD and MESRS through the Laboratoire Mixte International LMI BIOTEC VAL; A*Midex Foundation (Aix-Marseille University), project BIOBLUE (AMX-22-IN1-23); Provence-Alpes-Côte d'Azur Region, project GRILLON.
