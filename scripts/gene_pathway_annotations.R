## Gene annotations and functional enrichment analysis of Triticum aestivum miRNA - Human gene target interactions

## *************************************
#R version 4.3.2

## Settings ----------------------------
set.seed(1)

rm(list=ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(clusterProfiler, org.Hs.eg.db, AnnotationDbi, dplyr, ggplot2, ggbreak, tidyr, ReactomePA,rrvgo, readr)  


## Gene annotations ----------------------------

data <- read.delim("data/mirnas_targets_original.tsv",  sep = "\t", header = TRUE)
fname = "mirna_annotated_targets.tsv"

if(!file.exists(fname)){
  split_target <- strsplit(as.character(data$target), "\\|")
  data$transcript_id <- sapply(split_target, `[`, 1)
  data$gene_symbol <- sapply(split_target, `[`, 2)
  symbols <- data$gene_symbol
  annotations <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys = symbols,
                                     columns = c("ENSEMBL", "GENENAME"),
                                     keytype = "SYMBOL")
                                     #"ENTREZID", 

### keep only the first match 
  annotations_unique <- annotations[!duplicated(annotations$SYMBOL), ]
  data_annotated <- merge(data, annotations_unique, by.x = "gene_symbol", by.y = "SYMBOL", all.x = TRUE)
    # full table with annotated genes
  write.table(data_annotated, file = "tables/mirna_annotated_targets.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  }else{
    data_annotated <- read.delim(fname)
  } 


#get the total number of mirna by mfe threshold
data_annotated$mfe <- as.numeric(gsub(",", ".", data_annotated$mfe))

## Summary data analysis by MFE thresholds ----------------------------

thresholds <- c(-25, -30, -35)

count_mirnas <- function(threshold) {
  subset <- data_annotated[data_annotated$mfe <= threshold, ]
  miRNA = length(unique(subset$miRNA)) #conta unique! -> correct!
  gene_symbol = length(unique(subset$gene_symbol)) 
  no_rows = dim(subset)[1]
  cat("Threshold:", threshold, " - miRNA no:", miRNA, " - gene_symbol no:", gene_symbol, " - row_no:", no_rows, "\n")
  
  data.frame(
    threshold   = threshold,
    miRNA       = miRNA,
    gene_symbol = gene_symbol,
    no_rows     = no_rows,
    stringsAsFactors = FALSE
  )
}

counts <- do.call(rbind, lapply(thresholds, count_mirnas))

write.table(counts, file = "tables/mirna_counts_by_mfe_cutoff.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


### Unique targets (genes) per miRNA by thresholds
#### PLOT histogramm
df25 <- read.delim("data/tables/targets_per_miRNA_mfe-25.tsv")
df30 <- read.delim("data/tables/targets_per_miRNA_mfe-30.tsv")
df35 <- read.delim("data/tables/targets_per_miRNA_mfe-35.tsv")

df25$threshold <- "MFE ≤ -25"
df30$threshold <- "MFE ≤ -30"
df35$threshold <- "MFE ≤ -35"

df_all <- bind_rows(df25, df30, df35)


library(ggplot2)

p <- ggplot(df_all, aes(x = n_targets)) +
  geom_histogram(bins = 30, color = "black", fill = "slategray3", alpha = 0.8) +
  facet_wrap(~threshold, scales = "fixed") + #scales = "free_y" -> same scale for comparison
  labs(
    title = "",
    x = "Number of targets",
    y = "Frequency"
  ) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 17),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 11)
  )
#, face = "bold"

ggsave("Number_of_targets_per_miRNA_pub_1.png", plot = p, width = 18, height = 13, dpi = 600, units = "cm")


thresholds <- c(-25, -30, -35)

for (thr in thresholds) {
  subset <- data_annotated[data_annotated$mfe <= thr, ]
  targets_per_miRNA <- subset %>%
    group_by(miRNA) %>%
    summarise(
      n_targets = n_distinct(gene_symbol),
      target_genes = paste(sort(unique(gene_symbol)), collapse = ";"), 
      .groups = "drop") %>%
    arrange(desc(n_targets))
  out_file <- paste0("data/tables/targets_per_miRNA_with_genes_mfe", thr, ".tsv")
  write.table(targets_per_miRNA, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
}


### Gene lists by thresholds
for (thr in thresholds) {
  genes <- sort(unique(
    data_annotated$gene_symbol[data_annotated$mfe <= thr]
  ))
  
  out_file <- sprintf("tables/gene_list_trheshold%s.tsv", thr)   # ➜ genes_trheshold-25.tsv …
  write.table(
    data.frame(gene_symbol = genes, stringsAsFactors = FALSE),
    file      = out_file,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE          # omit header; set to TRUE if you want “gene_symbol”
  )
  
  message("Saved ", length(genes), " genes → ", out_file)
}


### mirna lists by thresholds
for (thr in thresholds) {
  mirnas <- sort(unique(
    data_annotated$miRNA[data_annotated$mfe <= thr]
  ))
  
  out_file <- sprintf("tables/miRNA_list_trheshold%s.tsv", thr)   # ➜ genes_trheshold-25.tsv …
  write.table(
    data.frame(miRNA = mirnas, stringsAsFactors = FALSE),
    file      = out_file,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE         
  )
  
  message("Saved ", length(mirnas), " mirnas → ", out_file)
}


### Annotated data subsets by thresholds
for (thr in thresholds) {
  test <- subset(data_annotated, data_annotated$mfe <= thr)
  
  write.table(test,
              file = paste0("tables/mirnas_targets_mfe_", abs(thr), "_full_annotated.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}


## GO Enrichment by mfe Thresholds ----------------------------

enrich_by_threshold <- function(threshold, ontology = "BP") {
  genes <- data_annotated %>%
    filter(mfe <= threshold) %>%
    pull(gene_symbol) %>% ### su intera lista geni, non lista geni per miRNA!
    unique()
  
  # Map gene symbols to Entrez IDs
  entrez_ids <- bitr(genes,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
  
  if (nrow(entrez_ids) == 0) return(NULL)
  

  # Run GO enrichment
  enrichGO(gene = entrez_ids$ENTREZID,
           OrgDb = org.Hs.eg.db,
           ont = ontology,
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2,
           readable = TRUE)
}


ontologies <- c("BP", "CC", "MF")

# GO folder all GO terms
if (!dir.exists("GO_results")) dir.create("GO_results")

for (thresh in thresholds) {
  for (ont in ontologies) {
    
    cat("Running GO enrichment for threshold", thresh, "and ontology", ont, "\n")
    
    result <- enrich_by_threshold(thresh, ontology = ont)
    
    if (!is.null(result)) {
      filename <- paste0("GO_results/GO_enrichment_mfe", thresh, "_", ont, "_1.csv")
      write.csv(result@result, file = filename, row.names = FALSE)
    }
  }
}



## GO Enrichment REDUCTIONS (REVIGO) by mfe Thresholds ----------------------------

files <- list.files("GO_results", pattern = "^GO_enrichment_mfe-.*_(BP|CC)_significant.csv$", full.names = TRUE)
log_skipped <- file("GO_results/skipped_files.log", open = "wt")

for (file in files) {
  name_parts <- strsplit(basename(file), "_|\\.")[[1]]
  mfe <- gsub("mfe-", "", name_parts[3])
  ont <- name_parts[4]
  message("Processing: ", file)

  test <- read.delim(file, sep = ",")
  test_p <- subset(test, p.adjust < 0.05)

  if (nrow(test_p) < 2) {
    msg <- paste("Skipping", file, ": not enough significant GO terms.")
    message(msg)
    writeLines(msg, log_skipped)
    next
  }

  if (nrow(test_p) > 300) {
    message("Reduced to top 300 terms based on p.adjust.")
    test_p <- test_p[order(test_p$p.adjust), ][1:300, ]
  }

  prefix <- paste0("GO_results/revigo_mfe-", mfe, "_", ont)

  simMatrix <- tryCatch({
    calculateSimMatrix(test_p$ID, orgdb = "org.Hs.eg.db", ont = ont, method = "Rel")
  }, error = function(e) {
    msg <- paste("calculateSimMatrix failed for", file, ":", e$message)
    message(msg)
    writeLines(msg, log_skipped)
    return(NULL)
  })

  if (is.null(simMatrix) || !is.matrix(simMatrix) || nrow(simMatrix) < 2 || all(is.na(simMatrix))) {
    msg <- paste("Skipping file due to invalid or empty simMatrix (possibly no terms found for", ont, "ontology):", file)
    message(msg)
    writeLines(msg, log_skipped)
    next
  }

  reducedTerms <- tryCatch({
    scores <- setNames(-log10(as.numeric(test_p$qvalue)), test_p$ID)
    reduceSimMatrix(simMatrix, scores, threshold = 0.7, orgdb = "org.Hs.eg.db")
  }, error = function(e) {
    msg <- paste("reduceSimMatrix failed for", file, ":", e$message)
    message(msg)
    writeLines(msg, log_skipped)
    return(NULL)
  })

  if (is.null(reducedTerms)) {
    msg <- paste("Skipping file due to failed reduceSimMatrix:", file)
    message(msg)
    writeLines(msg, log_skipped)
    next
  }

  saveRDS(simMatrix, file = paste0(prefix, "_simMatrix.rds"))

  saveRDS(reducedTerms, file = paste0(prefix, "_reducedTerms.rds"))

  tryCatch({
    png(paste0(prefix, "_heatmap.png"), height = 40, width = 50, res = 600, units = "cm")
    heatmapPlot(simMatrix, reducedTerms,
                annotateParent = TRUE,
                annotationLabel = "parentTerm",
                clusterRows = FALSE,
                clusterCols = FALSE,
                fontsize = 6)
    dev.off()
  }, error = function(e) {
    msg <- paste("Heatmap failed for", file, ":", e$message)
    message(msg)
    writeLines(msg, log_skipped)
  })

  tryCatch({
    png(paste0(prefix, "_scatterplot.png"), height = 40, width = 40, res = 600, units = "cm")
    p <- scatterPlot(simMatrix, reducedTerms)
    print(p)
    dev.off()
  }, error = function(e) {
    msg <- paste("Scatterplot failed for", file, ":", e$message)
    message(msg)
    writeLines(msg, log_skipped)
  })
}

close(log_skipped)



GO25 <- readRDS("GO_results/revigo_mfe-25_BP_reducedTerms.rds")

GO25_selected <- c(
  "proteasome-mediated ubiquitin-dependent protein catabolic process",
  "regulation of T cell activation",
  "leukocyte migration",
  "cellular response to environmental stimulus",
  "wound healing",
  "small GTPase mediated signal transduction",
  "positive regulation of cell adhesion",
  "establishment or maintenance of cell polarity",
  "regulation of microtubule cytoskeleton organization",
  "response to salt",
  "viral process",
  "signal release",
  "response to nutrient levels",
  "positive regulation of protein localization",
  "response to mechanical stimulus",
  "stem cell differentiation",
  "morphogenesis of a branching structure",
  "epithelial tube morphogenesis",
  "actomyosin structure organization",
  "developmental maturation"
)

# retrieve the enriched GO terms by parentTerms after reduction, to retrieve the genes involved in all pathways associated to each parent term 
GO25_matched_terms <- GO25$term[GO25$parentTerm %in% GO25_selected]

# read significant GO terms and select TERMS
GO_enrich25 <- read.delim("GO_results/GO_enrichment_mfe-25_BP_significant.csv", sep=',')
GO_enrich25_selected <- GO_enrich25[GO_enrich25$Description %in% GO25_matched_terms,]

#### PLOT
GO_enrich25_selected <- GO_enrich25_selected[order(GO_enrich25_selected$p.adjust), ]

GO_enrich25_selected$Description <- factor(GO_enrich25_selected$Description, levels = rev(GO_enrich25_selected$Description[order(GO_enrich25_selected$p.adjust)]))
GO_enrich25_selected$GeneRatio <- sapply(strsplit(as.character(GO_enrich25_selected$GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
GO_enrich25_selected$Count <- as.numeric(as.character(GO_enrich25_selected$Count))
# no plots!!!


## mfe -30
GO30 <- readRDS("GO_results/revigo_mfe-30_BP_reducedTerms.rds")

GO30_selected <- c(
  "positive regulation of lymphocyte activation",
  "positive regulation of cell adhesion",
  "cellular response to environmental stimulus",
  "wound healing",
  "mononuclear cell differentiation",
  "small GTPase mediated signal transduction",
  "cellular response to external stimulus",
  "signal release",
  "macroautophagy",
  "cellular response to mechanical stimulus",
  "establishment or maintenance of cell polarity",
  "protein autophosphorylation",
  "viral process",
  "developmental maturation",
  "positive regulation of kinase activity",
  "peptidyl-tyrosine phosphorylation",
  "lymph vessel development",
  "biological process involved in intraspecies interaction between organisms",
  "homeostasis of number of cells",
  "regulation of carbohydrate metabolic process"
)

GO30_matched_terms <- GO30$term[GO30$parentTerm %in% GO30_selected]

GO_enrich30 <- read.delim("GO_results/GO_enrichment_mfe-30_BP_significant.csv", sep=',')
GO_enrich30_selected <- GO_enrich30[GO_enrich30$Description %in% GO30_matched_terms,]

#### PLOT
GO_enrich30_selected <- GO_enrich30_selected[order(GO_enrich30_selected$p.adjust), ]

GO_enrich30_selected$Description <- factor(GO_enrich30_selected$Description, levels = rev(GO_enrich30_selected$Description[order(GO_enrich30_selected$p.adjust)]))
GO_enrich30_selected$GeneRatio <- sapply(strsplit(as.character(GO_enrich30_selected$GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
GO_enrich30_selected$Count <- as.numeric(as.character(GO_enrich30_selected$Count))

### mfe 35
GO35 <- readRDS("GO_results/revigo_mfe-35_BP_reducedTerms.rds")
GO35_matched_terms <- GO35$term

GO_enrich35 <- read.delim("GO_results/GO_enrichment_mfe-35_BP_significant.csv", sep=',')
GO_enrich35_selected <- GO_enrich35[GO_enrich35$Description %in% GO35_matched_terms,]

#### PLOT
GO_enrich35_selected <- GO_enrich35_selected[order(GO_enrich35_selected$p.adjust), ]

GO_enrich35_selected$Description <- factor(GO_enrich35_selected$Description, levels = rev(GO_enrich35_selected$Description[order(GO_enrich35_selected$p.adjust)]))
GO_enrich35_selected$GeneRatio <- sapply(strsplit(as.character(GO_enrich35_selected$GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
GO_enrich35_selected$Count <- as.numeric(as.character(GO_enrich35_selected$Count))


#annotate GO results with miRNA info
annotate_GO_with_miRNAs_by_symbol <- function(go_df, mirna_file) {
  #miRNA annotation table
  mirna_df <- read_tsv(mirna_file, show_col_types = FALSE) %>%
    select(gene_symbol, miRNA) %>%
    filter(!is.na(gene_symbol), !is.na(miRNA)) %>%
    distinct()
  go_annotated <- go_df %>%
    separate_rows(geneID, sep = "/") %>%
    rename(gene_symbol = geneID) %>%
    left_join(mirna_df, by = "gene_symbol") %>%
    group_by(ID, Description) %>%
    summarise(
      geneID = paste(unique(gene_symbol), collapse = "/"),
      miRNA = paste(unique(na.omit(miRNA)), collapse = ";"),
      .groups = "drop"
    ) %>%
    right_join(go_df, by = c("ID", "Description", "geneID"))

  return(go_annotated)
}

# all thresholds
GO_enrich25_annotated <- annotate_GO_with_miRNAs_by_symbol(
  GO_enrich25_selected, "tables/mirnas_targets_mfe_25_full_annotated.tsv"
)

GO_enrich30_annotated <- annotate_GO_with_miRNAs_by_symbol(
  GO_enrich30_selected, "tables/mirnas_targets_mfe_30_full_annotated.tsv"
)

GO_enrich35_annotated <- annotate_GO_with_miRNAs_by_symbol(
  GO_enrich35_selected, "tables/mirnas_targets_mfe_35_full_annotated.tsv"
)

write_csv(GO_enrich25_annotated, "GO_results/GO_BP_25_annotated.csv")
write_csv(GO_enrich30_annotated, "GO_results/GO_BP_30_annotated.csv")
write_csv(GO_enrich35_annotated, "GO_results/GO_BP_35_annotated.csv")

writeLines(shared_mirnas_25, "GO_results/GO-BP_mfe_-25_shared_mirnas.txt")
writeLines(shared_mirnas_30, "GO_results/GO-BP_mfe_-30_shared_mirnas.txt")
writeLines(shared_mirnas_35, "GO_results/GO-BP_mfe_-35_shared_mirnas.txt")


## GO summary by mfe Thresholds ----------------------------

go_summary <- data.frame()
for (thresh in thresholds) {
  for (ont in ontologies) {
    
    cat("Summarizing GO enrichment for threshold", thresh, "and ontology", ont, "\n")
    
    result <- enrich_by_threshold(thresh, ontology = ont)
    
    if (!is.null(result)) {
      sig_count <- sum(result@result$p.adjust < 0.05)
    } else {
      sig_count <- 0
    }
    
    go_summary <- rbind(go_summary, data.frame(
      mfe_threshold = thresh,
      ontology = ont,
      significant_terms = sig_count
    ))
  }
}

write.csv(go_summary, "GO_results/GO_summary_counts.csv", row.names = FALSE)


#plot summary
go_summary <- read.delim("data/GO_results/GO_summary_counts.csv", sep=',')
ggplot(go_summary, aes(x = factor(mfe_threshold), y = significant_terms, fill = ontology)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Significant GO terms per MFE threshold",
       x = "MFE Threshold",
       y = "Number of significant GO terms",
       fill = "GO domain") +
  theme_minimal(base_size = 12) +
  scale_y_break(c(350, 1400), scales = 0.5)  # collapse the gap

ggsave("GO_results/GO_summary_plot.png", width = 18, height = 13, dpi = 600, units="cm")



## KEGG Enrichment by mfe Thresholds ----------------------------

enrich_kegg_by_threshold <- function(threshold) {
  genes <- data_annotated %>%
    filter(mfe <= threshold) %>% #ok, include everithing below 
    pull(gene_symbol) %>%
    unique()
  
  entrez_ids <- bitr(genes,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
  
  if (nrow(entrez_ids) == 0) return(NULL)
  
  enrichKEGG(gene = entrez_ids$ENTREZID,
             organism = 'hsa',    # 'hsa' for human
             pvalueCutoff = 0.05)
}


if (!dir.exists("KEGG_results")) dir.create("KEGG_results")

kegg_summary <- data.frame()

for (thresh in thresholds) {
  kegg_result <- enrich_kegg_by_threshold(thresh)
    if (!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
    write.csv(kegg_result@result, file = paste0("KEGG_results/KEGG_mfe_", thresh, "_14082025.csv"), row.names = FALSE)
    sig_kegg <- sum(kegg_result@result$p.adjust < 0.05) #conta i TRUE!
  } else {
    sig_kegg <- 0
  }
  
    kegg_summary <- rbind(kegg_summary, data.frame(
    mfe_threshold = thresh,
    significant_KEGG = sig_kegg
  ))
}

write.csv(kegg_summary, "KEGG_results/summary_KEGG_significant_terms.csv", row.names = FALSE)
