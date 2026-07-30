## Annotations and functional enrichment analysis of Triticum aestivum miRNA - Human gene target interactions

## *************************************
#R version 4.3.2

## Settings ----------------------------
set.seed(1)

rm(list=ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(clusterProfiler, org.Hs.eg.db, AnnotationDbi, dplyr, ggplot2, ggbreak, tidyr, ReactomePA,rrvgo, readr, ggvenn)  


# miRNAs prioritisation ----------------------------

## mfe_-30
selected_pathways_mfe30 <- c("Platelet activation, signaling and aggregation","Death Receptor Signaling","VEGFA-VEGFR2 Pathway","Signaling by VEGF","Costimulation by the CD28 family","CD28 co-stimulation","Signaling by Erythropoietin","GPVI-mediated activation cascade","Platelet homeostasis","VEGFR2 mediated vascular permeability","Glycosaminoglycan metabolism","Interleukin-3, Interleukin-5 and GM-CSF signaling","Signaling by Interleukins")
reactome30 <- read.delim("Reactome_results/Reactome_mfe_-30.csv", sep=',')
reactome30_selected <- subset(reactome30, reactome30$Description %in% selected_pathways_mfe30)
subset30 <- read.delim("tables/mirnas_targets_mfe_30_full_annotated.tsv", sep="\t")

subset30_clean <- subset30 %>%
  mutate(ENTREZID = as.character(as.integer(ENTREZID))) %>%
  select(ENTREZID, ENSEMBL, gene_symbol, miRNA) %>%
  distinct()

reactome_expanded <- reactome30_selected %>%
  separate_rows(geneID, sep = "/") %>%
  rename(ENTREZID = geneID) %>%
  mutate(ENTREZID = as.character(ENTREZID)) %>%
  left_join(subset30_clean, by = "ENTREZID")

reactome_annotated <- reactome_expanded %>%
  group_by(ID, Description) %>%
  summarise(
    ENSEMBL = paste(unique(na.omit(ENSEMBL)), collapse = ";"),
    gene_symbol = paste(unique(na.omit(gene_symbol)), collapse = ";"),
    miRNA = paste(unique(na.omit(miRNA)), collapse = ";"),
    .groups = "drop"
  )

reactome30_selected_annotated <- reactome30_selected %>%
  left_join(reactome_annotated, by = c("ID", "Description"))

write.csv(reactome30_selected_annotated, "Reactome_results/Reactome_mfe_-30_associated_mirnas.csv", row.names = FALSE)

mirna_lists <- reactome30_selected_annotated %>%
  select(Description, miRNA) %>%
  filter(!is.na(miRNA)) %>%
  mutate(miRNA = strsplit(miRNA, ";"))

shared_mirnas <- Reduce(intersect, mirna_lists$miRNA)

print(shared_mirnas)

write.csv(shared_mirnas, "Reactome_results/reactome_-30_shared_mirnas.csv", row.names = FALSE)


## mfe_-25
selected_pathways_mfe25 <- c("RHO GTPase cycle","Signaling by VEGF","VEGFA-VEGFR2 Pathway","Death Receptor Signaling","Platelet activation, signaling and aggregation","Toll Like Receptor 4 (TLR4) Cascade","Toll Like Receptor TLR1:TLR2 Cascade","Toll Like Receptor 2 (TLR2) Cascade","Toll Like Receptor TLR6:TLR2 Cascade","Toll Like Receptor 10 (TLR10) Cascade","Toll Like Receptor 5 (TLR5) Cascade","Cell-Cell communication","Glycosaminoglycan metabolism","Costimulation by the CD28 family","TRAF6 mediated induction of NFkB and MAP kinases upon TLR7/8 or 9 activation","Class I MHC mediated antigen processing & presentation","Heparan sulfate/heparin (HS-GAG) metabolism","Toll Like Receptor 7/8 (TLR7/8) Cascade","MyD88 dependent cascade initiated on endosome","Antigen processing: Ubiquitination & Proteasome degradation","Neutrophil degranulation","Basigin interactions","Toll Like Receptor 3 (TLR3) Cascade")
reactome25 <- read.delim("Reactome_results/Reactome_mfe_-25.csv", sep=',')
reactome25_selected <- subset(reactome25, reactome25$Description %in% selected_pathways_mfe25)

subset25 <- read.delim("tables/mirnas_targets_mfe_25_full_annotated.tsv", sep="\t")

subset25_clean <- subset25 %>%
  mutate(ENTREZID = as.character(as.integer(ENTREZID))) %>%
  select(ENTREZID, ENSEMBL, gene_symbol, miRNA) %>%
  distinct()

reactome_expanded <- reactome25_selected %>%
  separate_rows(geneID, sep = "/") %>%
  rename(ENTREZID = geneID) %>%
  mutate(ENTREZID = as.character(ENTREZID)) %>%
  left_join(subset25_clean, by = "ENTREZID")

reactome_annotated <- reactome_expanded %>%
  group_by(ID, Description) %>%
  summarise(
    ENSEMBL = paste(unique(na.omit(ENSEMBL)), collapse = ";"),
    gene_symbol = paste(unique(na.omit(gene_symbol)), collapse = ";"),
    miRNA = paste(unique(na.omit(miRNA)), collapse = ";"),
    .groups = "drop"
  )

reactome25_selected_annotated <- reactome25_selected %>%
  left_join(reactome_annotated, by = c("ID", "Description"))

write.csv(reactome25_selected_annotated, "Reactome_results/Reactome_mfe_-25_associated_mirnas.csv", row.names = FALSE)


mirna_lists <- reactome25_selected_annotated %>%
  select(Description, miRNA) %>%
  filter(!is.na(miRNA)) %>%
  mutate(miRNA = strsplit(miRNA, ";"))

shared_mirnas <- Reduce(intersect, mirna_lists$miRNA)

print(shared_mirnas)

write.csv(shared_mirnas, "Reactome_results/reactome_-25_shared_mirnas.csv", row.names = FALSE)


## mirnas from KEGG enriched pathways ----------------------------
# MFE -25

set.seed(1)
rm(list=ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(clusterProfiler, org.Hs.eg.db, AnnotationDbi, dplyr, ggplot2, ggbreak, tidyr, ReactomePA,rrvgo)  


annotate_and_find_shared_mirnas <- function(kegg_file, mirna_file) {
  pacman::p_load(dplyr, tidyr, readr, stringr)

  kegg_df <- read_csv(kegg_file, show_col_types = FALSE)
  mirna_df <- read_tsv(mirna_file, show_col_types = FALSE)

  selected_kegg_mfe30 <- c("Focal adhesion","Signaling pathways regulating pluripotency of stem cells",
                           "Cell adhesion molecules","Cytokine-cytokine receptor interaction",
                           "Cushing syndrome","Type II diabetes mellitus","Insulin resistance",
                           "AGE-RAGE signaling pathway in diabetic complications","Gastric acid secretion",
                           "C-type lectin receptor signaling pathway","Fc gamma R-mediated phagocytosis",
                           "T cell receptor signaling pathway","Chemokine signaling pathway",
                           "Th17 cell differentiation","Th1 and Th2 cell differentiation","Platelet activation",
                           "Fc epsilon RI signaling pathway","B cell receptor signaling pathway",
                           "Leukocyte transendothelial migration","Inflammatory mediator regulation of TRP channels")

  selected_kegg_mfe25 <- c("Cell adhesion molecules","ECM-receptor interaction",
                           "Cytokine-cytokine receptor interaction","Neuroactive ligand-receptor interaction",
                           "Viral protein interaction with cytokine and cytokine receptor",
                           "Platelet activation","T cell receptor signaling pathway",
                           "Th1 and Th2 cell differentiation","Th17 cell differentiation",
                           "C-type lectin receptor signaling pathway","Fc epsilon RI signaling pathway",
                           "Fc gamma R-mediated phagocytosis","Leukocyte transendothelial migration",
                           "Leukocyte transendothelial migration - Homo sapiens (human)","Rheumatoid arthritis",
                           "Primary immunodeficiency","Gastric acid secretion","Protein digestion and absorption",
                           "Cholesterol metabolism","IL-17 signaling pathway","Hematopoietic cell lineage",
                           "Intestinal immune network for IgA production",
                           "Inflammatory mediator regulation of TRP channels","Cushing syndrome",
                           "AGE-RAGE signaling pathway in diabetic complications","Insulin resistance",
                           "Type II diabetes mellitus","Alcoholic liver disease")

  if (grepl("-25", kegg_file)) {
    kegg_df <- kegg_df %>% filter(Description %in% selected_kegg_mfe25)
  } else if (grepl("-30", kegg_file)) {
    kegg_df <- kegg_df %>% filter(Description %in% selected_kegg_mfe30)
  }

  mirna_clean <- mirna_df %>%
    mutate(ENTREZID = as.character(as.integer(ENTREZID))) %>%
    select(ENTREZID, ENSEMBL, gene_symbol, miRNA) %>%
    distinct()

  kegg_expanded <- kegg_df %>%
    separate_rows(geneID, sep = "/") %>%
    rename(ENTREZID = geneID) %>%
    mutate(ENTREZID = as.character(ENTREZID)) %>%
    left_join(mirna_clean, by = "ENTREZID")

  kegg_annotated <- kegg_expanded %>%
    group_by(ID, Description) %>%
    summarise(
      ENSEMBL = paste(unique(na.omit(ENSEMBL)), collapse = ";"),
      gene_symbol = paste(unique(na.omit(gene_symbol)), collapse = ";"),
      miRNA = paste(unique(na.omit(miRNA)), collapse = ";"),
      .groups = "drop"
    )

  kegg_final <- kegg_df %>%
    left_join(kegg_annotated, by = c("ID", "Description"))

  # shared miRNAs
  mirna_lists <- kegg_final %>%
    filter(!is.na(miRNA)) %>%
    mutate(miRNA = strsplit(miRNA, ";"))

  shared_mirnas <- Reduce(intersect, mirna_lists$miRNA)
  all_mirnas <- unique(unlist(mirna_lists$miRNA))

  base_name <- str_replace(kegg_file, ".csv", "")
  write_csv(kegg_final, paste0(base_name, "_annotated.csv"))
  writeLines(shared_mirnas, paste0(base_name, "_shared_mirnas.txt"))

  if (length(shared_mirnas) == 0 && length(all_mirnas) > 0) {
    writeLines(all_mirnas, paste0(base_name, "_all_mirnas.txt"))
  }

  return(shared_mirnas)
}

# Run for -25
shared_25 <- annotate_and_find_shared_mirnas("KEGG_results/KEGG_mfe_-25.csv", "tables/mirnas_targets_mfe_25_full_annotated.tsv")

# Run for -30
shared_30 <- annotate_and_find_shared_mirnas("KEGG_results/KEGG_mfe_-30.csv", "tables/mirnas_targets_mfe_30_full_annotated.tsv")

# Run for -35
shared_35 <- annotate_and_find_shared_mirnas("KEGG_results/KEGG_mfe_-35.csv", "tables/mirnas_targets_mfe_35_full_annotated.tsv")


### VENN diagramm mirnas

react25_mirnas <- read.delim("data/Reactome_results/reactome_-25_shared_mirnas.csv", sep=',', header=F)
react25_mirnas <- react25_mirnas$V1[react25_mirnas$V1 != "x"]
react30_mirnas <- read.delim("data/Reactome_results/reactome_-30_shared_mirnas.csv", sep=',', header=F)
react30_mirnas <- react30_mirnas$V1[react30_mirnas$V1 != "x"]

kegg30_mirnas <- read.delim("data/KEGG_results/KEGG_mfe_-30_shared_mirnas.txt", sep='\t')
kegg30_mirnas  <- kegg30_mirnas[[1]] 
kegg25_mirnas <- read.delim("KEGG_results/KEGG_mfe_-25_shared_mirnas.txt", sep='\t')
kegg25_mirnas  <- kegg25_mirnas[[1]] 

go_bp_30_mirnas <- read.delim("data/GO_results/GO-BP_mfe_-30_shared_mirnas.txt", sep='\t')
go_bp_30_mirnas  <- go_bp_30_mirnas[[1]] 
go_bp_25_mirnas <- read.delim("data/GO_results/GO-BP_mfe_-25_shared_mirnas.txt", sep='\t')
go_bp_25_mirnas  <- go_bp_25_mirnas[[1]]

### -30

venn_plot <- ggvenn(
  list(
    Reactome = react25_mirnas,
    KEGG = kegg25_mirnas,
    GO_BP = go_bp_25_mirnas
  ),
  fill_color = c("#3498DB", "Gold2", "red")
)


venn_plot <- ggvenn(
  list(
    Reactome = react25_mirnas,
    KEGG = kegg25_mirnas,
    GO_BP = go_bp_25_mirnas
  ),
  fill_color = c("#3498DB", "Gold2", "red")
) +
  theme_void() +  # removes axes and background
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 20, hjust = 0.5)) +
  ggtitle("mfe = -25")


ggplot2::ggsave(
  filename = "Venn_shared_miRNAs_Reactome_KEGG_GO-BP_-25.png",
  plot = venn_plot,
  dpi = 300,
  width = 12,
  height = 9,
  units = "in"
)


venn_plot <- ggvenn(
  list(
    Reactome = react30_mirnas,
    KEGG = kegg30_mirnas,
    GO_BP = go_bp_30_mirnas
  ),
  fill_color = c("#3498DB", "Gold2", "red")) +
  theme_void() + 
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 20, hjust = 0.5)) +
  ggtitle("mfe = -30")

ggplot2::ggsave(
  filename = "Venn_shared_miRNAs_Reactome_KEGG_GO-BP_-30.png",
  plot = venn_plot,
  dpi = 300,
  width = 12,
  height = 9,
  units = "in"
)

length(react25_mirnas)
length(kegg25_mirnas)
length(go_bp_25_mirnas)

shared_all_25 <- Reduce(intersect, list(react25_mirnas, kegg25_mirnas, go_bp_25_mirnas))
writeLines(shared_all_25, "tables/shared_miRNAs_Reactome_KEGG_GO-BP_-25.txt")
shared_all_30 <- Reduce(intersect, list(react30_mirnas, kegg30_mirnas, go_bp_30_mirnas))
writeLines(shared_all_30, "tables/shared_miRNAs_Reactome_KEGG_GO-BP_-30.txt")



#### Check target number distribution  ----------------------------

set.seed(1)

rm(list=ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, dplyr, openxlsx)  

setwd("/Users/mariangelasantorsola/Documents/Documents - santorsola-macbook/anca_macovei")


files <- c(
  "data/tables/SupplementaryFile3_targets_per_miRNA_with_genes_mfe-35.tsv",
  "data/tables/SupplementaryFile2_targets_per_miRNA_with_genes_mfe-30.tsv",
  "data/tables/SupplementaryFile1_targets_per_miRNA_with_genes_mfe-25.tsv"
)


for (file in files) {
  
  # get MFE threshold from filename
  mfe <- gsub(".*mfe-|\\.tsv", "", file)
  
    data <- read.table(
    file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  data_summary <- data %>%
    mutate(percent_targets = n_targets / sum(n_targets) * 100) %>%
    arrange(desc(n_targets))

wb <- createWorkbook()

# Complete summary table
addWorksheet(wb, "All_miRNAs")
writeData(wb, "All_miRNAs", data_summary)

# Top 20 miRNAs
addWorksheet(wb, "Top20_miRNAs")
writeData(wb, "Top20_miRNAs", head(data_summary, 20))

# dominance summary
dominance <- data_summary %>%
  summarise(
    total_targets = sum(n_targets),
    top1_percent = sum(percent_targets[1]),
    top2_percent = sum(percent_targets[1:2]),
    top5_percent = sum(percent_targets[1:5])
  )

addWorksheet(wb, "Dominance")
writeData(wb, "Dominance", dominance)

saveWorkbook(
  wb,
  file = paste0("target_frequency_summary_mfe-", mfe, ".xlsx"),
  overwrite = TRUE
)

write.table(
  data_summary,
  file = paste0("target_frequency_summary_mfe-", mfe, ".tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  head(data_summary, 20),
  file = paste0("top20_target_frequency_mfe-", mfe, ".tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
  print(paste("MFE threshold:", mfe))
  print(head(data_summary[, c("miRNA", "n_targets", "percent_targets")], 10))
  
  
### Distribution of target numbers per miRNA -----------------------------
  
  p1 <- ggplot(data_summary, aes(x = n_targets)) +
    geom_histogram(bins = 50) +
    labs(
      x = "Number of predicted target genes per tae-miRNA",
      y = "Number of tae-miRNAs",
      title = paste0(
        "Distribution of predicted targets per tae-miRNA (MFE ≤ -",
        mfe,
        " kcal/mol)"
      )
    ) +
    theme_classic()
  
  
  ggsave(
    filename = paste0("target_distribution_mfe-", mfe, ".png"),
    plot = p1,
    width = 8,
    height = 6,
    dpi = 600
  )
  
  
### Top 20 miRNAs by number of targets -----------------------------
  top20 <- data_summary %>%
    slice_head(n = 20)
  
  
  p2 <- ggplot(top20, aes(x = reorder(miRNA, n_targets), y = n_targets)) +
    geom_col() +
    coord_flip() +
    labs(
      x = "tae-miRNA",
      y = "Number of predicted target genes",
      title = paste0(
        "Top 20 tae-miRNAs by predicted target number (MFE ≤ -",
        mfe,
        " kcal/mol)"
      )
    ) +
    theme_classic()
  
  
  ggsave(
    filename = paste0("top20_target_miRNAs_mfe-", mfe, ".png"),
    plot = p2,
    width = 8,
    height = 6,
    dpi = 600
  )
  
  

### Cumulative contribution of top miRNAs -----------------------------
  
  data_summary <- data_summary %>%
    mutate(cumulative_percent = cumsum(percent_targets))
  
  
  print(head(
    data_summary[, c("miRNA", "n_targets", "percent_targets", "cumulative_percent")],
    10
  ))
  
}



