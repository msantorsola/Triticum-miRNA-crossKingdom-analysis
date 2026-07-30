## Plot enriched pathways of interest

## *************************************
#R version 4.3.2

## Settings ----------------------------
set.seed(1)

rm(list=ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(clusterProfiler, org.Hs.eg.db, AnnotationDbi, dplyr, ggplot2, ggbreak, tidyr, ReactomePA,rrvgo, readr)  


## Plot enriched pathways ----------------------------

kegg_summary <- read.delim("data/KEGG_results/summary_KEGG_significant_terms.csv", sep=',')

p_kegg <- ggplot(kegg_summary, aes(x = factor(mfe_threshold), y = significant_KEGG)) +
  geom_bar(stat = "identity", fill = "Gold2") +
  labs(title = "",
       x = "MFE Threshold",
       y = "Number of KEGG Pathways") +
  theme_bw(base_size = 20)+
     theme(
   panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # bordo nero
   axis.title = element_text(size = 22, color = "black"),                   # titoli assi
   axis.text = element_text(size = 15, color = "black")                       # titoli assi     axis.text = element_text(size = 12, color = "black"),                    # numeri assi
    #plot.title = element_text(size = 16, color = "black", hjust = 0.5)      # titolo centrale
   )


ggsave("data/KEGG_results/plot_KEGG_summary.png", plot = p_kegg, width = 7, height = 5, dpi = 300)


## KEGG bubble plots by mfe Thresholds ----------------------------

kegg <- read.delim("data/KEGG_results/KEGG_mfe_-25.csv", sep=',')
kegg <- subset(kegg, kegg$p.adjust < 0.05)
#top20 <- kegg[order(kegg$p.adjust), ][1:20, ]

selected_kegg_mfe25 <- c("Cell adhesion molecules","ECM-receptor interaction","Cytokine-cytokine receptor interaction","Neuroactive ligand-receptor interaction","Viral protein interaction with cytokine and cytokine receptor","Platelet activation","T cell receptor signaling pathway","Th1 and Th2 cell differentiation","Th17 cell differentiation","C-type lectin receptor signaling pathway","Fc epsilon RI signaling pathway","Fc gamma R-mediated phagocytosis","Leukocyte transendothelial migration","Leukocyte transendothelial migration - Homo sapiens (human)","Rheumatoid arthritis","Primary immunodeficiency","Gastric acid secretion","Protein digestion and absorption","Cholesterol metabolism","IL-17 signaling pathway","Hematopoietic cell lineage","Intestinal immune network for IgA production","Inflammatory mediator regulation of TRP channels","Cushing syndrome","AGE-RAGE signaling pathway in diabetic complications","Insulin resistance","Type II diabetes mellitus","Alcoholic liver disease")

kegg_sub <- subset(kegg, kegg$Description %in% selected_kegg_mfe25)
kegg_sub <- kegg_sub[order(kegg_sub$p.adjust), ]

kegg_sub$Description <- factor(kegg_sub$Description, levels = rev(kegg_sub$Description[order(kegg_sub$p.adjust)]))
kegg_sub$GeneRatio <- sapply(strsplit(as.character(kegg_sub$GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))

png("KEGG_results/KEGG_mfe-25_selected_22082025.png", height = 20, width = 35, res=600, units='cm')
    
ggplot(kegg_sub, aes(x = GeneRatio, y = Description), p.adjust) +
  geom_point(aes(size=Count, colour = p.adjust))  +
  scale_size_continuous(range = c(5, 12)) +
    theme_bw() +
     theme(#panel.grid.major.y = element_blank(),
          #panel.grid.minor.x=element_blank(),  #blank grid
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          axis.title.x=element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20)) +
          # scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
          ylab("") +
          scale_color_gradient(low = "darkorange2", high = "gold2", name = "p.adjust")
          ggtitle("Enriched KEGG pathways (mfe: -25)") 
                    #theme(strip.background = element_blank(), strip.text = element_blank()) 
          #theme(legend.position="top")


dev.off() 

## MFE -30

kegg <- read.delim("data/KEGG_results/KEGG_mfe_-30.csv", sep=',')
kegg <- subset(kegg, kegg$p.adjust < 0.05)

selected_kegg_mfe30 <- c("Focal adhesion","Signaling pathways regulating pluripotency of stem cells","Cell adhesion molecules","Cytokine-cytokine receptor interaction","Cushing syndrome","Type II diabetes mellitus","Insulin resistance","AGE-RAGE signaling pathway in diabetic complications","Gastric acid secretion","C-type lectin receptor signaling pathway","Fc gamma R-mediated phagocytosis","T cell receptor signaling pathway","Chemokine signaling pathway","Th17 cell differentiation","Th1 and Th2 cell differentiation","Platelet activation","Fc epsilon RI signaling pathway","B cell receptor signaling pathway","Leukocyte transendothelial migration","Inflammatory mediator regulation of TRP channels")

kegg_sub <- subset(kegg, kegg$Description %in% selected_kegg_mfe30)
kegg_sub <- kegg_sub[order(kegg_sub$p.adjust), ]
kegg_sub$Description <- factor(kegg_sub$Description, levels = rev(kegg_sub$Description[order(kegg_sub$p.adjust)]))

#to plot number instead of the ratio!
kegg_sub$GeneRatio <- sapply(strsplit(as.character(kegg_sub$GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))


kegg_sub$Description_wrapped <- str_wrap(kegg_sub$Description, width = 60)

# factor levels
kegg_sub$Description_wrapped <- factor(
  kegg_sub$Description_wrapped,
  levels = str_wrap(levels(kegg_sub$Description), width = 60)
)

png("data/KEGG_results/KEGG_mfe-30_selected_pub_2.png", height = 20, width = 35, res=600, units='cm')
    
ggplot(kegg_sub, aes(x = GeneRatio, y = Description_wrapped)) +
  geom_point(aes(size=Count, colour = p.adjust)) +
  scale_size_continuous(range = c(5, 12)) +
    theme_bw() +
     theme(#panel.grid.major.y = element_blank(),
          #panel.grid.minor.x=element_blank(),  #blank grid
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.title = element_text(size = 21, color = "black"),
          legend.text = element_text(size = 21, color = "black"),
          axis.title.x=element_text(size = 25, color = "black"),
          axis.text.x = element_text(size = 19, color = "black"),
          axis.text.y = element_text(size = 25, color = "black")) +
          # scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
          ylab("") +
          scale_color_gradient(low = "darkorange2", high = "gold2", name = "p.adjust") +
          ggtitle("") 
          #theme(strip.background = element_blank(), strip.text = element_blank()) 
          #theme(legend.position="top")


dev.off() 

## MFE -35
kegg <- read.delim("KEGG_results/KEGG_mfe_-35.csv", sep=',')
kegg_sub <- subset(kegg, kegg$p.adjust < 0.05)

#order plot by padjusted
kegg_sub$Description <- factor(kegg_sub$Description, levels = rev(kegg_sub$Description[order(kegg_sub$p.adjust)]))
kegg_sub$GeneRatio <- sapply(strsplit(as.character(kegg_sub$GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))

png("KEGG_results/KEGG_mfe-35_selected.png", height = 20, width = 35, res=600, units='cm')
    
ggplot(kegg_sub, aes(x = GeneRatio, y = Description), p.adjust) +
  geom_point(aes(size=Count, colour = p.adjust))  +
  scale_size_continuous(range = c(5, 12)) +
    theme_bw() +
     theme(#panel.grid.major.y = element_blank(),
          #panel.grid.minor.x=element_blank(),  #blank grid
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          axis.title.x=element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20)) +
          # scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
          ylab("") +
          scale_color_gradient(low = "darkorange2", high = "gold2", name = "p.adjust") +
          ggtitle("Enriched KEGG pathways (mfe: -35)") 
          #theme(strip.background = element_blank(), strip.text = element_blank()) 
          #theme(legend.position="top")


dev.off() 



### GO and KEGG summary

go_summary_wide <- go_summary %>%
  pivot_wider(names_from = ontology, values_from = significant_terms,
              names_prefix = "significant_GO_")

# Merge with KEGG
final_summary <- merge(go_summary_wide, kegg_summary, by = "mfe_threshold")

# Save the final summary table
write.csv(final_summary, "tables/summary_GO_KEGG_significant_terms_13052025.csv", row.names = FALSE)

print(final_summary)


##### ReactomePA

enrich_reactome_by_threshold <- function(threshold) {
  genes <- data_annotated %>%
    filter(mfe <= threshold) %>%
    pull(gene_symbol) %>%
    unique()
  
  entrez_ids <- bitr(genes,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
  
  if (nrow(entrez_ids) == 0) return(NULL)
  
  enrichPathway(gene = entrez_ids$ENTREZID,
             organism = 'human',
             pvalueCutoff = 0.05)
}


if (!dir.exists("Reactome_results")) dir.create("Reactome_results")

reactome_summary <- data.frame()

# for each threshold, run and record Reactome enrichment
for (thresh in thresholds) {
  reactome_result <- enrich_reactome_by_threshold(thresh)
  
  if (!is.null(reactome_result) && nrow(reactome_result@result) > 0) {
    write.csv(reactome_result@result, file = paste0("Reactome_results/Reactome_mfe_", thresh, "_13052025.csv"), row.names = FALSE)
    sig_reactome <- sum(reactome_result@result$p.adjust < 0.05)
  } else {
    sig_reactome <- 0
  }
  
  # Add to summary
  reactome_summary <- rbind(reactome_summary, data.frame(
    mfe_threshold = thresh,
    significant_Reactome = sig_reactome
  ))
}

write.csv(reactome_summary, "Reactome_results/summary_Reactome_significant_terms_13052025.csv", row.names = FALSE)

## plot
reactome_summary <- read.delim("Reactome_results/summary_Reactome_significant_terms_13052025.csv", sep=',')

p_reactome <- ggplot(reactome_summary, aes(x = factor(mfe_threshold), y = significant_Reactome)) +
  geom_bar(stat = "identity", fill = "#3498DB") +
  labs(title = "Significant Reactome Pathways by MFE Threshold",
       x = "MFE Threshold",
       y = "Number of Significant Reactome Pathways") +
  theme_bw(base_size = 12)

ggsave("Reactome_results/plot_Reactome_summary.png", plot = p_reactome, width = 7, height = 5, dpi = 300)


### bubble plots

# MFE -25

reactome <- read.delim("Reactome_results/Reactome_mfe_-25.csv", sep=',')
reactome_sub <- subset(reactome, reactome$p.adjust < 0.05)
selected_pathways_mfe25 <- c("RHO GTPase cycle","Signaling by VEGF","VEGFA-VEGFR2 Pathway","Death Receptor Signaling","Platelet activation, signaling and aggregation","Toll Like Receptor 4 (TLR4) Cascade","Toll Like Receptor TLR1:TLR2 Cascade","Toll Like Receptor 2 (TLR2) Cascade","Toll Like Receptor TLR6:TLR2 Cascade","Toll Like Receptor 10 (TLR10) Cascade","Toll Like Receptor 5 (TLR5) Cascade","Cell-Cell communication","Glycosaminoglycan metabolism","Costimulation by the CD28 family","TRAF6 mediated induction of NFkB and MAP kinases upon TLR7/8 or 9 activation","Class I MHC mediated antigen processing & presentation","Heparan sulfate/heparin (HS-GAG) metabolism","Toll Like Receptor 7/8 (TLR7/8) Cascade","MyD88 dependent cascade initiated on endosome","Antigen processing: Ubiquitination & Proteasome degradation","Neutrophil degranulation","Basigin interactions","Toll Like Receptor 3 (TLR3) Cascade")

reactome_sub <- subset(reactome_sub, reactome_sub$Description %in% selected_pathways_mfe25)
reactome_sub <- reactome_sub[order(reactome_sub$p.adjust), ]
reactome_sub$Description <- factor(reactome_sub$Description, levels = rev(reactome_sub$Description[order(reactome_sub$p.adjust)]))
reactome_sub$GeneRatio <- sapply(strsplit(as.character(reactome_sub$GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))

png("Reactome_results/Reactome_mfe-25_selected.png", height = 20, width = 35, res=600, units='cm')
    
ggplot(reactome_sub, aes(x = GeneRatio, y = Description), p.adjust) +
  geom_point(aes(size=Count, colour = p.adjust))  +
  scale_size_continuous(range = c(5, 12)) +
    theme_bw() +
     theme(#panel.grid.major.y = element_blank(),
          #panel.grid.minor.x=element_blank(),  #blank grid
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18),
          axis.title.x=element_text(size = 18),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18)) +
          # scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
          ylab("") +
          scale_color_gradient(low = "#2F9EC2", high = "cadetblue3", name = "p.adjust") +
          ggtitle("Enriched Reactome pathways (mfe: -25)") 
          #theme(strip.background = element_blank(), strip.text = element_blank()) 
          #theme(legend.position="top")


dev.off() 


## MFE 30
reactome <- read.delim("data/Reactome_results/Reactome_mfe_-30_13052025.csv", sep=',')
reactome_sub <- subset(reactome, reactome$p.adjust < 0.05)
selected_pathways_mfe30 <- c("Platelet activation, signaling and aggregation","Death Receptor Signaling","VEGFA-VEGFR2 Pathway","Signaling by VEGF","Costimulation by the CD28 family","CD28 co-stimulation","Signaling by Erythropoietin","GPVI-mediated activation cascade","Platelet homeostasis","VEGFR2 mediated vascular permeability","Glycosaminoglycan metabolism","Interleukin-3, Interleukin-5 and GM-CSF signaling","Signaling by Interleukins")

reactome_sub <- subset(reactome_sub, reactome_sub$Description %in% selected_pathways_mfe30)
reactome_sub <- reactome_sub[order(reactome_sub$p.adjust), ]

reactome_sub$Description <- factor(reactome_sub$Description, levels = rev(reactome_sub$Description[order(reactome_sub$p.adjust)]))
reactome_sub$GeneRatio <- sapply(strsplit(as.character(reactome_sub$GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
reactome_sub$Description_wrapped <- str_wrap(reactome_sub$Description, width = 60)
reactome_sub$Description_wrapped <- factor(
  reactome_sub$Description_wrapped,
  levels = str_wrap(levels(reactome_sub$Description), width = 60)
)

png("data/Reactome_results/Reactome_mfe-30_selected_pub_2.png", height = 15, width = 30, res=600, units='cm')
    
ggplot(reactome_sub, aes(x = GeneRatio, y = Description_wrapped)) +
  geom_point(aes(size=Count, colour = p.adjust)) +
  scale_size_continuous(range = c(5, 12)) +
    theme_bw() +
     theme(#panel.grid.major.y = element_blank(),
          #panel.grid.minor.x=element_blank(),  #blank grid
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.title = element_text(size = 20, color = "black"),
          legend.text = element_text(size = 20, color = "black"),
          axis.title.x=element_text(size = 23, color = "black"),
          axis.text.x = element_text(size = 18, color = "black"),
          axis.text.y = element_text(size = 23, color = "black")) +
          # scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
          ylab("") +
          scale_color_gradient(low = "#2F9EC2", high = "cadetblue3", name = "p.adjust") +
          ggtitle("") 
          #theme(strip.background = element_blank(), strip.text = element_blank()) 
          #theme(legend.position="top")


dev.off() 


