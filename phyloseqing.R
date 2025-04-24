##################################################
## Project: PE462
## Script purpose: Importing raw data & cleaning
## Date: 04-01-2023
## Author: E. Zeghal
##################################################


# load libraries ----------------------------------------------------------

library("magrittr")
library("tidyverse")
library("phyloseq")
library("mia")
library("microbiome")
library("Hmisc")

# raw data import ---------------------------------------------------------

tax <- read.csv2("./data/functional data/functional data/16S_tax.csv",
                  row.names = 1, na.strings = c("NA", " ", "")) %>% 
  mutate( across(.cols = everything(), ~str_replace_all(.,c(".__"="")))) %>% 
  mutate_if(is.character,str_squish) %>% 
  mutate(across(.cols = everything(), ~ifelse(as.character(.)!="", ., "NA"))) %>% 
  as.matrix() %>% 
  tax_table()

otu <- as.matrix(read.csv2("./data/functional data/functional data/16S_asv.csv", row.names = 1)) %>% 
  otu_table(taxa_are_rows = TRUE)

map <- sample_data(read.csv("./data/functional data/functional data/mapfile.csv", 
                            row.names = 1, 
                            fileEncoding= "windows-1252",
                            na.strings = c("NA", "", " ")))
physeq_object <-  merge_phyloseq(otu, tax, map)            

# data cleaning 1 ----------------------------------------------------------
summarize_phyloseq(physeq_object)

# Delete poorly assigned ASVs
get_taxa_unique(physeq_object, "Kingdom") 
physeq_object <- subset_taxa(physeq_object, !is.na(Kingdom) & Kingdom %in% c("Bacteria", "Archaea", "Eukaryota"))


# Unify unassigned --------------------------------------------------------

only_unassigned <- function(x) {
  if_else(x %in% c("uncultured","Uncultured","metagenome","Metagenome","unknown","Unknown","NA"),"unassigned", x)
}

taxo <- as.data.frame(physeq_object@tax_table)%>% 
  replace(is.na(.),"unassigned") %>% 
  mutate( across(.cols = everything(), ~only_unassigned(.)))

taxo <- tax_table(as.matrix(taxo))

physeq_object <- merge_phyloseq(physeq_object@otu_table, taxo, map)

# Delete chloroplast and Mitochodria
any((get_taxa_unique(physeq_object, "Order") == "Chloroplast"))
any((get_taxa_unique(physeq_object, "Family") == "Mitochondria"))
physeq_object <- physeq_object %>% 
  subset_taxa(!Order %in% c("Chloroplast")) %>%
  subset_taxa(!Family %in% c("Mitochondria"))

any((get_taxa_unique(physeq_object, "Phylum") == "unassigned"))


# Alpha div ---------------------------------------------------------------

alpha_tab <-microbiome::alpha(physeq_object, index = "all")
alpha_tab <- alpha_tab %>% cbind(data.frame(physeq_object@sam_data))
# write.csv(alpha_tab, file = "./data/functional data/alpha_div_indexes_microbiome_package.csv")
# #phyloseq package
# a_div <- phyloseq::estimate_richness(physeq_object)
# write.csv(a_div, file = "../alpha_div_phyloseq.csv")


# Melt data all -----------------------------------------------------------

source("./data/tidy_psmelt.R")
tidy_physeq_asv <- tidy_psmelt(physeq_object)

tidy_physeq_asv$Species <- if_else(
  (!tidy_physeq_asv$Genus=="unassigned" & tidy_physeq_asv$Species=="unassigned"),
  str_c(tidy_physeq_asv$Genus," sp."), 
  str_c(tidy_physeq_asv$Species))

# Merge replicates and calculate relative abundances ---------------------

t3 <- tidy_physeq_asv  %>% group_by(Sample_ID) %>% mutate(Sample_rel_abund = Abundance / sum(Abundance)) %>% #relative abundance of each otu per sample
  ungroup() %>%

  #Kingdom_section
  group_by(Sample_ID, Kingdom) %>% 
  mutate(Kingdom_rel_abund_Sample = sum(Sample_rel_abund)) %>%  #Kingdom relative abundance per sample 
  ungroup() %>% 
  
  #Phylum_section
  group_by(Sample_ID, Phylum) %>% 
  mutate(Phylum_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  
  #Class_section
  group_by(Sample_ID, Class) %>% 
  mutate(Class_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  
  #Order_section
  group_by(Sample_ID, Order) %>% 
  mutate(Order_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  
  #Family_section
  group_by(Sample_ID, Family) %>% 
  mutate(Family_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  
  #Genus_section
  group_by(Sample_ID, Genus) %>% 
  mutate(Genus_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  
  #Species_section
  group_by(Sample_ID, Species) %>% 
  mutate(Species_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup()  

# write.csv(t3, "./data/tidy_data.csv")
