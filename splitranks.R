##################################################
## Project: KD
## Script purpose: get separate tables with 
## relative abundances for each taxonomic rank
## Date: 2025-04-23
## Author: E. Zeghal
##################################################


######libraries#####
library(tidyverse)
library(Hmisc)

# import of final file created in script 1_data_import.R ------------------
t3 <- read_csv("./data/tidy_data.csv")
t3 <- t3 %>% 
  dplyr::left_join(functions, by = "ASV")%>% 
  mutate(Functions = replace_na(Functions, "unknown"))

# Kingdom -----------------------------------------------------------------
Kingdom <- t3  %>%  
  select(c("Kingdom", "Kingdom_rel_abund_Sample", names(map)))%>% 
  distinct() 
# write_csv(Kingdom, "../Kingdom_PE462.csv")

# Phylum ------------------------------------------------------------------
Phyla <- t3 %>%
  select(c("Phylum", "Phylum_rel_abund_Sample" , names(map)))%>% 
  distinct() 
# write_csv(Phyla, "../phyla_PE462.csv")


# Class -------------------------------------------------------------------
Class <-t3 %>% 
  select(c(
         "Class", "Class_rel_abund_Sample", names(map)) )%>% 
  distinct()  
# write_csv(Class, "../class_PE462.csv")


# Order -------------------------------------------------------------------
Order <- t3 %>% 
  select(c(
           "Order", "Order_rel_abund_Sample", names(map)) )%>% 
  distinct() 
# write_csv(Order, "../order_PE462.csv")



# Family ------------------------------------------------------------------
Family <- t3 %>% 
  select(c("Family", "Family_rel_abund_Sample", names(map)) )%>% 
  distinct()
# write_csv(Family, "../family_PE462.csv")

# Genus -------------------------------------------------------------------
Genus <- t3%>% 
  select(c("Genus", "Genus_rel_abund_Sample", names(map)) )%>%
  distinct()
# write_csv(Genus, "../genus_PE462.csv")


# Species -----------------------------------------------------------------
Species <- t3 %>% 
  select(c("Species", "Species_rel_abund_Sample", names(map)))%>%
  distinct()
# write_csv(Species, "../species_PE462.csv")