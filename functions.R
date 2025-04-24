
func <- read_csv2("./data/functional data/functional data/16S_functions.csv")

names(func)[1] <- "ASV"
p <- 
pivot_longer(func, !ASV, 
             names_to = "f", 
             values_to = "val") %>% 
filter(val == 1) %>% 
select(ASV,f) %>%   
mutate(ASV = as.factor(ASV))

functions <- p %>%
  group_by(ASV) %>%
  summarise(Functions = paste(f, collapse = ", ")) %>% 
  mutate(Functions = replace_na(Functions, "unknown"))
