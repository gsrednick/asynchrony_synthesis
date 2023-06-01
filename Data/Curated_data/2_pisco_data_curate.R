# Try this with CI database that I have generated

# packages
library(tidyverse)


# survey data accessible at: https://search.dataone.org/view/doi%3A10.6085%2FAA%2FPISCO_kelpforest.1.6 
# upload full pisco dataset: this is a wide-format aggregation of "PISCO_kelpforest_fish.1.3.csv", "PISCO_kelpforest_swath.1.2.csv", and "PISCO_kelpforest_upc.1.2.csv"
pisco_data_full<-read.csv("./Data/PISCO/PISCO_kelpforest_data_aggregated.csv")


# additional CSVs for trophic designation
spp_table<-read.csv("./Data/PISCO/PISCO_kelpforest_taxon_table.1.2_revised_V2.csv")
spp_table_long<-read.csv("./Data/PISCO/PISCO_kelpforest_taxon_table.1.2.csv")
trophic_networks<-read.csv("./Data/PISCO/trophic_network.csv")
network_trophic_position<-read.csv("./Data/PISCO/network_trophic_position.csv")

# filter only for basal prey and consumers of them
basal_networks<-network_trophic_position %>% filter(resource_PA =="1") # only herbivore networks
algae_spp<-spp_table %>% filter(sample_subtype == "ALGAE")
updated_networks<-trophic_networks %>% filter(resource_ID %in% algae_spp$classcode)



# filter sites with greater than 9 years of data
ts_Pisco_test<-pisco_data_full %>% 
  group_by(site) %>%
  summarize(year_count = n()) %>%
  filter(year_count > 9)


pisco_full_long<-pisco_data_full %>% 
  select(-transect) %>%
  filter(site %in% ts_Pisco_test$site) %>% # filter sites with less than 9 years of data
  pivot_longer(-c(site,year),names_to = "species", values_to = "cover") # make long


# spp list
full_spp_list<-data.frame(species = unique(pisco_full_long$species))

cons_spp<-spp_table %>% 
  filter(classcode %in% updated_networks$consumer_ID) %>% 
  select(classcode, sample_subtype) %>% 
  mutate(mode = paste(sample_subtype, "herbivore", sep = "_"),
         mode = tolower(mode),
         species = classcode) %>%
  select(species,mode)

full_spp_list_anno<-merge(full_spp_list,cons_spp,all= T)


spp_table_long$species<-spp_table_long$classcode

spp_table_long_red<-spp_table_long[-c(1,14,16:38)]
dim(spp_table_long)
spp_table_long_nodups<-spp_table_long_red[!duplicated(spp_table_long_red), ]
dim(spp_table_long_nodups)

# have to annotate this dataset
full_spp_list_anno_df<-merge(spp_table_long_nodups,full_spp_list_anno)

write.csv(full_spp_list_anno_df,"./Data/PISCO/spp_for_annotation.csv",row.names = F)
# now bring back in 

full_spp_list_done<-read.csv("./Data/PISCO/spp_annotated.csv")

full_spp_list_red<-full_spp_list_done %>% filter(keep == "yes")

full_mode_list<-full_spp_list_red %>% select(species,mode)
full_mode_list<-full_mode_list[!duplicated(full_mode_list), ]


pisco_full_long_reduced<-pisco_full_long %>% filter(species %in% full_spp_list_red$species)
dim(pisco_full_long_reduced)
dim(full_mode_list)

# something isnt right here....
pisco_data_ready<-merge(pisco_full_long_reduced,full_mode_list)
dim(pisco_data_ready)

pisco_data_ready$system <- "temperate"
pisco_data_ready$dataset <- "PISCO"


write.csv(pisco_data_ready,"./Data/Curated_data/2_pisco_data.csv",row.names = F)


length(unique(pisco_data_ready$site))
unique(pisco_data_ready$year)

pisco_data_ready %>%                    
  group_by(site) %>%          
  summarise(Unique_Elements = n_distinct(year)) %>%
  summarise(mean_year = mean(Unique_Elements),
            sd_year = sd(Unique_Elements))


pisco_data_ready %>%                    
  group_by(mode) %>%          
  summarise(Unique_Elements = n_distinct(species))

# END #