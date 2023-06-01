# Try doing this with MCR dataset --- or even better....do it with a bunch of LTER datasets

library(tidyverse)
library(codyn)
library(rfishbase)


# MCR data load
MCR_fish<-read.csv("./Data/MCR/MCR_LTER_Annual_Fish_Survey_20220119.csv") # http://mcrlter.msi.ucsb.edu/cgi-bin/showDataset.cgi?docid=knb-lter-mcr.6
MCR_algae<-read.csv("./Data/MCR/MCR_LTER_Annual_Survey_Benthic_Cover_20220311.csv") # http://mcrlter.msi.ucsb.edu/cgi-bin/showDataset.cgi?docid=knb-lter-mcr.8
MCR_invert<-read.csv("./Data/MCR/MCR_LTER_Annual_Survey_Herbiv_Invert_20220315.csv") # http://mcrlter.msi.ucsb.edu/cgi-bin/showDataset.cgi?docid=knb-lter-mcr.7
MCR_invert_mode<-read.csv("./Data/MCR/ivert_mode.csv")

invert_mode_ready<-MCR_invert_mode %>% filter(keep == "Yes")

unique(MCR_fish$Fine_Trophic)
MCR_invert_list<-unique(MCR_invert$Taxonomy)


# Curate fish
MCR_herbivores<-MCR_fish %>% 
  #filter(Fine_Trophic == "Herbivore/Detritivore") %>% 
  filter(Coarse_Trophic == "Primary Consumer",
         !Taxonomy == "No data") %>% 
  group_by(Year,Site,Habitat,Taxonomy) %>%
  dplyr::summarise(sum = sum(Count,na.rm=T)) %>%
  unite(site_hab, c("Habitat", "Site"))

MCR_herbivores$Year<-as.numeric(MCR_herbivores$Year)
unique(MCR_fish$Coarse_Trophic)

# Curate algae
not_of_interest<-c("No data","Ascidian","Sand","Coral","Soft Coral","Shell Debris","Bare Space","Coral Rubble","Sponge","")

MCR_algae_upd<-MCR_algae %>% 
  mutate(Habitat = recode(Habitat,"Fringing" = "FR",
                                      "Outer 10" = "FO",
                                      "Backreef" = "BA")) %>%
  mutate(Site = recode(Site,"LTER 1" = "1",
                       "LTER 2" = "2",
                       "LTER 3" = "3",
                       "LTER 4" = "4",
                       "LTER 5" = "5",
                       "LTER 6" = "6"))

MCR_algae_red<-MCR_algae_upd %>% 
  filter(!Taxonomy_Substrate_Functional_Group %in% not_of_interest,
         !str_detect(Habitat, "Outer 17")) %>% 
  group_by(Year,Site,Habitat,Taxonomy_Substrate_Functional_Group) %>%
  dplyr::summarise(mean_cov = mean(Percent_Cover,na.rm=T)) %>% 
  unite(site_hab, c("Habitat", "Site"))


MCR_algae_red$Year<-as.numeric(MCR_algae_red$Year)

unique(MCR_algae_red$Taxonomy_Substrate_Functional_Group)


# Inverts
MCR_invert_upd<-MCR_invert %>% 
  mutate(Habitat = recode(Habitat,"Fringing" = "FR",
                          "Outer 10" = "FO",
                          "Backreef" = "BA")) %>%
  mutate(Site = recode(Site,"LTER 1" = "1",
                       "LTER 2" = "2",
                       "LTER 3" = "3",
                       "LTER 4" = "4",
                       "LTER 5" = "5",
                       "LTER 6" = "6"))

MCR_invert_red<-MCR_invert_upd %>% 
  filter(Taxonomy %in% invert_mode_ready$Taxonomy,
         !str_detect(Habitat, "Outer 17")) %>% 
  group_by(Year,Site,Habitat,Taxonomy) %>%
  dplyr::summarise(count = sum(Count,na.rm=T)) %>% 
  unite(site_hab, c("Habitat", "Site"))

MCR_invert_red$Year<-as.numeric(MCR_invert_red$Year)







# make wide

MCR_herbivores_wide<-MCR_herbivores %>% 
  pivot_wider(names_from = Taxonomy,values_from = sum, values_fill = 0)

MCR_herbivores_filled<-MCR_herbivores %>% 
  pivot_wider(names_from = Taxonomy,values_from = sum, values_fill = 0) %>%
  pivot_longer(!c(Year,site_hab),names_to = "Taxonomy", values_to = "sum")

#MCR_herbivores_wide[rowSums(MCR_herbivores_wide[-c(1:2)] < 1), ]
assess<-MCR_herbivores_wide[which(rowSums(MCR_herbivores_wide[-c(1:2)]) < 2), ]





# can generate interaction network from Fishbase; want to know what each herbivore eats...

# function that filters fishes for only herbivores
herbivore_list<-data.frame(Taxonomy = unique(MCR_herbivores$Taxonomy))
herbivore_list[!grep("unidentified", herbivore_list$Taxonomy), ]
herbivore_list_cleaned<-herbivore_list %>% filter(!str_detect(Taxonomy, "unidentified"))


# Check which of these are ACTUALLY herbivores...
herb_spp<-nrow(herbivore_list_cleaned)
real_herbivores <- data.frame()

for (i in 1:herb_spp){
  
  diet_df<-fooditems(herbivore_list_cleaned[i,])
  
  plants_present <- any(str_detect(diet_df$FoodI, "plants"))
  
  print(plants_present)
  
  real_herbivores <- rbind(real_herbivores, data.frame(plants_present = plants_present,
                                                       species = herbivore_list_cleaned[i,]))
}

True_herbivores<-real_herbivores %>% filter(plants_present == "TRUE")

upd_herbivores<-MCR_herbivores_filled %>% filter(Taxonomy %in% True_herbivores$species)




# Rename dataframe columns for merge
MCR_herbfish_ready<-upd_herbivores
MCR_alg_ready<-MCR_algae_red
MCR_invert_ready<-MCR_invert_red

names(MCR_invert_ready) <- c("year","site","species","cover")
MCR_invert_ready$mode <- "invert_herbivore"
MCR_invert_ready$dataset <- "MCR"


names(MCR_herbfish_ready) <- c("year","site","species","cover")
MCR_herbfish_ready$mode <- "fish_herbivore"
MCR_herbfish_ready$dataset <- "MCR"

names(MCR_alg_ready) <- c("year","site","species","cover")
MCR_alg_ready$mode <- "algae"
MCR_alg_ready$dataset <- "MCR"

MCR_data_ready_list<-list(MCR_herbfish_ready,MCR_alg_ready,MCR_invert_ready)
MCR_data_ready <- do.call(rbind, MCR_data_ready_list)

MCR_data_ready$system <- "tropical"

write.csv(MCR_data_ready,"./Data/Curated_data/1_MCR_data.csv", row.names = F)


length(unique(MCR_data_ready$site))

MCR_data_ready %>%                    
  group_by(site) %>%          
  summarise(Unique_Elements = n_distinct(year)) %>%
  summarise(mean_year = mean(Unique_Elements),
            sd_year = sd(Unique_Elements))

MCR_data_ready %>%                    
  group_by(mode) %>%          
  summarise(Unique_Elements = n_distinct(species))

# END # 