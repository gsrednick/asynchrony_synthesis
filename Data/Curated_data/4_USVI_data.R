# Code in support of "Understanding diversity-synchrony-stability relationships in multitrophic communities"
# in Nature Ecology & Evolution 2024
# Griffin Srednick and Stephen Swearer

# ====== Part B - Synthesis of Long-term Marine Datasets - Dataset 4, Tropical - US Virgin Islands Territorial Coral Reef Monitoring Program (TCRMP) ======

# Try doing this with IMOS dataset 

library(tidyverse)
library(codyn)
library(rfishbase)


# Notes
# Takes a while to run due to the high volume of sites and species to filter

# USVI data load
# fish: https://docs.google.com/spreadsheets/d/1-FvhfMAFs5hNEqmYm0EJd2P23BBYeV7_/edit#gid=1781574307
USVI_fish<-read.csv("./Data/USVI/TCRMP_Master_Fish_Census_Aug2022.xlsx - FishData.csv") # skip metadata rows
USVI_fish_meta<-read.csv("./Data/USVI/TCRMP_Master_Fish_Census_Aug2022.xlsx - FishMetadata.csv") # skip metadata rows

# benthic: https://docs.google.com/spreadsheets/d/1If0CsxbG469kv6U_GyZ0t9ewAtlm03zx/edit#gid=178619352
USVI_benthic<-read.csv("./Data/USVI/TCRMP_Master_Benthic_Cover_Feb2022.xlsx - BenthicData.csv") # skip metadata rows
USVI_benthic_meta<-read.csv("./Data/USVI/TCRMP_Master_Benthic_Cover_Feb2022.xlsx - BenthicCodes.csv") # skip metadata rows


# Datetime is good



# curate fish and get mean abundance across transects per site
names(USVI_fish)

USVI_herb<-USVI_fish %>% filter(TrophicGroup == "herb",
                                Metric == "Abundance")

# fix inconsistent naming
USVI_herb_updated<-USVI_herb %>% mutate(ScientificName = recode(ScientificName, 
                                                                "scarus vetula" = "Scarus vetula",
                                                                "sparisoma aurofrenatum" = "Sparisoma aurofrenatum",
                                                                "sparisoma viride" = "Sparisoma viride"))
USVI_herb_reduced<-USVI_herb_updated %>% 
  group_by(Location,Year,ScientificName) %>% 
  summarise_if(is.numeric,mean)



length(unique(USVI_herb_reduced$Year)) # 19 year dataset


USVI_herb_ready<-USVI_herb_reduced %>% 
  ungroup() %>%
  select(Location,Year,ScientificName,Total) %>%
  mutate(mode = "fish_herbivore",
         dataset = "TCRMP") 

names(USVI_herb_ready) <- c("site","year","species","cover","mode", "dataset")


# curate benthic cover and get mean abundance across transects per site
names(USVI_benthic_meta)
USVI_algal_spp<-USVI_benthic_meta %>% filter(Category == "Macroalgae")


USVI_algal_data<-USVI_benthic[names(USVI_benthic) %in% USVI_algal_spp$Code]

USVI_algal_data<-cbind(USVI_benthic[c(1:9)],USVI_algal_data)


USVI_algae_reduced<-USVI_algal_data %>% 
  group_by(Location,SampleYear) %>% 
  summarise_if(is.numeric,mean, na.rm = T)



# replace NaNs -- this is for taxa that werent sampled...

USVI_algae_reduced<-USVI_algae_reduced %>%
  mutate_all(~replace(., is.nan(.), NA))

 
# now make algae long 
USVI_algae_ready<-USVI_algae_reduced %>% 
  select(-c(SampleMonth,Transect,NoPts)) %>% 
  pivot_longer(-c(Location,SampleYear),
               names_to = "species",
               values_to = "cover") %>%
  mutate(mode = "algae",
         dataset = "TCRMP")



names(USVI_algae_ready) <- c("site","year","species","cover","mode", "dataset")


USVI_data_almost_ready<-rbind(USVI_algae_ready,USVI_herb_ready)

names(USVI_herb_ready)
# filter for sites only with more than 9 years
USVI_ts<-USVI_data_almost_ready %>% 
  group_by(site,year) %>%
  summarize_all(mean) %>%
  summarize(year_count = n()) %>%
  filter(year_count > 9)

dim(USVI_ts)

USVI_data_ready<-USVI_data_almost_ready %>% 
  filter(site %in% USVI_ts$site)
  
USVI_data_ready$system = "tropical"

USVI_data_ready %>% 
  ggplot() +
  geom_point(aes(x = year, y = site))


# data ready 
write.csv(USVI_data_ready,"./Data/Curated_data/4_USVI_data.csv", row.names = F)

length(unique(USVI_data_ready$site))
length(unique(USVI_data_ready$year))

USVI_data_ready %>%                    
  group_by(site) %>%          
  summarise(Unique_Elements = n_distinct(year)) %>%
  summarise(mean_year = mean(Unique_Elements),
            sd_year = sd(Unique_Elements))

USVI_data_ready %>%                    
  group_by(mode) %>%          
  summarise(Unique_Elements = n_distinct(species))

# END #
