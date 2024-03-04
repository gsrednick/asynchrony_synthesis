# Code in support of "Understanding diversity-synchrony-stability relationships in multitrophic communities"
# in Nature Ecology & Evolution 2024
# Griffin Srednick and Stephen Swearer

# ====== Part B - Synthesis of Long-term Marine Datasets - Dataset 3, Temperate - Reef Life Survey (RLS), Australian Temperate Reef Collaboration (ATRC)  ======

library(tidyverse)
library(codyn)
library(rfishbase)



# Notes
# Takes a while to run due to the high volume of sites and species to filter

# IMOS data load
# benthic: https://catalogue-imos.aodn.org.au/geonetwork/srv/eng/catalog.search#/metadata/ec424e4f-0f55-41a5-a3f2-726bc4541947
IMOS_benthic<-read.csv("./Data/IMOS/IMOS_-_National_Reef_Monitoring_Network_Sub-Facility_-_Benthic_cover_data_(in_situ_surveys).csv",
                       skip = 71) # skip metadata rows
# fish: https://catalogue-imos.aodn.org.au/geonetwork/srv/eng/catalog.search#/metadata/b273fafa-03d6-4fc2-9acf-39d8c06581e5
IMOS_fish<-read.csv("./Data/IMOS/IMOS_-_National_Reef_Monitoring_Network_Sub-Facility_-_Global_reef_fish_abundance_and_biomass.csv",
                    skip = 71) # skip metadata rows

# invert: https://portal.aodn.org.au/search?uuid=48cf3cb9-caa9-4633-9baa-8bba3c4d904a
IMOS_invert<-read.csv("./Data/IMOS/IMOS_-_National_Reef_Monitoring_Network_Sub-Facility_-_Global_mobile_macroinvertebrate_abundance.csv",
                    skip = 71) # skip metadata rows

# Site list - these are the sites of interest
sites_for_filter<-read.csv("./Data/IMOS/IMOS_sites.csv")

# Invert data
IMOS_invert_mode<-read.csv("./Data/IMOS/IMOS_ivert_mode.csv")
IMOS_invert_herbs <- IMOS_invert_mode %>% filter(keep == "Yes")

unique(IMOS_invert$species_name)
invert_filtering<-IMOS_invert %>% filter(site_code %in% sites_for_filter,
                                         !phylum %in% c("Chordata", "Cnidaria","Platyhelminthes","Echiura"))


unique(invert_filtering$phylum)
ivert_herbs_of_interest<-data.frame(species_name = unique(invert_filtering$species_name))

ivert_herbs_of_interest_V2<-merge(ivert_herbs_of_interest,IMOS_invert_herbs, all = T)
write.csv(ivert_herbs_of_interest_V2,"./Data/IMOS/ivert_herbs_of_interest.csv", row.names = F)

# filter from manual annotation
ivert_filter_df<-read.csv("./Data/IMOS/ivert_herbs_of_interest_import.csv") %>% 
  filter(keep == "Yes")


IMOS_invert_filtered<-IMOS_invert %>% filter(species_name %in% ivert_filter_df$species_name)

IMOS_invert_filtered$mode = "invert_herbivore"

# Bring in a variable for system type
system_list<-read.csv("./Data/IMOS/system_list.csv")

IMOS_fish<-merge(IMOS_fish,system_list)
IMOS_benthic<-merge(IMOS_benthic,system_list)



# datetime conversion
IMOS_benthic$survey_date <- as.Date(IMOS_benthic$survey_date, format = "%Y-%m-%d")
IMOS_fish$survey_date <- as.Date(IMOS_fish$survey_date, format = "%Y-%m-%d")
IMOS_invert_filtered$survey_date <- as.Date(IMOS_invert_filtered$survey_date, format = "%Y-%m-%d")

IMOS_benthic_upd<-IMOS_benthic
IMOS_fish_upd<-IMOS_fish
IMOS_invert_upd<-IMOS_invert_filtered

extractdate <- function(date) {
  day <- format(date, format="%d")
  month <- format(date, format="%m")
  year <- format(date, format="%Y")
  
  cbind(day, month, year)
}

IMOS_benthic_upd<-cbind(IMOS_benthic_upd, extractdate(IMOS_benthic_upd$survey_date))
IMOS_fish_upd<-cbind(IMOS_fish_upd, extractdate(IMOS_fish_upd$survey_date))
IMOS_invert_upd<-cbind(IMOS_invert_upd, extractdate(IMOS_invert_upd$survey_date))



# Summarize to species x site x depth x survey date -- multiple survey dates across years -- then take mean across depths and then annual mean

# fish
IMOS_fish_reduced<-IMOS_fish_upd %>% 
  select(site_code,site_name,survey_date,reporting_name,species_name,depth,total,biomass) %>%
  group_by(site_code,site_name,survey_date,reporting_name,depth,species_name) %>% 
  summarize_all(sum)

length(unique(IMOS_fish_reduced$site_code)) # bigger

# algae
IMOS_benthic_reduced<-IMOS_benthic_upd %>% 
  filter(str_detect(phylum, "phyta") | phylum == "Algae") %>% # filter for algae only
  select(site_code,site_name,survey_date,reporting_name,species_name,depth,total) %>%
  group_by(site_code,site_name,survey_date,reporting_name,depth,species_name) %>% 
  summarize_all(sum) %>% 
  mutate(total_cover = total/10)

length(unique(IMOS_benthic_reduced$site_code)) # smaller; match to me
unique(IMOS_benthic_upd$phylum)

# invert
IMOS_invert_reduced<-IMOS_invert_upd %>%
  #filter(species_name %in% IMOS_invert_herbs$species_name) %>% # filter for algae only
  select(site_code,site_name,survey_date,species_name,depth,mode,total) %>%
  group_by(site_code,site_name,survey_date,depth,species_name,mode) %>% 
  summarize_all(sum) 



# match sites
IMOS_fish_reduced_matched<- subset(IMOS_fish_reduced, site_code %in% IMOS_benthic_reduced$site_code)
IMOS_invert_reduced_matched<- subset(IMOS_invert_reduced, site_code %in% IMOS_benthic_reduced$site_code)


IMOS_fish_reduced_matched$mode = "fish_herbivore"

length(unique(IMOS_fish_reduced_matched$site_code)) # 611 -- updated 606
length(unique(IMOS_invert_reduced_matched$site_code)) # 600 


# filter for fish that are herbivores
fish_list<-data.frame(species = unique(IMOS_fish_reduced_matched$species_name))

write.csv(fish_list,"./Data/IMOS/IMOS_fish_list.csv",row.names = F)

lit_herbivores<-read.csv("./Data/IMOS/IMOS_herb_list.csv")

herb_test_list<-fish_list %>% filter(species %in% lit_herbivores$species) # having issues finding enough herbivores...

fish_spp<-nrow(fish_list)
real_herbivores <- data.frame()

for (i in 1:fish_spp){
  
  diet_df<-fooditems(fish_list[i,])
  #diet_df<-fooditems("Acanthaluteres vittiger")
  
  plants_present <- any(str_detect(diet_df$FoodI, "plants"))
  
  #print(plants_present)
  
  real_herbivores <- rbind(real_herbivores, data.frame(plants_present = plants_present,
                                                       species = fish_list[i,]))
}

True_herbivores<-real_herbivores %>% filter(plants_present == "TRUE")


# manually annotate
herb_manual<-real_herbivores %>% filter(is.na(plants_present)) 
herb_manual <- herb_manual[!duplicated(herb_manual), ]
write.csv(herb_manual,"./Data/IMOS/IMOS_manual.csv",row.names = F)

# bring back in 
herb_annotated<-read.csv("./Data/IMOS/IMOS_fish_annotated.csv")
herb_annotated_filtered<-herb_annotated %>% filter(mode == "herbivore")

#True_herbivores %>% filter(species %in% lit_herbivores$species)
True_herbivores$mode <- "herbivore"
True_herbivores_V2<-rbind(True_herbivores,herb_annotated_filtered)

# save this for quicker processing later
write.csv(True_herbivores,"./Data/IMOS/IMOS_herbivores.csv", row.names = F)

IMOS_herbivores<-IMOS_fish_reduced_matched %>% filter(species_name %in% True_herbivores_V2$species)
dim(IMOS_herbivores)
dim(IMOS_benthic_reduced)

# get means across depths
IMOS_herbivores_site<-IMOS_herbivores %>% 
  ungroup() %>%
  select(site_code,survey_date,species_name,total) %>%
  group_by(site_code,survey_date,species_name) %>%
  summarize_all(mean)

IMOS_algae_site<-IMOS_benthic_reduced %>% 
  ungroup() %>%
  select(site_code,survey_date,species_name,total_cover) %>%
  group_by(site_code,survey_date,species_name) %>%
  summarize_all(mean)

IMOS_invert_site<-IMOS_invert_reduced %>% 
  ungroup() %>%
  select(site_code,survey_date,species_name,total) %>%
  group_by(site_code,survey_date,species_name) %>%
  summarize_all(mean)

# deal with date_time
IMOS_herbivores_site<-cbind(IMOS_herbivores_site, extractdate(IMOS_herbivores_site$survey_date))
IMOS_algae_site<-cbind(IMOS_algae_site, extractdate(IMOS_algae_site$survey_date))
IMOS_invert_site<-cbind(IMOS_invert_site, extractdate(IMOS_invert_site$survey_date))




IMOS_herbivores_site_merge<-IMOS_herbivores_site %>% 
  select(site_code,survey_date,year,month,day,species_name,total) %>%
  mutate(mode = "fish_herbivore")

IMOS_invert_site_merge<-IMOS_invert_site %>% 
  select(site_code,survey_date,year,month,day,species_name,total) %>%
  mutate(mode = "invert_herbivore")

IMOS_algae_site_merge<-IMOS_algae_site %>% 
  mutate(total = total_cover) %>%
  select(site_code,survey_date,year,month,day,species_name,total) %>%
  mutate(mode = "algae")


IMOS_near_complete_list<-list(IMOS_herbivores_site_merge,IMOS_invert_site_merge,IMOS_algae_site_merge)
IMOS_near_complete <- do.call(rbind, IMOS_near_complete_list)

IMOS_nearer_complete<-IMOS_near_complete %>% 
  group_by(site_code,month,day) %>%
  mutate(time_id = as.numeric(format(survey_date, "%j"))/365,
         year_id = as.numeric(year) + time_id)

ts_test<-IMOS_nearer_complete %>% 
  select(-c(mode)) %>%
  pivot_wider(names_from = species_name,
              values_from = total,
              values_fill = 0) %>%
  group_by(site_code) %>%
  summarize(year_count = n()) %>%
  filter(year_count > 9)


IMOS_complete<-IMOS_nearer_complete %>% 
  mutate(year = year_id,
         species = species_name,
         site = site_code,
         dataset = "IMOS",
         cover = total) %>%
  select(site,year,species,dataset, mode,cover)

IMOS_complete$month <-NULL
IMOS_complete$day <-NULL
IMOS_complete$site_code <-NULL

IMOS_data_ready<-IMOS_complete

IMOS_data_rep_check<-IMOS_data_ready %>% 
  group_by(site,mode) %>% 
  summarize(num_spp = n_distinct(species)) %>% ## look for replicates with only one species
  filter(num_spp > 1)

IMOS_data_nearly<-IMOS_complete %>% 
  filter(site %in% IMOS_data_rep_check$site) # only use replicates with more than one species

spplist_check<-data.frame(species = unique(IMOS_data_nearly$species))
# check sites and years for completeness
dup_test<-data.frame(duplicated(IMOS_data_nearly[-6]))

IMOS_data_ready_wide<-IMOS_data_nearly %>%
  select(-c(mode,dataset)) %>%
  pivot_wider(names_from = species,
              values_from = cover,
              values_fill = 0)

ts_test<-IMOS_data_ready_wide %>% 
  group_by(site) %>%
  summarize(year_count = n()) %>%
  filter(year_count > 9)

# assess time series 
IMOS_data_ready_wide %>% 
  filter(site %in% ts_test$site) %>%
  ggplot() +
  geom_point(aes(x = year, y = site))

IMOS_data_ready<-IMOS_data_nearly %>% filter(site %in% ts_test$site)
IMOS_data_ready$system <-"temperate"

site_check<-IMOS_data_ready %>% filter(site == "MIR-S3", mode == "algae")
# yes, there are a number of sites with very high algal richness (> 100 taxa)
unique(site_check$species)

# final check -- remove all "Unidentified xxx" taxa and remove where no species were present
IMOS_data_ready<-IMOS_data_ready %>% filter(!str_detect(species,"Unidentified"),
                                            !str_detect(species,"No species found"))

site_check<-IMOS_data_ready %>% filter(site == "MIR-S3", mode == "algae")
sites_for_filter<-unique(IMOS_data_ready$site)


# data ready 
write.csv(IMOS_data_ready,"./Data/Curated_data/3_IMOS_data.csv", row.names = F)


length(unique(IMOS_data_ready$site))

min(IMOS_data_ready$year)
max(IMOS_data_ready$year)

IMOS_data_ready %>%                    
  group_by(site) %>%          
  summarise(Unique_Elements = n_distinct(year)) %>%
  summarise(mean_year = mean(Unique_Elements),
            sd_year = sd(Unique_Elements))

IMOS_data_ready %>%                    
  group_by(mode) %>%          
  summarise(Unique_Elements = n_distinct(species))

# END #