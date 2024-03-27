# Code in support of "Understanding diversity-synchrony-stability relationships in multitrophic communities"
# in Nature Ecology & Evolution 2024
# Griffin Srednick and Stephen Swearer

# ====== Part B - Synthesis of Long-term Marine Datasets - Dataset 5, Tropical - AIMS' Long-Term Great Barrier Reef Monitoring Program  (AIMS GBR) ======



# Data provided courtesy of AIMS GBR LTMP - Data owner and primary contact: Michael Emslie
# Genus fish and benthic resolution (data analyzed here) must be sourced from M. Emslie

load("./Data/GBR/orig.ltmp.RData")
load("./Data/GBR/all.benthic.22.RData")
load("./Data/GBR/com_.desc.RData")

# will have to filter benthic for just "MA_XX" functional groups
ltmp_algae<-all.benthic.22 %>% filter(str_detect(COMP_2021, "MA_")| COMP_2021 %in% c("TA","CA")) # filter for algae only
dim(ltmp_algae)

unique(ltmp_algae$COMP_2021) # algal functional groups  

ltmp_algae_red<-ltmp_algae %>% 
  select(REEF_NAME,REPORT_YEAR,SITE_NO,COMP_2021,COVER) %>%
  group_by(REEF_NAME,REPORT_YEAR,SITE_NO,COMP_2021) %>%
  summarize_all(mean,na.rm=T)

ltmp_algae_nearlyready<-ltmp_algae_red %>% 
  select(REEF_NAME,REPORT_YEAR,SITE_NO,COMP_2021,COVER) %>% 
  mutate(site=paste0( REEF_NAME,'_', SITE_NO))

length(unique(ltmp_algae_nearlyready$site))

LTMP_ts<-ltmp_algae_nearlyready %>% 
  group_by(site,REPORT_YEAR) %>%
  summarize_all(mean) %>%
  summarize(year_count = n()) %>%
  filter(year_count > 9)

length(unique(LTMP_ts$site))

ltmp_algae_nearlyready$mode <- "algae"

LTMP_algae_ready<-ltmp_algae_nearlyready %>% 
  filter(site %in% LTMP_ts$site)

length(unique(LTMP_algae_ready$site))


names(LTMP_algae_ready)<-c("REEF_NAME","REPORT_YEAR","SITE_NO","species","cover","site","mode")

# will have to annotate fishes -- actually, just use herbivore functional group
ltmp_fish<- orig.ltmp %>% filter(FUNC_GRP == "herbivore")
dim(ltmp_fish)

unique(ltmp_fish$P_CODE)

# Check if both methods were used
ltmp_method_check<-ltmp_fish %>% 
  group_by(REEF_NAME,SAMPLE_TYPE,REPORT_YEAR) %>%
  reframe() %>%
  group_by(REEF_NAME,REPORT_YEAR) %>%
  dplyr::summarise(count = n()) # only an issue at HAYMAN ISLAND REEF in the first year



ltmp_fish_red<-ltmp_fish %>% 
  select(REEF_NAME,REPORT_YEAR,SITE_NO,TAXA,ABUNDANCE) %>%
  group_by(REEF_NAME,REPORT_YEAR,SITE_NO,TAXA) %>%
  summarize_all(mean,na.rm=T)

ltmp_fish_nearlyready<-ltmp_fish_red %>% 
  select(REEF_NAME,REPORT_YEAR,SITE_NO,TAXA,ABUNDANCE) %>% 
  mutate(site=paste0(REEF_NAME,'_', SITE_NO))

LTMP_fish_ready<-ltmp_fish_nearlyready %>% 
  filter(site %in% LTMP_ts$site)

LTMP_fish_ready$mode <- "fish_herbivore"

names(LTMP_fish_ready)<-c("REEF_NAME","year","SITE_NO","species","cover","site","mode")
names(LTMP_algae_ready)<-c("REEF_NAME","year","SITE_NO","species","cover","site","mode")
#names(LTMP_algae_ready)<-c("site","year","species","cover","mode")
# no invert data



# Bring all together 
LTMP_data_almostready<-rbind(LTMP_fish_ready,LTMP_algae_ready)

LTMP_data_almostready$system = "tropical"
LTMP_data_almostready$dataset = "GBR-LTMP"



# Filter out first year at HAYMAN ISLAND REEF --> only 1 of 2 methods used
LTMP_data_almostready_fishmethod_adj<-LTMP_data_almostready %>% filter(!(str_detect(site,"HAYMAN ISLAND REEF") & year == "1995"))
dim(LTMP_data_almostready_test)
dim(LTMP_data_almostready)

LTMP_data_almostready_test
length(unique(LTMP_data_almostready$site))

LTMP_data_ready<-LTMP_data_almostready_fishmethod_adj %>% ungroup() %>% select(-c(REEF_NAME,SITE_NO))



# data ready 
write.csv(LTMP_data_ready,"./Data/Curated_data/5_GBR_data.csv", row.names = F)

length(unique(LTMP_data_ready$site))
length(unique(LTMP_data_ready$year))

LTMP_data_ready %>%                    
  group_by(site) %>%          
  summarise(Unique_Elements = n_distinct(year)) %>%
  summarise(mean_year = mean(Unique_Elements),
            sd_year = sd(Unique_Elements))

LTMP_data_ready %>%                    
  group_by(mode) %>%          
  summarise(Unique_Elements = n_distinct(species))

LTMP_data_ready %>%
  ggplot(aes(x = year, y = site, color = mode)) +
           geom_point() +
  facet_grid(~mode)



LTMP_species<-data.frame(species = unique(LTMP_data_ready$species))

min(LTMP_data_ready$year)
max(LTMP_data_ready$year)

# END #