# Code in support of "Understanding diversity-synchrony-stability relationships in multitrophic communities"
# in Nature Ecology & Evolution 2024
# Griffin Srednick and Stephen Swearer
# DOI: 10.1038/s41559-024-02419-3


# ====== Part B - Synthesis of Long-term Marine Datasets ======



# General pacakges
library(tidyverse)
library(lavaan)
library(cowplot)
library(tidygraph)
library(patchwork)
library(ggExtra)
library(codyn)
library(car)
library(vegan)

# Linear mixed models
library(lme4)
library(lmerTest)
library(ggeffects)
library(kableExtra)

# Lavaan model performance
library(performance)
library(dynamic)
library(ggraph)
library(igraph)

# Data import and merge ####

MCR_data_ready<-read.csv("./Data/Curated_data/1_MCR_data.csv")
pisco_data_ready<-read.csv("./Data/Curated_data/2_pisco_data.csv")
IMOS_data_ready<-read.csv("./Data/Curated_data/3_IMOS_data.csv")
USVI_data_ready<-read.csv("./Data/Curated_data/4_USVI_data.csv")
GBR_data_ready<-read.csv("./Data/Curated_data/5_GBR_data.csv")


names(MCR_data_ready)
names(pisco_data_ready)
names(IMOS_data_ready)
names(USVI_data_ready)
names(GBR_data_ready)

## Merge Datasets together here ####
complete_data_ready_list<-list(pisco_data_ready,MCR_data_ready,IMOS_data_ready,USVI_data_ready,GBR_data_ready)


complete_data_ready <- do.call(rbind.data.frame, complete_data_ready_list)

# DO NOT DELETE
#complete_data_ready_OLD_DONTDEELTE<-complete_data_ready
#write.csv(complete_data_ready_OLD_DONTDEELTE,"/Users/icarus3/Desktop/complete_data_DONTDEELTE.csv")

# ignore warnings --- this is just a summary of sites
site_table<-complete_data_ready %>% 
  group_by(site,dataset,system) %>% 
  summarize_all(mean) %>% 
  select(site,dataset,system)

#OLD_site_check<-complete_data_ready %>%                    
#  group_by(system,dataset) %>%          
#  summarise(Unique_Elements = n_distinct(site))

pre_analysis_site_check<-complete_data_ready %>%  #
  group_by(system,dataset) %>%          
  summarise(Unique_Elements = n_distinct(site))

print(pre_analysis_site_check)

# Compute Hierarchal measures ####

# Loop function to compute:
# Hierarchical species diversity (richness and H') - averaged over time
# Hierarchical species synchrony (Loreau)
# Hierarchical temporal stability (1/CV)

# Function computes this over time for each site; runs for (1) algae, (2) herbivores, and (3) combinations of them
# Output is then curated to bring all estimates in line for each [site] repliate

CR_sync_stab<-lapply(site_list,function(x){
  tryCatch({  
    
    #site_interest = "19131S_1" # for testing
    
    site_interest<-unique(x$site)
    #print(site_interest)
    print(paste0("******CURRENTLY ANALYZING ", site_interest,"******"))
    
    
    cons_df_pre <- complete_data_ready %>% 
      filter(site == site_interest, # filter for site x
             !mode == "") %>% # remove species with no mode: a few algal categories that arent of interest
      group_by(species) %>%
      filter(any(!is.na(cover))) # remove species with all NAs for given site
    
    spp_mode<-cons_df_pre %>% group_by(species,mode) %>% summarise() # make species/mode list
    
    cons_df_filled<- cons_df_pre %>%
      select(-mode) %>%
      pivot_wider(names_from = species, values_from = cover) %>% # make wide to fill NAs
      select(where(function(x) sum(!is.na(x)) > 1)) %>%  # Keep species with more than one observation
      mutate(across(everything(), ~replace(., is.na(.), 0))) %>% # fill species with some NAs with 0; This is only okay because we have removed taxa that are all zeros from the dataset.
      pivot_longer(cols = -c(site,dataset,year,system),  # make long again for analyses 
                   names_to = "species", 
                   values_to = "cover") 
    
    
    cons_df<-merge(cons_df_filled,spp_mode) # merge analysis df with mode data
    
    spp <- unique(cons_df$species) # species list
    
    
    sync_stab_results <- list() # make list for correlations 
    mode<-unique_combination_strings[[1]] # test
    
    for (mode in unique_combination_strings) { # this takes every mode combination for a given site and runs calculations
      tryCatch({  
        
        print(mode) # print mode so we know its working
        
        modes <- unlist(strsplit(mode, "-")) # unlist the combination string so each mode is liste; allows for filtering in following calculations
        
        spp_df<-cons_df %>%  
          filter(mode %in% modes) %>%
          select(species) %>% 
          reframe(species = unique(species))
        
        
        stability<-cons_df %>%       
          filter(mode %in% modes) %>%
          community_stability(.,
                              time.var = "year",
                              abundance.var = "cover") %>%
          as.data.frame() %>%
          mutate(estimate_mode = paste0(mode,"_stability")) 
        
        colnames(stability) <- c("value", "estimate_mode")
        
        # richness 
        mode_rich_TG <- cons_df %>% 
          filter(mode %in% modes) %>%
          group_by(species) %>% 
          summarize_if(is.numeric,sum, na.rm =T) %>%
          filter(cover > 0) %>% 
          summarize(value = n()) %>%
          mutate(estimate_mode = paste0(mode,"_RCH"))
        
        annual_RCH <- cons_df %>%
          filter(mode %in% modes) %>%
          select(-mode) %>%
          group_by(species, year) %>%
          spread(species, cover, fill = 0) %>%
          select(-c(site,system,dataset)) %>% 
          ungroup() %>%
          #select(-c(site, system, dataset)) %>%
          group_by(year) %>%
          mutate(value = vegan::specpool(across(where(is.numeric)))$Species,
                 estimate_mode = paste0(mode, "_AnRCH")) %>%
          distinct(year, .keep_all = TRUE) %>%
          select(year,value,estimate_mode) %>%
          group_by(estimate_mode) %>%
          summarise(value = mean(value, na.rm = TRUE))
        
        
        

        # diversity
        div_df<-cons_df %>% 
          filter(mode %in% modes) %>%
          select(-mode) %>%
          group_by(species) %>%
          summarise(sum_cover = sum(cover,na.rm =T)) %>% 
          spread(species, sum_cover, fill = 0)
        
        diversity<-data.frame(value = vegan::diversity(div_df,"shannon"),
                              estimate_mode = paste0(mode,"_DIV"))
        
        annual_div<-cons_df %>% 
          filter(mode %in% modes) %>%
          select(-mode) %>%
          group_by(species) %>%
          spread(species, cover, fill = 0) %>%
          select(-c(site,system,dataset)) %>% 
          ungroup() %>%
          group_by(year) %>%
          mutate(value = vegan::diversity(across(where(is.numeric)), "shannon"),
                 estimate_mode = paste0(mode,"_AnDIV")) %>%
          distinct(year, .keep_all = TRUE) %>%
          select(year,value,estimate_mode) %>%
          group_by(estimate_mode) %>%
          summarize(value = mean(value,na.rm =T))
        
    
        # VR and Loreau -- VR just for testing 
        mode_VR <- cons_df %>% 
          filter(mode %in% modes) %>%
          variance_ratio(.,
                         time = "year",
                         species.var = "species",
                         abundance.var = "cover",
                         bootnumber = boots) %>%
          mutate(value = VR,
                 estimate_mode = paste0(mode,"_VR")) %>%
          select(estimate_mode,value)
        
        loreau <- cons_df %>% 
          filter(mode %in% modes) %>%
          synchrony(.,
                    time = "year",
                    species.var = "species",
                    abundance.var = "cover",
                    metric = "Loreau") %>%
          as.data.frame()
        
        loreau$estimate_mode = paste0(mode,"_phi")
        colnames(loreau)<-c("value","estimate_mode")
        
        
        
        product<-rbind(stability,mode_rich_TG,mode_VR,diversity,loreau,annual_div,annual_RCH)
        product$mode = mode
        
        sync_stab_results[[mode]] <- product
      }, error = function(e) {sync_stab_results[[mode]] <- NA # put NA for metrics/sites that cant be calculated due to no species present
      })
    }
    
    # Combine the results into a single data frame
    sync_stab_results_df <- do.call(rbind, sync_stab_results)
    
    wide_res<-sync_stab_results_df %>%
      select(-mode) %>%
      pivot_wider(names_from = estimate_mode,values_from = value)   
    
    # Add replicate attributes
    wide_res$site = unique(x$site)
    wide_res$dataset = unique(x$dataset)
    wide_res$system = unique(x$system)
    
    
    synch_stab_by_TG<-cbind(wide_res,mode_rich) # bind it all together
    
    return(synch_stab_by_TG)
    
  }, error = function(e) return(NULL)) # ignore sites for which this cannot be computed; shouldn't actually be any
  
}
)

# combine these lists into a df
CR_sync_stab_df<-do.call(bind_rows, CR_sync_stab)
head(CR_sync_stab_df) # check for completeness
length(unique(CR_sync_stab_df$site)) # check: 497


# Curate into individual dataframes to be merged later

# Temporal stability
stability_df<-CR_sync_stab_df %>% 
  select(contains("_stability"),dataset,system,site) %>%
  pivot_longer(cols = contains("stability"), 
               names_to = "mode_est", 
               values_to = "stability")

# Species richness
RCH_df<-CR_sync_stab_df %>% 
  select(contains("_RCH"),dataset,system,site) %>%
  pivot_longer(cols = contains("_RCH"), 
               names_to = "mode_est", 
               values_to = "RCH")

# Variance ratio
VR_df<-CR_sync_stab_df %>% 
  select(contains("VR"),dataset,system,site) %>%
  pivot_longer(cols = contains("VR"), 
               names_to = "mode_est", 
               values_to = "VR")

# Diversity
DIV_df<-CR_sync_stab_df %>% 
  select(contains("_DIV"),dataset,system,site) %>%
  pivot_longer(cols = contains("_DIV"), 
               names_to = "mode_est", 
               values_to = "DIV")

# Time averaged diversity
AnDIV_df<-CR_sync_stab_df %>% 
  select(contains("AnDIV"),dataset,system,site) %>%
  pivot_longer(cols = contains("AnDIV"), 
               names_to = "mode_est", 
               values_to = "AnDIV")

# Time averaged richness
AnRCH_df<-CR_sync_stab_df %>% 
  select(contains("AnRCH"),dataset,system,site) %>%
  pivot_longer(cols = contains("AnRCH"), 
               names_to = "mode_est", 
               values_to = "AnRCH")

# Species synchrony
loreau_df<-CR_sync_stab_df %>% 
  select(contains("_phi"),dataset,system,site) %>%
  pivot_longer(cols = contains("_phi"), 
               names_to = "mode_est", 
               values_to = "loreau")

stability_df$mode_est <- gsub("_stability", "", stability_df$mode_est)
RCH_df$mode_est <- gsub("_RCH", "", RCH_df$mode_est)
VR_df$mode_est <- gsub("_VR", "", VR_df$mode_est)
DIV_df$mode_est <- gsub("_DIV", "", DIV_df$mode_est)
AnDIV_df$mode_est <- gsub("_AnDIV", "", AnDIV_df$mode_est)
loreau_df$mode_est <- gsub("_phi", "", loreau_df$mode_est)
AnRCH_df$mode_est <- gsub("_AnRCH", "", AnRCH_df$mode_est)


# Merge these together
metrics_list <- list(stability_df,RCH_df,VR_df,DIV_df,AnDIV_df,loreau_df,AnRCH_df)
metrics_df<-Reduce(function(x, y) merge(x, y, all=TRUE), metrics_list, accumulate=FALSE)

unique(metrics_df$mode_est)

# New Levels
metrics_df <- metrics_df %>%
  mutate(mode_est_code = dplyr::recode(mode_est,  
                                "invert_herbivore-fish_herbivore" = "CC",
                                "algae-invert_herbivore-fish_herbivore" = "CR_B",
                                "algae-fish_herbivore" = "CR_F",
                                "algae-invert_herbivore" = "CR_I",
                                "fish_herbivore" = "CC",
                                "invert_herbivore" = "CC",
                                "algae" = "RR"))


metrics_df <- metrics_df %>%
  mutate(mode_agg = dplyr::recode(mode_est,  
                           "invert_herbivore-fish_herbivore" = "CC",
                           "algae-invert_herbivore-fish_herbivore" = "CR",
                           "algae-fish_herbivore" = "CR",
                           "algae-invert_herbivore" = "CR",
                           "fish_herbivore" = "CC",
                           "invert_herbivore" = "CC",
                           "algae" = "RR"))

metrics_df$mode_agg <- factor(metrics_df$mode_agg,
                              levels = c("RR","CC","CR"))
# Other levels
metrics_df <- metrics_df %>%
  mutate(mode_written = dplyr::recode(mode_agg,  
                               "CC" = "monotrophic herbivores",
                               "RR" = "monotrophic algae",
                               "CR" = "multitrophic"))






# more curation #

multi_full_metrics<-metrics_df %>% filter(mode_est == "algae-invert_herbivore-fish_herbivore")
algal_metrics<-metrics_df %>% filter(mode_est_code == "RR")
consumer_metrics<-metrics_df %>% filter(mode_est == "invert_herbivore-fish_herbivore")


# Rename levels for formatting dataframe
consumer_metrics_merge <- consumer_metrics %>%
  rename(
    consumer_rich = RCH,
    consumer_Anrich = AnRCH,
    consumer_div = DIV,
    consumer_Andiv = AnDIV,
    consumer_sync = loreau) %>%
  select(site,dataset,system,consumer_div,consumer_Andiv,consumer_rich,consumer_Anrich,consumer_sync)

algae_metrics_merge <- algal_metrics %>%
  rename(
    algae_rich = RCH,
    algae_Anrich = AnRCH,
    algae_div = DIV,
    algae_Andiv = AnDIV,
    algae_sync = loreau) %>%
  select(site,dataset,system,algae_div,algae_Andiv,algae_rich,algae_Anrich,algae_sync)


multi_metrics_merge <- multi_full_metrics %>%
  select(site,dataset,system,stability)

# Check
head(multi_full_metrics)


multi_1st_merge<-merge(consumer_metrics_merge,algae_metrics_merge)
multi_SEM_df<-merge(multi_1st_merge,multi_metrics_merge)
names(multi_SEM_df)

# Plot to see if time aggregated diversity and stepwise diversity are similar
plot(multi_SEM_df$consumer_Andiv~multi_SEM_df$consumer_div) # aggregated is slightly skewed greater for herbivores
plot(multi_SEM_df$algae_Andiv~multi_SEM_df$algae_div) # aggregated is slightly greater for algae

# Plot to see if time aggregated diversity and richness are similar
plot(multi_SEM_df$consumer_Andiv~multi_SEM_df$consumer_Anrich) # loose decay fit
plot(multi_SEM_df$algae_Andiv~multi_SEM_df$algae_Anrich) # loose decay fit

## Assumptions
# Normality
qqnorm(multi_SEM_df$consumer_Andiv)
qqnorm(log(multi_SEM_df$consumer_Andiv + 1))

qqnorm(multi_SEM_df$algae_Andiv)
qqnorm(log(multi_SEM_df$algae_Andiv + 1))

qqnorm(multi_SEM_df$consumer_sync)
qqnorm(log(multi_SEM_df$consumer_sync+1))

qqnorm(multi_SEM_df$algae_sync)
qqnorm(log(multi_SEM_df$algae_sync + 1))

qqnorm(multi_SEM_df$stability)
qqnorm(log(multi_SEM_df$stability + 1))

qqnorm(multi_SEM_df$algae_Anrich)
qqnorm(log(multi_SEM_df$algae_Anrich + 1))

qqnorm(multi_SEM_df$consumer_Anrich)
qqnorm(log(multi_SEM_df$consumer_Anrich + 1))

# H. of variances
# Violates homogeneity of variances for diversity between systems -- log transformation satisfies
leveneTest(consumer_Andiv ~ dataset, data=multi_SEM_df)
leveneTest(algae_Andiv ~ dataset, data=multi_SEM_df)


# Provide transformed variables for testing -- this could be the fit issue
multi_SEM_df$consumer_Andiv_log<-log(multi_SEM_df$consumer_Andiv+1)
multi_SEM_df$algae_Andiv_log<-log(multi_SEM_df$algae_Andiv+1)

multi_SEM_df$consumer_Anrich_log<-log(multi_SEM_df$consumer_Anrich+1)
multi_SEM_df$algae_Anrich_log<-log(multi_SEM_df$algae_Anrich+1)


multi_SEM_df$consumer_sync_log<-log(multi_SEM_df$consumer_sync+1)
multi_SEM_df$algae_sync_log<-log(multi_SEM_df$algae_sync+1)

multi_SEM_df$stability_log<-log(multi_SEM_df$stability+1)



# Structural equation modelling ####

# Count levels for summary
# Replication
multi_SEM_df %>% 
  na.omit() %>%
  group_by(dataset,site,system) %>% 
  dplyr::summarise() %>%
  group_by(system,dataset) %>%
  dplyr::summarize(count = n())


# CHECK THIS!!!!!
# system    dataset  count
# <chr>     <chr>    <int>
#   1 temperate IMOS       130
# 2 temperate PISCO       93
# 3 tropical  GBR-LTMP   138
# 4 tropical  MCR         18
# 5 tropical  TCRMP       33

# Dataset temporal resolution
multi_SEM_df_nonas<-na.omit(multi_SEM_df) # We renove NAs for formal analysis -- replication reflects this

TS_summary<-complete_data_ready %>% 
  filter(site %in% multi_SEM_df_nonas$site) %>%
  group_by(dataset,site,system,year) %>% 
  dplyr::summarise() %>%
  group_by(dataset,site,system) %>%
  dplyr::summarize(TS_range = max(year)-min(year)) %>%
  group_by(dataset,system) %>%
  dplyr::summarise(mean = mean(TS_range),
                   SD = sd(TS_range)) %>%
  dplyr::mutate_if(is.numeric, format, 2)


year_range<-complete_data_ready %>% 
  filter(site %in% multi_SEM_df_nonas$site) %>%
  group_by(dataset,site,system,year) %>% 
  dplyr::summarise() %>%
  group_by(dataset,system) %>%
  dplyr::summarize(max = max(year),
                   min = min(year))


# Species count per dataset

species_summary<-complete_data_ready %>% 
  filter(site %in% multi_SEM_df_nonas$site) %>%
  group_by(dataset,system,species,mode) %>% 
  dplyr::summarise() %>%
  group_by(dataset,system,mode) %>%
  dplyr::summarize(count = n())



## Model structures ####

# Fully saturated, "just-identified" model
mult_SEM_model<-
  '
consumer_sync_log ~ algae_Andiv + consumer_Andiv
algae_sync_log ~ algae_Andiv + consumer_Andiv
algae_sync_log ~~ consumer_sync_log
stability_log ~ consumer_sync_log + algae_sync_log + algae_Andiv + consumer_Andiv

'


# Seperate for tropical
trop_multi_SEM_df<-multi_SEM_df_nonas %>% filter(system == "tropical")

# Scale diversity within system
trop_multi_SEM_df$consumer_Andiv_scale<-scale(trop_multi_SEM_df$consumer_Andiv)
trop_multi_SEM_df$algae_Andiv_scale<-scale(trop_multi_SEM_df$algae_Andiv)

trop_multi_SEM_df$consumer_Andiv_scale<-trop_multi_SEM_df$consumer_Andiv_scale[,1]
trop_multi_SEM_df$algae_Andiv_scale<-trop_multi_SEM_df$algae_Andiv_scale[,1]

#leveneTest(consumer_Andiv_scale ~ dataset, data=trop_multi_SEM_df) # good
#leveneTest(algae_Andiv_scale ~ dataset, data=trop_multi_SEM_df) # good



# Main model
multi_model_trop<-sem(mult_SEM_model,
                      data = trop_multi_SEM_df, # changed to use a dataframe with no NAs, either df will produce same result
                      estimator = "ML")  

summary(multi_model_trop)

multi_model_trop_fit<-data.frame(model = lavInspect(multi_model_trop, "fit")[c("chisq","df", "pvalue", "rmsea","cfi","tli")])
multi_model_trop_fit
multi_model_trop_df<-as.data.frame(multi_model_trop_fit)
modindices(multi_model_trop, standardized = FALSE, minimum.value = 5)


equivTest(multi_model_trop) # is model a good fit?
varTable(multi_model_trop)
#modificationindices(multi_model_trop, minimum.value = 2) 
trop_model_perf<-as.data.frame(model_performance(multi_model_trop)) 
plot_matrix(residuals(multi_model_trop, type="cor")$cov)




multi_model_trop_df_plot<-multi_model_trop_df %>% mutate(variables = row.names(.)) %>% filter(variables %in% c("chisq","df","pvalue","rmsea"))

param_summary_all_trop<-standardizedSolution(multi_model_trop, type="std.all") # standardized
SEM_trop_plotting<-param_summary_all_trop %>% filter(!rhs == lhs,
                                                     !lhs == rhs)
SEM_trop_plotting<-na.omit(SEM_trop_plotting) # remove covariances that arent of interest for plotting
SEM_trop_plotting$from = SEM_trop_plotting$rhs
SEM_trop_plotting$to = SEM_trop_plotting$lhs

SEM_trop_rsq<-parameterEstimates(multi_model_trop,rsquare = TRUE) %>% filter(op == "r2")

# Rename variables
multi_model_trop_df_plot$model <- sprintf("%.2f", multi_model_trop_df_plot$model)
combined_annotation_trop <- paste(multi_model_trop_df_plot$variables, "=", multi_model_trop_df_plot$model, collapse = "; \n")


trop_covariance<-SEM_trop_plotting %>% 
  dplyr::filter(from == "algae_sync_log" & to == "consumer_sync_log") %>% 
  mutate(to = fct_recode(to, "algae_sync_log" = "consumer_sync_log"),
         from = fct_recode(from, "consumer_sync_log" = "algae_sync_log"))


SEM_trop_plotting <- SEM_trop_plotting %>% 
  mutate(from = fct_recode(from,
                           "herbivore synchrony" = "consumer_sync_log",
                           "algal synchrony" = "algae_sync_log",
                           "algal diversity" = "algae_Andiv",
                           "herbivore diversity" = "consumer_Andiv"),
         to = fct_recode(to,
                         "herbivore synchrony" = "consumer_sync_log",
                         "algal synchrony" = "algae_sync_log",
                         "multitrophic stability" = "stability_log"))




SEM_trop_plotting<-SEM_trop_plotting %>% mutate(direction = ifelse(est.std < 0, "-", "+"),
                                                sig = ifelse(pvalue < 0.05, "yes","no"))

SEM_trop_table<-SEM_trop_plotting[c("to","from","est.std","se","pvalue")]
write.csv(SEM_trop_table,"./Tables/SEM_trop.csv", row.names = F)

hs_graph_trop <- as_tbl_graph(SEM_trop_plotting, directed = T)

coord_positions <- data.frame(name = c("algal diversity","herbivore diversity",
                                       "herbivore synchrony","algal synchrony",
                                       "multitrophic stability"),
                              y = c(15,15,
                                    10,10,
                                    5),
                              x = c(2,3,
                                    3,2,
                                    2.5))

trop_coord_ordered<-data.frame(name = row.names(data.frame(hs_graph_trop[4])))
trop_coord<-merge(trop_coord_ordered,coord_positions, sort = F)

sem_path_tropical<-
  ggraph(hs_graph_trop,trop_coord) + 
  theme_bw() +
  removeGrid() +
  geom_edge_link(aes(edge_color = ifelse(!is.na(sig) & sig == "yes", est.std, NA),
                     alpha = ifelse(!is.na(sig) & sig == "yes", 1,0.5),
                     label = ifelse(!is.na(sig) & sig == "yes", round(est.std,2), NA)),
                 arrow = arrow(length = unit(3, 'mm'),angle = 20,
                               ends = "last",
                               type = "open"),
                 start_cap = circle(12, 'mm'),
                 end_cap = square(25, 'mm'),
                 width = 3,
                 angle_calc = 'along',
                 label_colour = "black",
                 label_size = 2.6) +
  labs(edge_color = "Coef. estimate", linetype = "none", fill = "variable", title = "Tropical") +
  geom_node_point(aes(fill = name), fill = "black",size = 30, pch =21) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(.6,0.07),
        legend.direction = "horizontal",
        legend.key.width=unit(0.8, 'cm'),
        legend.background = element_blank()) +
  guides(alpha = "none", fill = "none") +
  scale_edge_linetype(guide = "none") +
  scale_edge_alpha(guide = 'none') +
  scale_edge_colour_gradientn(colors = c(low = "blue",mid = "white",high = "red"), 
                              values = scales::rescale(c(-0.5, -0.15,0, 0.15,0.5)),
                              limits = c(-0.5, 0.5), oob = scales::squish) +
  geom_node_text(aes(label = stringr::str_wrap(name, 8)), color = "white") +
  lims(x = c(1.5,3.5),
       y = c(2,17))


# temperate
temper_multi_SEM_df<-multi_SEM_df_nonas %>% filter(system == "temperate")

multi_model_temper<-sem(mult_SEM_model,
                        data = temper_multi_SEM_df)


summary(multi_model_temper)
multi_model_temper_fit<-data.frame(model = lavInspect(multi_model_temper, "fit")[c("chisq","df", "pvalue", "rmsea","cfi","tli")])
multi_model_temper_df<-as.data.frame(multi_model_temper_fit)
multi_model_temper_df
multi_model_temper_df_plot<-multi_model_temper_df %>% mutate(variables = row.names(.)) %>% filter(variables %in% c("chisq","df","pvalue","rmsea"))


multi_model_temper_df_plot$model <- sprintf("%.2f", multi_model_temper_df_plot$model)
combined_annotation_temper <- paste(multi_model_temper_df_plot$variables, "=", multi_model_temper_df_plot$model, collapse = "; \n")


param_summary_all_temper<-standardizedSolution(multi_model_temper, type="std.all") # standardized


SEM_temper_rsq<-parameterEstimates(multi_model_temper,rsquare = TRUE) %>% filter(op == "r2") %>% mutate(est.std = est) %>% select(-est)

SEM_temper_plotting<-param_summary_all_temper %>% filter(!rhs == lhs,
                                                         !lhs == rhs)

SEM_temper_plotting<-na.omit(SEM_temper_plotting) # remove covariances that arent of interest for plotting


SEM_temper_plotting$from = SEM_temper_plotting$rhs
SEM_temper_plotting$to = SEM_temper_plotting$lhs



temper_covariance<-SEM_temper_plotting %>% 
  filter(from == "algae_sync_log" & to == "consumer_sync_log") %>% 
  mutate(to = fct_recode(to, "algae_sync_log" = "consumer_sync_log"),
         from = fct_recode(from, "consumer_sync_log" = "algae_sync_log"))

SEM_temper_plotting<-rbind(SEM_temper_plotting,temper_covariance)

SEM_temper_plotting <- SEM_temper_plotting %>% 
  mutate(from = fct_recode(from,
                           "herbivore synchrony" = "consumer_sync_log",
                           "algal synchrony" = "algae_sync_log",
                           "algal diversity" = "algae_Andiv",
                           "herbivore diversity" = "consumer_Andiv"),
  to = fct_recode(to,
                  "herbivore synchrony" = "consumer_sync_log",
                  "algal synchrony" = "algae_sync_log",
                  "multitrophic stability" = "stability_log"))


SEM_temper_plotting<-SEM_temper_plotting %>% mutate(direction = ifelse(est.std < 0, "-", "+"),
                                                    sig = ifelse(pvalue < 0.05, "yes","no"))

# Write table to CSV
SEM_temper_table<-SEM_temper_plotting[c("to","from","est.std","se","pvalue")]
write.csv(SEM_temper_table,"./Tables/SEM_temper.csv", row.names = F)

hs_graph_temper <- as_tbl_graph(SEM_temper_plotting, directed = T)

# extracted but manually written in Illustrator
temper_coord_ordered<-data.frame(name = row.names(data.frame(hs_graph_temper[4])))


temper_coord<-merge(temper_coord_ordered,coord_positions, sort = F)


sem_path_temperate<-
  ggraph(hs_graph_temper,temper_coord) + 
  theme_bw() +
  removeGrid() +
  geom_edge_link(aes(edge_color = ifelse(!is.na(sig) & sig == "yes", est.std, NA), 
                     alpha = ifelse(!is.na(sig) & sig == "yes", 1,0.5),
                     label = ifelse(!is.na(sig) & sig == "yes", round(est.std,2), NA)),
                 arrow = arrow(length = unit(3, 'mm'), angle = 20, ends = "last", type = "open"),
                 start_cap = circle(12, 'mm'),
                 end_cap = square(25, 'mm'),
                 width = 3,
                 angle_calc = 'along',
                 label_colour = "black",
                 label_size = 2.6) +
  labs(edge_color = "Coef. estimate", linetype = "none", fill = "variable", title = "Temperate") +
  geom_node_point(aes(fill = name), fill = "black",size = 30, pch =21) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.6,0.07),
        legend.direction = "horizontal",
        legend.key.width=unit(0.8, 'cm'),
        legend.background = element_blank()) +
  guides(alpha = "none", fill = "none") +
  scale_edge_linetype(guide = "none") +
  scale_edge_alpha(guide = 'none') +
  scale_edge_colour_gradientn(colors = c(low = "blue",mid = "white",high = "red"), 
                              values = scales::rescale(c(-0.5, -0.15,0, 0.15,0.5)),
                              limits = c(-0.5, 0.5), oob = scales::squish) +
  geom_node_text(aes(label = stringr::str_wrap(name, 8)), color = "white") +
  lims(x = c(1.5,3.5),
       y = c(2,17))



## Figure 4 ####
multi_sem_plot<- sem_path_tropical + sem_path_temperate + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("./Figures/SEM_plot_R2.pdf",
       plot = multi_sem_plot,
       width = 10,
       height =6)









# Linear mixed models ####



## Figure 3 - boxplots showing difference in measures between trophic groups and systems  ####

diversity_comp_df<-multi_SEM_df %>% 
  na.omit() %>%
  dplyr::select(site,dataset,system,consumer_Andiv,algae_Andiv) %>%
  pivot_longer(cols = -c(site,dataset,system),  # make long again for analyses 
               names_to = "group", 
               values_to = "diversity") 

richness_comp_df<-multi_SEM_df %>% 
  na.omit() %>%
  dplyr::select(site,dataset,system,consumer_rich,algae_rich) %>%
  pivot_longer(cols = -c(site,dataset,system),  # make long again for analyses 
               names_to = "group", 
               values_to = "richness") 


synchrony_comp_df<-multi_SEM_df %>% 
  na.omit() %>%
  dplyr::select(site,dataset,system,consumer_sync_log,algae_sync_log) %>%
  pivot_longer(cols = -c(site,dataset,system),  # make long again for analyses 
               names_to = "group", 
               values_to = "synchrony") 

stability_comp_df<-multi_SEM_df %>% 
  na.omit() %>%
  dplyr::select(site,dataset,system,stability,stability_log)


diversity_comp_df$group <- gsub("_Andiv", "", diversity_comp_df$group)
synchrony_comp_df$group <- gsub("_sync_log", "", synchrony_comp_df$group)


## Models ####

# diversity across groups and datasets
lmer_div<-lmer(diversity~group * system + (1|dataset), diversity_comp_df)

summary(lmer_div)
anova(lmer_div)
ranova(lmer_div)
performance(lmer_div)


div_predict<-ggpredict(lmer_div, terms = c("group","system"))

# synchrony across groups and datasets
lmer_sync<-lmer(synchrony~group * system + (1|dataset), synchrony_comp_df)
summary(lmer_sync)
anova(lmer_sync)
ranova(lmer_sync)
performance(lmer_sync)


sync_predict<-ggpredict(lmer_sync, terms = c("group","system"))

# Stability across groups and datasets
lmer_stab<-lmer(log(stability+1)~system + (1|dataset), stability_comp_df)
summary(lmer_stab)
anova(lmer_stab)
ranova(lmer_stab)
performance(lmer_stab)

stab_predict<-ggpredict(lmer_stab, terms = c("system"))



richness_plot<-ggplot(richness_comp_df,aes(x = group, y= log(richness+1), fill = system)) +
  geom_boxplot(alpha = 0.4) +
  #geom_point() +
  scale_fill_manual(values = c("purple","yellow3")) +
  scale_color_manual(values = c("purple","yellow3")) +  
  theme_bw() +
  labs(x = "Trophic group", y = "richness") +
  scale_x_discrete(labels = c("algae","herbivores"))+
  removeGrid()


diversity_plot<-ggplot(diversity_comp_df,aes(x = group, y= diversity, fill = system)) +
  geom_boxplot(alpha = 0.4) +
  #geom_point() +
  scale_fill_manual(values = c("purple","yellow3")) +
  scale_color_manual(values = c("purple","yellow3")) +  
  theme_bw() +
  labs(x = "Trophic group", y = "H'") +
  scale_x_discrete(labels = c("algae","herbivores"))+
  removeGrid()




synchrony_plot<-ggplot(synchrony_comp_df,aes(x = group, y= synchrony, fill = system)) +
  geom_boxplot(alpha = 0.4) +
  #geom_point() +
  scale_fill_manual(values = c("purple","yellow3")) +
  scale_color_manual(values = c("purple","yellow3")) +  
  theme_bw() +
  labs(x = "Trophic group", y = "Synchrony (log transformed)") +
  scale_x_discrete(labels = c("algae","herbivores"))+
  removeGrid()

stability_plot<-ggplot(stability_comp_df,aes(x = system, y= log(stability+1), fill = system)) +
  geom_boxplot(alpha = 0.4) +
  #geom_point() +
  scale_fill_manual(values = c("purple","yellow3")) +
  scale_color_manual(values = c("purple","yellow3")) +  
  theme_bw() +
  labs(x = "System", y = "1/CV (log transformed)") +
  removeGrid()


comparison_plot<-diversity_plot + synchrony_plot + stability_plot + plot_layout(guides = "collect") & theme(legend.position = "top")

ggsave("./Figures/comparison_plot_V2.pdf",
       plot = comparison_plot,
       width = 9,
       height =3.8)


## Linear Mixed Effects Tables ####
library(htmltools)
library(kableExtra)

round_p_value <- function(p_value) {
  if (p_value < 0.001) {
    return("; p < 0.001")
  } else {
    return(paste0("; p = ", formatC(p_value, format = "f", digits = 3)))
  }
}

makeStars <- function(x){
  stars <- c("****", "***", "**", "*", "ns")
  vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}

# (1) Diversity across groups and systems
lmer_div_table<-data.frame(anova(lmer_div))
lmer_div_table$variable <-rownames(lmer_div_table)
rownames(lmer_div_table) <-c("trophic group (TG)","system","system:TG")
lmer_div_table$variable <-c("trophic group (TG)","system","system:TG")
#lmer_div_table_export<-lmer_div_table

# Add a column with subscript text using the variable
DenDF_div_est <- sapply(round(lmer_div_table$DenDF,2), function(x) HTML(paste0("<i>F</i><sub>", x, "</sub>")))
lmer_div_table$rounded_p_value <- sapply(lmer_div_table$Pr..F., round_p_value)
lmer_div_table$p_star <- makeStars(lmer_div_table$Pr..F.)
#stab_rich_table$label <- paste0("</p>", DenDF_values, "</p>", "  =  ",round(stab_rich_table$F.value,2),stab_rich_table$rounded_p_value)
lmer_div_table$label <- paste0("</p>", DenDF_div_est, "</p>", "  =  ",round(lmer_div_table$F.value,2),"<sup>", lmer_div_table$p_star, "</sup>")
lmer_div_table$Diversity <-lmer_div_table$label

# (2) Synchrony across groups and systems
lmer_sync_table<-data.frame(anova(lmer_sync))
lmer_sync_table$variable <-rownames(lmer_sync_table)
rownames(lmer_sync_table) <-c("trophic group (TG)","system","system:TG")
lmer_sync_table$variable <-c("trophic group (TG)","system","system:TG")
#lmer_div_table_export<-lmer_div_table

# Add a column with subscript text using the variable
DenDF_sync_est <- sapply(round(lmer_sync_table$DenDF,2), function(x) HTML(paste0("<i>F</i><sub>", x, "</sub>")))
lmer_sync_table$rounded_p_value <- sapply(lmer_sync_table$Pr..F., round_p_value)
lmer_sync_table$p_star <- makeStars(lmer_sync_table$Pr..F.)
#stab_rich_table$label <- paste0("</p>", DenDF_values, "</p>", "  =  ",round(stab_rich_table$F.value,2),stab_rich_table$rounded_p_value)
lmer_sync_table$label <- paste0("</p>", DenDF_sync_est, "</p>", "  =  ",round(lmer_sync_table$F.value,2),"<sup>", lmer_sync_table$p_star, "</sup>")
lmer_sync_table$Synchrony <-lmer_sync_table$label

# (3) Stability across groups and systems
lmer_stab_table<-data.frame(anova(lmer_stab))
lmer_stab_table$variable <-rownames(lmer_stab_table)
rownames(lmer_stab_table) <-c("system")
lmer_stab_table$variable <-c("system")
#lmer_div_table_export<-lmer_div_table

# Add a column with subscript text using the variable
DenDF_stab_est <- sapply(round(lmer_stab_table$DenDF,2), function(x) HTML(paste0("<i>F</i><sub>", x, "</sub>")))
lmer_stab_table$rounded_p_value <- sapply(lmer_stab_table$Pr..F., round_p_value)
lmer_stab_table$p_star <- makeStars(lmer_stab_table$Pr..F.)
#stab_rich_table$label <- paste0("</p>", DenDF_values, "</p>", "  =  ",round(stab_rich_table$F.value,2),stab_rich_table$rounded_p_value)
lmer_stab_table$label <- paste0("</p>", DenDF_stab_est, "</p>", "  =  ",round(lmer_stab_table$F.value,2),"<sup>", lmer_stab_table$p_star, "</sup>")
lmer_stab_table$Stability <-lmer_stab_table$label



# Combined 
lmer_div_merge<-lmer_div_table %>% dplyr::select(variable,Diversity)
lmer_sync_merge<-lmer_sync_table %>% dplyr::select(variable,Synchrony)
lmer_stab_merge<-lmer_stab_table %>% dplyr::select(variable,Stability)

lmer_div_merge$id  <- 1:nrow(lmer_div_merge)
lmer_sync_merge$id  <- 1:nrow(lmer_sync_merge)

completed_table_list<-list(lmer_div_merge,lmer_sync_merge,lmer_stab_merge)
completed_table<-Reduce(function(x, y) merge(x, y, all=TRUE), completed_table_list, accumulate=FALSE)

# finished table
completed_table<-completed_table[order(completed_table$id), ]
completed_table
row.names(completed_table)<-NULL
completed_table$id<-NULL

completed_table$Stability <- replace(completed_table$Stability, is.na(completed_table$Stability), "-")

completed_table %>% 
  kbl(format = "html", escape = F, col.names = c("Fixed effects", "Diversity","Synchrony","Stability")) %>%
  kable_classic(html_font = "Times") %>%
  footnote(general = "Results from linear mixed effects modeling (LME) assessing variation in Shannon-Wiener diversity, species synchrony, and stability across system and trophic group (TG) as fixed effects. Dataset was treated as a random effect. Subscripts are denominator degrees of freedom from Satterhwaite approximation. Superscripts represent significance based on p-value: ns is p > 0.05; * p < 0.05; ** p < 0.01; *** p < 0.001; **** p < 0.0001",
           escape = T) %>%
  kable_styling() 


# Supplemental -- analyses and figures for marine synthesis ####


## SEM with species richness instead of diversity  ####
mult_SEM_RCH_model<-
  '
consumer_sync_log ~ algae_Anrich_log + consumer_Anrich_log
algae_sync_log ~ algae_Anrich_log + consumer_Anrich_log
algae_sync_log ~~ consumer_sync_log
stability_log ~ consumer_sync_log + algae_sync_log + algae_Anrich_log + consumer_Anrich_log

'

## Tropical - good model ####

multi_model_RCH_trop<-sem(mult_SEM_RCH_model,
                          data = trop_multi_SEM_df, # changed to use a dataframe with no NAs
                          estimator = "ML",
                          se = "robust") 

summary(multi_model_RCH_trop)
multi_model_trop_RCH_fit<-data.frame(model = lavInspect(multi_model_RCH_trop, "fit")[c("chisq","df", "pvalue", "rmsea","cfi","tli")])
multi_model_trop_RCH_fit 
modindices(multi_model_RCH_trop, standardized = FALSE, minimum.value = 5) 
AIC(multi_model_RCH_trop)


SEM_trop_rch_rsq<-parameterEstimates(multi_model_RCH_trop,rsquare = TRUE) %>% filter(op == "r2") %>% mutate(est.std = est) %>% select(-est)

multi_model_trop_df_rch<-as.data.frame(multi_model_trop_RCH_fit)
multi_model_trop_df_rch_plot<-multi_model_trop_df_rch %>% mutate(variables = row.names(.)) %>% filter(variables %in% c("chisq","df","pvalue","rmsea"))

param_summary_all_trop_rch<-standardizedSolution(multi_model_RCH_trop, type="std.all") # standardized
SEM_trop_plotting_rch<-param_summary_all_trop_rch %>% filter(!rhs == lhs,
                                                             !lhs == rhs)
SEM_trop_plotting_rch<-na.omit(SEM_trop_plotting_rch) # remove covariances that arent of interest for plotting
SEM_trop_plotting_rch$from = SEM_trop_plotting_rch$rhs
SEM_trop_plotting_rch$to = SEM_trop_plotting_rch$lhs



# Rename variables
multi_model_trop_df_rch_plot$model <- sprintf("%.2f", multi_model_trop_df_rch_plot$model)
combined_annotation_trop_RCH <- paste(multi_model_trop_df_rch_plot$variables, "=", multi_model_trop_df_rch_plot$model, collapse = "; \n")

trop_covariance_RCH<-SEM_trop_plotting_rch %>% 
  dplyr::filter(from == "algae_sync_log" & to == "consumer_sync_log") %>% 
  mutate(to = fct_recode(to, "algae_sync_log" = "consumer_sync_log"),
         from = fct_recode(from, "consumer_sync_log" = "algae_sync_log"))


SEM_trop_plotting_rch <- SEM_trop_plotting_rch %>% 
  mutate(from = fct_recode(from,
                           "herbivore synchrony" = "consumer_sync_log",
                           "algal synchrony" = "algae_sync_log",
                           "algal richness" = "algae_Anrich_log",
                           "herbivore richness" = "consumer_Anrich_log"),
  to = fct_recode(to,
                  "herbivore synchrony" = "consumer_sync_log",
                  "algal synchrony" = "algae_sync_log",
                  "multitrophic stability" = "stability_log"))




SEM_trop_plotting_rch<-SEM_trop_plotting_rch %>% mutate(direction = ifelse(est.std < 0, "-", "+"),
                                                        sig = ifelse(pvalue < 0.05, "yes","no"))

SEM_trop_table_rch<-SEM_trop_plotting_rch[c("to","from","est.std","se","pvalue")]
write.csv(SEM_trop_table_rch,"./Tables/SEM_trop_richness.csv", row.names = F)

hs_graph_trop_rch <- as_tbl_graph(SEM_trop_plotting_rch, directed = T)

coord_positions_rch <- data.frame(name = c("algal richness","herbivore richness",
                                           "herbivore synchrony","algal synchrony",
                                           "multitrophic stability"),
                                  y = c(15,15,
                                        10,10,
                                        5),
                                  x = c(2,3,
                                        3,2,
                                        2.5))

rch_trop_coord_ordered<-data.frame(name = row.names(data.frame(hs_graph_trop_rch[4])))
rch_trop_coord<-merge(rch_trop_coord_ordered,coord_positions_rch, sort = F)

sem_path_tropical_rch<-
  ggraph(hs_graph_trop_rch,rch_trop_coord) + 
  theme_bw() +
  removeGrid() +
  geom_edge_link(aes(edge_color = ifelse(!is.na(sig) & sig == "yes", est.std, NA), 
                     alpha = ifelse(!is.na(sig) & sig == "yes", 1,0.5),
                     label = ifelse(!is.na(sig) & sig == "yes", round(est.std,2), NA)),
                 arrow = arrow(length = unit(3, 'mm'),angle = 20,
                               ends = "last",
                               type = "open"),
                 start_cap = circle(12, 'mm'),
                 end_cap = square(25, 'mm'),
                 width = 3,
                 angle_calc = 'along',
                 label_colour = "black",
                 label_size = 2.6) +
  labs(edge_color = "Coef. estimate", linetype = "none", fill = "variable", title = "Tropical") +
  geom_node_point(aes(fill = name), fill = "black",size = 30, pch =21) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(.6,0.07),
        legend.direction = "horizontal",
        legend.key.width=unit(0.8, 'cm'),
        legend.background = element_blank()) +
  guides(alpha = "none", fill = "none") +
  scale_edge_linetype(guide = "none") +
  scale_edge_alpha(guide = 'none') +
  scale_edge_colour_gradientn(colors = c(low = "blue",mid = "white",high = "red"), 
                              values = scales::rescale(c(-0.5, -0.15,0, 0.15,0.5)),
                              limits = c(-0.5, 0.5), oob = scales::squish) +
  geom_node_text(aes(label = stringr::str_wrap(name, 8)), color = "white") +
  lims(x = c(1.5,3.5),
       y = c(2,17)) 








## Temperate  ####
multi_model_RCH_temper<-sem(mult_SEM_RCH_model,
                            estimator = "ML",
                            se = "robust",
                            data = temper_multi_SEM_df)

summary(multi_model_RCH_temper)
multi_model_temper_RCH_fit<-data.frame(model = lavInspect(multi_model_RCH_temper, "fit")[c("chisq","df", "pvalue", "rmsea","cfi","tli")])
multi_model_temper_RCH_fit 
multi_model_temper_rich_df<-as.data.frame(multi_model_temper_RCH_fit)
modindices(multi_model_RCH_temper, standardized = FALSE, minimum.value = 5)
AIC(multi_model_RCH_temper)


multi_model_temper_rich_df_plot<-multi_model_temper_rich_df %>% mutate(variables = row.names(.)) %>% filter(variables %in% c("chisq","df","pvalue","rmsea"))


multi_model_temper_rich_df_plot$model <- sprintf("%.2f", multi_model_temper_rich_df_plot$model)
combined_annotation_temper_rch <- paste(multi_model_temper_rich_df_plot$variables, "=", multi_model_temper_rich_df_plot$model, collapse = "; \n")


param_summary_all_temper_rch<-standardizedSolution(multi_model_RCH_temper, type="std.all") # standardized

SEM_temper_rch_rsq<-parameterEstimates(multi_model_RCH_temper,rsquare = TRUE) %>% filter(op == "r2") %>% mutate(est.std = est) %>% select(-est)


SEM_temper_plotting_rch<-param_summary_all_temper_rch %>% filter(!rhs == lhs,
                                                                 !lhs == rhs)
SEM_temper_plotting_rch<-na.omit(SEM_temper_plotting_rch) # remove covariances that arent of interest for plotting
SEM_temper_plotting_rch$from = SEM_temper_plotting_rch$rhs
SEM_temper_plotting_rch$to = SEM_temper_plotting_rch$lhs


temper_covariance_rch<-SEM_temper_plotting_rch %>% 
  filter(from == "algae_sync_log" & to == "consumer_sync_log") %>% 
  mutate(to = fct_recode(to, "algae_sync_log" = "consumer_sync_log"),
         from = fct_recode(from, "consumer_sync_log" = "algae_sync_log"))

SEM_temper_plotting_rch<-rbind(SEM_temper_plotting_rch,temper_covariance_rch)

SEM_temper_plotting_rch <- SEM_temper_plotting_rch %>% 
  mutate(from = fct_recode(from,
                           "herbivore synchrony" = "consumer_sync_log",
                           "algal synchrony" = "algae_sync_log",
                           "algal richness" = "algae_Anrich_log",
                           "herbivore richness" = "consumer_Anrich_log"),
  to = fct_recode(to,
                  "herbivore synchrony" = "consumer_sync_log",
                  "algal synchrony" = "algae_sync_log",
                  "multitrophic stability" = "stability_log"))


SEM_temper_plotting_rch<-SEM_temper_plotting_rch %>% mutate(direction = ifelse(est.std < 0, "-", "+"),
                                                            sig = ifelse(pvalue < 0.05, "yes","no"))

# Write table to CSV
SEM_temper_table_rch<-SEM_temper_plotting_rch[c("to","from","est.std","se","pvalue")]
write.csv(SEM_temper_table_rch,"./Tables/SEM_temper_richness.csv", row.names = F)

hs_graph_temper_rch <- as_tbl_graph(SEM_temper_plotting_rch, directed = T)

SEM_trop_table_rch
rch_temper_coord_ordered<-data.frame(name = row.names(data.frame(hs_graph_temper_rch[4])))
rch_temper_coord<-merge(rch_temper_coord_ordered,coord_positions_rch, sort = F)


sem_path_temperate_rch<-
  ggraph(hs_graph_temper_rch,rch_temper_coord) + 
  theme_bw() +
  removeGrid() +
  geom_edge_link(aes(edge_color = ifelse(!is.na(sig) & sig == "yes", est.std, NA), 
                     alpha = ifelse(!is.na(sig) & sig == "yes", 1,0.5),
                     label = ifelse(!is.na(sig) & sig == "yes", round(est.std,2), NA)),
                 arrow = arrow(length = unit(3, 'mm'), angle = 20, ends = "last", type = "open"),
                 start_cap = circle(12, 'mm'),
                 end_cap = square(25, 'mm'),
                 width = 3,
                 angle_calc = 'along',
                 label_colour = "black",
                 label_size = 2.6) +
  labs(edge_color = "Coef. estimate", linetype = "none", fill = "variable", title = "Temperate") +
  geom_node_point(aes(fill = name), fill = "black",size = 30, pch =21) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.6,0.07),
        legend.direction = "horizontal",
        legend.key.width=unit(0.8, 'cm'),
        legend.background = element_blank()) +
  guides(alpha = "none", fill = "none") +
  scale_edge_linetype(guide = "none") +
  scale_edge_alpha(guide = 'none') +
  scale_edge_colour_gradientn(colors = c(low = "blue",mid = "white",high = "red"), 
                              values = scales::rescale(c(-0.5, -0.15,0, 0.15,0.5)),
                              limits = c(-0.5, 0.5), oob = scales::squish) +
  geom_node_text(aes(label = stringr::str_wrap(name, 8)), color = "white") +
  lims(x = c(1.5,3.5),
       y = c(2,17)) #+
annotate("text", x = 1.5, y= 3, hjust = 0,label = combined_annotation_temper_rch)



multi_sem_plot_rch<- sem_path_tropical_rch + sem_path_temperate_rch + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("./Figures/SEM_plot_richness_R2.pdf",
       plot = multi_sem_plot_rch,
       width = 10,
       height =6)


# SEM - Decomposing the between trophic groups effects ####

algal_metrics<-metrics_df %>% filter(mode_est_code == "RR")
consumer_metrics<-metrics_df %>% filter(mode_est == "invert_herbivore-fish_herbivore")


consumer_metrics_merge_STAB <- consumer_metrics %>%
  rename(
    consumer_rich = RCH,
    consumer_Anrich = AnRCH,
    consumer_div = DIV,
    consumer_Andiv = AnDIV,
    consumer_sync = loreau,
    consumer_stability = stability) %>%
  dplyr::select(site,dataset,system,consumer_div,consumer_Andiv,consumer_rich,consumer_Anrich,consumer_sync,consumer_stability)

algae_metrics_merge_STAB <- algal_metrics %>%
  rename(
    algae_rich = RCH,
    algae_Anrich = AnRCH,
    algae_div = DIV,
    algae_Andiv = AnDIV,
    algae_sync = loreau,
    algae_stability = stability) %>%
  dplyr::select(site,dataset,system,algae_div,algae_Andiv,algae_rich,algae_Anrich,algae_sync,algae_stability)


trophic_SEM_df<-merge(consumer_metrics_merge_STAB,algae_metrics_merge_STAB)

# transformations
trophic_SEM_df$consumer_sync_log<-log(trophic_SEM_df$consumer_sync+1)
trophic_SEM_df$algae_sync_log<-log(trophic_SEM_df$algae_sync+1)

trophic_SEM_df$consumer_stability_log<-log(trophic_SEM_df$consumer_stability+1)
trophic_SEM_df$algae_stability_log<-log(trophic_SEM_df$algae_stability+1)

# Modified model for this -- these are exactly the same models for temperate and tropical 
trophic_SEM_model_trop<-
  '
consumer_sync_log ~ algae_Andiv + consumer_Andiv
algae_sync_log ~ algae_Andiv + consumer_Andiv
algae_sync_log ~~ consumer_sync_log
algae_stability_log ~ consumer_sync_log + algae_sync_log + algae_Andiv + consumer_Andiv
consumer_stability_log ~ consumer_sync_log + algae_sync_log + algae_Andiv + consumer_Andiv

'

trophic_SEM_model_temper<-
  '
consumer_sync_log ~ algae_Andiv + consumer_Andiv
algae_sync_log ~ algae_Andiv + consumer_Andiv
algae_sync_log ~~ consumer_sync_log
algae_stability_log ~ consumer_sync_log + algae_sync_log + consumer_Andiv + algae_Andiv 
consumer_stability_log ~ consumer_sync_log + algae_sync_log + consumer_Andiv + algae_Andiv 

'


# Trop model
trophic_SEM_df_trop<-trophic_SEM_df %>% filter(system == "tropical")

trophic_model_trop<-sem(trophic_SEM_model_trop,
                        data = trophic_SEM_df_trop, # changed to use a dataframe with no NAs
                        estimator = "ML")

summary(trophic_model_trop)

trophic_model_trop_fit<-data.frame(model = lavInspect(trophic_model_trop, "fit")[c("chisq","df", "pvalue", "rmsea","cfi","tli")])
trophic_model_trop_fit 
modindices(trophic_model_trop, standardized = FALSE, minimum.value = 5) 
AIC(trophic_model_trop)

trophic_model_trop_fit_plot<-trophic_model_trop_fit %>% mutate(variables = row.names(.)) %>% filter(variables %in% c("chisq","df","pvalue","rmsea"))


trophic_model_trop_fit_plot$model <- sprintf("%.2f", trophic_model_trop_fit_plot$model)
combined_annotation_trop_trophic <- paste(trophic_model_trop_fit_plot$variables, "=", trophic_model_trop_fit_plot$model, collapse = "; \n")


param_summary_trophic_model_trop<-standardizedSolution(trophic_model_trop, type="std.all") # standardized
trophic_model_trop_rsq<-parameterEstimates(trophic_model_trop,rsquare = TRUE) %>% filter(op == "r2") %>% mutate(est.std = est) %>% select(-est)



SEM_trop_plotting_trophic<-param_summary_trophic_model_trop %>% filter(!rhs == lhs,
                                                                       !lhs == rhs)
SEM_trop_plotting_trophic<-na.omit(SEM_trop_plotting_trophic) # remove covariances that arent of interest for plotting
SEM_trop_plotting_trophic$from = SEM_trop_plotting_trophic$rhs
SEM_trop_plotting_trophic$to = SEM_trop_plotting_trophic$lhs


trop_covariance_trophic<-SEM_trop_plotting_trophic %>% 
  filter(from == "algae_sync_log" & to == "consumer_sync_log") %>% 
  mutate(to = fct_recode(to, "algae_sync_log" = "consumer_sync_log"),
         from = fct_recode(from, "consumer_sync_log" = "algae_sync_log"))

SEM_trop_plotting_trophic<-rbind(SEM_trop_plotting_trophic,trop_covariance_trophic)

SEM_trop_plotting_trophic <- SEM_trop_plotting_trophic %>% 
  mutate(from = fct_recode(from,
                           "herbivore synchrony" = "consumer_sync_log",
                           "algal synchrony" = "algae_sync_log",
                           "algal diversity" = "algae_Andiv",
                           "herbivore diversity" = "consumer_Andiv",
                           "herbivore stability" = "consumer_stability_log"),
  to = fct_recode(to,
                  "herbivore synchrony" = "consumer_sync_log",
                  "algal synchrony" = "algae_sync_log",
                  "algal stability" = "algae_stability_log",
                  "herbivore stability" = "consumer_stability_log"))


SEM_trop_plotting_trophic<-SEM_trop_plotting_trophic %>% mutate(direction = ifelse(est.std < 0, "-", "+"),
                                                                sig = ifelse(pvalue < 0.05, "yes","no"))

# Write table to CSV
SEM_trop_plotting_trophic_table<-SEM_trop_plotting_trophic[c("to","from","est.std","se","pvalue")]
write.csv(SEM_trop_plotting_trophic_table,"./Tables/SEM_trop_trophic.csv", row.names = F)

hs_graph_trop_trophic <- as_tbl_graph(SEM_trop_plotting_trophic, directed = T)

# Coords
coord_positions_troph <- data.frame(name = c("algal diversity","herbivore diversity",
                                             "herbivore synchrony","algal synchrony",
                                             "herbivore stability", "algal stability"),
                                    y = c(15,15,
                                          10,10,
                                          5,5),
                                    x = c(2,3,
                                          3,2,
                                          2.8,2.2))

#SEM_trop_table_rch
troph_trop_coord_ordered<-data.frame(name = row.names(data.frame(hs_graph_trop_trophic[4])))
troph_trop_coord<-merge(troph_trop_coord_ordered,coord_positions_troph, sort = F)


sem_path_trop_troph<-
  ggraph(hs_graph_trop_trophic,troph_trop_coord) +
  theme_bw() +
  removeGrid() +
  geom_edge_link(aes(edge_color = ifelse(!is.na(sig) & sig == "yes", est.std, NA), 
                     alpha = ifelse(!is.na(sig) & sig == "yes", 1,0.5),
                     label = ifelse(!is.na(sig) & sig == "yes", round(est.std,2), NA)),
                 arrow = arrow(length = unit(3, 'mm'), angle = 20, ends = "last", type = "open"),
                 start_cap = circle(12, 'mm'),
                 end_cap = square(25, 'mm'),
                 width = 3,
                 angle_calc = 'along',
                 label_colour = "black",
                 label_size = 2.6) +
  labs(edge_color = "Coef. estimate", linetype = "none", fill = "variable", title = "Tropical") +
  geom_node_point(aes(fill = name), fill = "black",size = 30, pch =21) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.6,0.07),
        legend.direction = "horizontal",
        legend.key.width=unit(0.8, 'cm'),
        legend.background = element_blank()) +
  guides(alpha = "none", fill = "none") +
  scale_edge_linetype(guide = "none") +
  scale_edge_alpha(guide = 'none') +
  scale_edge_colour_gradientn(colors = c(low = "blue",mid = "white",high = "red"), 
                              values = scales::rescale(c(-0.5, -0.15,0, 0.15,0.5)),
                              limits = c(-0.5, 0.5), oob = scales::squish) +
  geom_node_text(aes(label = stringr::str_wrap(name, 8)), color = "white") +
  lims(x = c(1.5,3.5),
       y = c(2,17))


# Temperate model

trophic_SEM_df_temper<-trophic_SEM_df %>% filter(system == "temperate")

trophic_model_temper<-sem(trophic_SEM_model_temper,
                          data = trophic_SEM_df_temper,
                          estimator = "ML")

summary(trophic_model_temper)
trophic_model_temper_fit<-data.frame(model = lavInspect(trophic_model_temper, "fit")[c("chisq","df", "pvalue", "rmsea","cfi","tli")])
trophic_model_temper_fit 
modindices(trophic_model_temper, standardized = FALSE, minimum.value = 5) 
AIC(trophic_model_temper)


trophic_model_temper_fit_plot<-trophic_model_temper_fit %>% mutate(variables = row.names(.)) %>% filter(variables %in% c("chisq","df","pvalue","rmsea"))


trophic_model_temper_fit_plot$model <- sprintf("%.2f", trophic_model_temper_fit_plot$model)
combined_annotation_temper_trophic <- paste(trophic_model_temper_fit_plot$variables, "=", trophic_model_temper_fit_plot$model, collapse = "; \n")


param_summary_trophic_model_temper<-standardizedSolution(trophic_model_temper, type="std.all") # standardized

trophic_model_temper_rsq<-parameterEstimates(trophic_model_temper,rsquare = TRUE) %>% filter(op == "r2") %>% mutate(est.std = est) %>% select(-est)


SEM_temper_plotting_trophic<-param_summary_trophic_model_temper %>% filter(!rhs == lhs,
                                                                           !lhs == rhs)
SEM_temper_plotting_trophic<-na.omit(SEM_temper_plotting_trophic) # remove covariances that arent of interest for plotting
SEM_temper_plotting_trophic$from = SEM_temper_plotting_trophic$rhs
SEM_temper_plotting_trophic$to = SEM_temper_plotting_trophic$lhs


temper_covariance_trophic<-SEM_temper_plotting_trophic %>% 
  filter(from == "algae_sync_log" & to == "consumer_sync_log") %>% 
  mutate(to = fct_recode(to, "algae_sync_log" = "consumer_sync_log"),
         from = fct_recode(from, "consumer_sync_log" = "algae_sync_log"))

SEM_temper_plotting_trophic<-rbind(SEM_temper_plotting_trophic,temper_covariance_trophic)

SEM_temper_plotting_trophic <- SEM_temper_plotting_trophic %>% 
  mutate(from = fct_recode(from,
                           "herbivore synchrony" = "consumer_sync_log",
                           "algal synchrony" = "algae_sync_log",
                           "algal diversity" = "algae_Andiv",
                           "herbivore diversity" = "consumer_Andiv",
                           "herbivore stability" = "consumer_stability_log"),
  to = fct_recode(to,
                  "herbivore synchrony" = "consumer_sync_log",
                  "algal synchrony" = "algae_sync_log",
                  "algal stability" = "algae_stability_log",
                  "herbivore stability" = "consumer_stability_log"))


SEM_temper_plotting_trophic<-SEM_temper_plotting_trophic %>% mutate(direction = ifelse(est.std < 0, "-", "+"),
                                                                    sig = ifelse(pvalue < 0.05, "yes","no"))

# Write table to CSV
SEM_temper_plotting_trophic_table<-SEM_temper_plotting_trophic[c("to","from","est.std","se","pvalue")]
write.csv(SEM_temper_plotting_trophic_table,"./Tables/SEM_temper_trophic.csv", row.names = F)

hs_graph_temper_trophic <- as_tbl_graph(SEM_temper_plotting_trophic, directed = T)

# Coords
coord_positions_troph <- data.frame(name = c("algal diversity","herbivore diversity",
                                             "herbivore synchrony","algal synchrony",
                                             "herbivore stability", "algal stability"),
                                    y = c(15,15,
                                          10,10,
                                          5,5),
                                    x = c(2,3,
                                          3,2,
                                          2.8,2.2))

#SEM_trop_table_rch
troph_temper_coord_ordered<-data.frame(name = row.names(data.frame(hs_graph_temper_trophic[4])))
troph_temper_coord<-merge(troph_temper_coord_ordered,coord_positions_troph, sort = F)


sem_path_temperate_troph<-
  ggraph(hs_graph_temper_trophic,troph_temper_coord) + 
  theme_bw() +
  removeGrid() +
  geom_edge_link(aes(edge_color = ifelse(!is.na(sig) & sig == "yes", est.std, NA),
                     alpha = ifelse(!is.na(sig) & sig == "yes", 1,0.5),
                     label = ifelse(!is.na(sig) & sig == "yes", round(est.std,2), NA)),
                 arrow = arrow(length = unit(3, 'mm'), angle = 20, ends = "last", type = "open"),
                 start_cap = circle(12, 'mm'),
                 end_cap = square(25, 'mm'),
                 width = 3,
                 angle_calc = 'along',
                 label_colour = "black",
                 label_size = 2.6) +
  labs(edge_color = "Coef. estimate", linetype = "none", fill = "variable", title = "Temperate") +
  geom_node_point(aes(fill = name), fill = "black",size = 30, pch =21) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.6,0.07),
        legend.direction = "horizontal",
        legend.key.width=unit(0.8, 'cm'),
        legend.background = element_blank()) +
  guides(alpha = "none", fill = "none") +
  scale_edge_linetype(guide = "none") +
  scale_edge_alpha(guide = 'none') +
  scale_edge_colour_gradientn(colors = c(low = "blue",mid = "white",high = "red"), 
                              values = scales::rescale(c(-0.5, -0.15,0, 0.15,0.5)),
                              limits = c(-0.5, 0.5), oob = scales::squish) +
  geom_node_text(aes(label = stringr::str_wrap(name, 8)), color = "white") +
  lims(x = c(1.5,3.5),
       y = c(2,17))




trophic_sem_plot<- sem_path_trop_troph + sem_path_temperate_troph + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("./Figures/SEM_plot_trophic_R2.pdf",
       plot = trophic_sem_plot,
       width = 10,
       height =6)



# END #