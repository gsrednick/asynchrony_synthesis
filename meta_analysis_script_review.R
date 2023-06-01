# Script supporting Single-trophic level studies bias our understanding of synchrony and stability in complex communities
# Submitted to Ecology Letters - Synthesis

# Srednick and Swearer XXXX

# Script 2B - synthesis analysis and plotting 

library(tidyverse)
library(codyn)
library(ggpubr)
library(ggExtra)
library(patchwork)
#library(ggpmisc)
library(ggforce)
library(ggh4x)
library(ggtext)

# for lmer
library(lme4)
library(lmerTest)
library(ggeffects)
library(emmeans)

# for tables; not used in ms
library(gt)
library(htmltools)
library(kableExtra)


## ============ Data import and merge ===============

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


complete_data_ready <- do.call(rbind, complete_data_ready_list)

# ignore warnings --- this is just a summary of sites
site_table<-complete_data_ready %>% 
  group_by(site,dataset,system) %>% 
  summarize_all(mean) %>% 
  select(site,dataset,system)

complete_data_ready %>%                    
  group_by(system,dataset) %>%          
  summarise(Unique_Elements = n_distinct(site))



## ============ Calculating synchrony and stability ===============

# low to reduce time when calculating VR; change boot if necessary -- it isn't here, we don't need significance testing for VR estimate
boots = 1

# run variance ratio (synchrony) for herbivores at each site
complete_data_ready$year<-as.numeric(complete_data_ready$year)
site_list<-data.frame(site = unique(complete_data_ready$site))
VR_herb <- data.frame()
site_list_list<- split(site_list, site_list$site)



# For fish herbivores
VR_res<-lapply(site_list_list,function(x){
  tryCatch({  
    cons_df <- complete_data_ready %>% 
      filter(site == x$site,
             mode == "fish_herbivore")
    
    spp <- unique(cons_df$species)
    
    spp_rich<-cons_df %>% group_by(species) %>% 
      summarize_if(is.numeric,sum) %>%
      filter(cover > 0) %>% 
      nrow()
    
    cons_VR<-cons_df %>% 
      filter(species %in% spp & cover != 0) %>%
      variance_ratio(.,
                     time = "year",
                     species.var = "species",
                     replicate.var = "site",
                     abundance.var = "cover",
                     bootnumber = boots)
    
    cons_CV<-cons_df %>% 
      filter(species %in% spp & cover != 0) %>%
      community_stability(.,
                          time = "year",
                          replicate.var = "site",
                          abundance.var = "cover")
    
    
    cons_effect_df <- complete_data_ready %>% 
      filter(site == x$site,
             !mode == "fish_herbivore")
    
    spp_eff <- unique(cons_effect_df$species)
    
    cons_effect_VR<-cons_effect_df %>% 
      filter(species %in% spp_eff & cover != 0) %>%
      variance_ratio(.,
                     time = "year",
                     species.var = "species",
                     replicate.var = "site",
                     abundance.var = "cover",
                     bootnumber = boots)
    
    cons_effect_CV<-cons_effect_df %>% 
      filter(species %in% spp_eff & cover != 0) %>%
      community_stability(.,
                          time = "year",
                          replicate.var = "site",
                          abundance.var = "cover")
    
    VR_herb <- rbind(VR_herb, data.frame(fish_herb_sync = cons_VR$VR, 
                                         site = x$site,
                                         fish_herb_stability=cons_CV$stability, 
                                         fish_herb_rich = spp_rich, 
                                         fish_herb_sync_eff = cons_effect_VR$VR, 
                                         fish_herb_stability_eff = cons_effect_CV$stability))
    
    return(VR_herb)
    
  }, error = function(e) return(NULL))
  
})

# rbind list
VR_res_df<-do.call(rbind.data.frame, VR_res)
head(VR_res_df)




# For invert herbivores
VR_inv_herb <- data.frame()
VR_res_invert<-lapply(site_list_list,function(x){
  tryCatch({  
    cons_df <- complete_data_ready %>% 
      filter(site == x$site,
             mode == "invert_herbivore")
    
    spp <- unique(cons_df$species)
    
    spp_rich<-cons_df %>% group_by(species) %>% 
      summarize_if(is.numeric,sum) %>%
      filter(cover > 0) %>% 
      nrow()
    
    cons_VR<-cons_df %>% 
      filter(species %in% spp & cover != 0) %>%
      variance_ratio(.,
                     time = "year",
                     species.var = "species",
                     replicate.var = "site",
                     abundance.var = "cover",
                     bootnumber = boots)
    
    cons_CV<-cons_df %>% 
      filter(species %in% spp & cover != 0) %>%
      community_stability(.,
                          time = "year",
                          replicate.var = "site",
                          abundance.var = "cover")
    
    
    cons_effect_df <- complete_data_ready %>% 
      filter(site == x$site,
             !mode == "invert_herbivore")
    
    spp_eff <- unique(cons_effect_df$species)
    
    cons_effect_VR<-cons_effect_df %>% 
      filter(species %in% spp_eff & cover != 0) %>%
      variance_ratio(.,
                     time = "year",
                     species.var = "species",
                     replicate.var = "site",
                     abundance.var = "cover",
                     bootnumber = boots)
    
    cons_effect_CV<-cons_effect_df %>% 
      filter(species %in% spp_eff & cover != 0) %>%
      community_stability(.,
                          time = "year",
                          replicate.var = "site",
                          abundance.var = "cover")
    
    VR_inv_herb <- rbind(VR_inv_herb, data.frame(invert_herb_sync = cons_VR$VR, 
                                                 site = x$site,
                                                 invert_herb_stability=cons_CV$stability, 
                                                 invert_herb_rich = spp_rich,
                                                 invert_herb_sync_eff = cons_effect_VR$VR, 
                                                 invert_herb_stability_eff = cons_effect_CV$stability))
    
    return(VR_inv_herb)
    
  }, error = function(e) return(NULL))
  
})

# rbind list
VR_res_invert_df<-do.call(rbind.data.frame, VR_res_invert)
head(VR_res_invert_df)


# run variance ratio for algae at each site
VR_alg <- data.frame()

VR_alg_res<-lapply(site_list_list,function(x){
  tryCatch({  # this is here to treat replicates that have all zeros as NULL -- these reps are not included in analysis
    res_df <- complete_data_ready %>% 
      filter(site == x$site,
             mode == "algae")
    
    spp <- unique(res_df$species)
    spp_rich <- length(unique(res_df$species))
    count_comp_presnt = n_distinct(res_df$species > 0)
    
    spp_rich<-res_df %>% group_by(species) %>% 
      summarize_if(is.numeric,sum) %>%
      filter(cover > 0) %>% 
      nrow()
    
    res_VR<-res_df %>% 
      filter(species %in% spp & cover != 0) %>%
      variance_ratio(.,
                     time = "year",
                     species.var = "species",
                     #replicate.var = "site",
                     abundance.var = "cover",
                     bootnumber = boots)
    
    res_CV<-res_df %>% 
      filter(species %in% spp & cover != 0) %>%
      community_stability(.,
                          time = "year",
                          replicate.var = "site",
                          abundance.var = "cover")
    
    res_effect_df <- complete_data_ready %>% 
      filter(site == x$site,
             !mode == "algae")
    
    spp_eff <- unique(res_effect_df$species)
    
    res_effect_VR<-res_effect_df %>% 
      filter(species %in% spp_eff & cover != 0) %>%
      variance_ratio(.,
                     time = "year",
                     species.var = "species",
                     replicate.var = "site",
                     abundance.var = "cover",
                     bootnumber = boots)
    
    res_effect_CV<-res_effect_df %>% 
      filter(species %in% spp_eff & cover != 0) %>%
      community_stability(.,
                          time = "year",
                          replicate.var = "site",
                          abundance.var = "cover")
    #VR_alg <- rbind(VR_alg, data.frame(alg_sync = res_VR$VR, CI_up_alg =res_VR$upperCI, CI_dn_alg = res_VR$lowerCI, site = x$site,dataset = unique(res_df$dataset), alg_stability = res_CV$stability, alg_rich = spp_rich))
    VR_alg <- rbind(VR_alg, data.frame(alg_sync = res_VR$VR, 
                                       site = x$site,
                                       dataset = unique(res_df$dataset), 
                                       alg_stability = res_CV$stability, 
                                       alg_rich = spp_rich,
                                       alg_sync_eff = res_effect_VR$VR, 
                                       alg_stability_eff = res_effect_CV$stability))
    
    return(VR_alg)
  }, error = function(e) return(NULL))
  
})

# rbind list
VR_alg_res_df<-do.call(rbind.data.frame, VR_alg_res)
head(VR_alg_res_df)


# merge the synch
comm_sync_list<-list(VR_res_invert_df,VR_alg_res_df,VR_res_df)
comm_sync<-Reduce(function(x, y) merge(x, y, all=TRUE), comm_sync_list)



### Synchrony/stability contribution approach from Firkowski et al. 2022
VR_all <- data.frame()

VR_all_res<-lapply(site_list_list,function(x){
  tryCatch({  
    res_df <- complete_data_ready %>% 
      filter(site == x$site)
    
    spp <- unique(res_df$species)
    
    all_sync<-res_df %>% 
      filter(species %in% spp & cover != 0) %>%
      variance_ratio(.,
                     time = "year",
                     species.var = "species",
                     #replicate.var = "site",
                     abundance.var = "cover",
                     bootnumber = boots)
    
    spp_rich<-res_df %>% 
      group_by(species) %>% 
      summarize_if(is.numeric,sum) %>%
      filter(cover > 0) %>% 
      nrow()
    
    
    
    all_CV<-res_df %>% 
      filter(species %in% spp & cover != 0) %>%
      community_stability(.,
                          time = "year",
                          replicate.var = "site",
                          abundance.var = "cover")
    
    VR_all <- rbind(VR_all, data.frame(multi_sync = all_sync$VR, CI_up_all =all_sync$upperCI, CI_dn_all = all_sync$lowerCI, site = x$site,multi_stability = all_CV$stability, all_rich = spp_rich))
    return(VR_all)
  }, error = function(e) return(NULL))
  
})

# rbind list
VR_all_res_df<-do.call(rbind.data.frame, VR_all_res)



all_res_sync<-merge(comm_sync,VR_all_res_df)


# Calculate contributions
multi_sync_df<-all_res_sync %>% mutate(alg_effect = multi_sync - alg_sync_eff, # calculates the effect of the lack of algae on multitrophic synchrony
                                       fish_effect = multi_sync - fish_herb_sync_eff, # calculates the effect of the lack of herbivores on multitrophic synchrony
                                       invert_effect = multi_sync - invert_herb_sync_eff, # calculates the effect of the lack of herbivores on multitrophic synchrony
                                       alg_CV_eff = multi_stability - alg_stability_eff, # calculates the effect of the lack of algae on multitrophic stability
                                       fish_CV_eff = multi_stability - fish_herb_stability_eff,
                                       invert_CV_eff = multi_stability - invert_herb_stability_eff) # calculates the effect of the lack of herbivores on multitrophic stability

names(all_res_sync)
names(multi_sync_df)


# merge synchrony/stability with site table
multi_sync_df<-merge(multi_sync_df,site_table)

write.csv(multi_sync_df,"./Data/Curated_data/data_for_analysis.csv", row.names = F)



# ============ Curate estimates ===============
# Need to remove TCRMP and GBR for this remove NA step because they lack invert_herbivores
invertmissing_data<-multi_sync_df %>% filter(dataset %in% c("GBR-LTMP","TCRMP"))
remaining<-multi_sync_df %>% filter(!dataset %in% c("GBR-LTMP","TCRMP"))

full_sites<-na.omit(remaining)
dim(full_sites)

com_sync_stab<-rbind(invertmissing_data,full_sites)


com_sync_stab %>%                    
  group_by(system) %>%          
  summarise(Unique_Elements = n_distinct(site))

com_sync_stab %>%                    
  group_by(system,dataset) %>%          
  summarise(Unique_Elements = n_distinct(site))


# ============ Plotting and analyses ===============


# tagging function
tag_facet_NP<-function (p, open = "(", close = ")", tag_pool = letters, x = -Inf, 
                        y = Inf, hjust = -0.5, vjust = 1.5, fontface = 2, family = "", 
                        ...) 
{
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(tag_pool[lay$PANEL]), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), 
                ..., hjust = hjust, vjust = vjust, fontface = fontface, 
                family = family, inherit.aes = FALSE) + theme(strip.text = element_blank(), 
                                                              strip.background = element_blank())
}

tag_facet_Wtop<-function (p, open = "(", close = ")", tag_pool = letters, x = -Inf, 
                        y = Inf, hjust = -0.5, vjust = 1.5, fontface = 2, family = "", 
                        ...) 
{
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(tag_pool[lay$PANEL]), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), 
                ..., hjust = hjust, vjust = vjust, fontface = fontface, 
                family = family, inherit.aes = FALSE)
}



# make long
com_sync_stab_long<-com_sync_stab %>% 
  pivot_longer(-c(site,dataset,system), names_to = "metric", values_to = "value")

com_sync_stab_long_test<-com_sync_stab_long %>% 
  filter(metric %in% c("multi_sync","fish_effect","alg_effect","invert_effect"))




# make an aggregated figure for this for use in the main text --- above can be in suppliment to show dataset diffs.
com_assess_agg<-com_sync_stab_long %>% 
  mutate(period = if_else(str_detect(metric,"eff"), "Effect","Estimate")) %>%
  filter(!str_detect(metric,paste(c("CI", "distance", "rich"),collapse = '|'))) %>%
  mutate(color = case_when(
    str_detect(metric, "invert") ~ "blue",
    str_detect(metric, "fish") ~ "red",
    str_detect(metric,"alg") ~ "green",
    TRUE~ "black")) 

# change order
com_assess_agg$period = factor(com_assess_agg$period, levels=c('Estimate','Effect')) 



## ============ Figure 3 ===============
## synchrony

sync_df <-com_assess_agg %>%
  filter(!str_detect(metric,paste(c("stability","CV"),collapse = '|')),
         !metric %in% c("alg_sync_eff","fish_herb_sync_eff","invert_herb_sync_eff"))

## stability

stab_df <- com_assess_agg %>%
  filter(str_detect(metric,paste(c("stability","CV"),collapse = '|')),
         !metric %in% c("alg_stability_eff","fish_herb_stability_eff","invert_herb_stability_eff"))

sync_df_plot <- sync_df %>% 
  mutate(metric = recode_factor(metric,
                                alg_sync = "algae",
                                fish_herb_sync = "fish",
                                invert_herb_sync = "invert",
                                multi_sync = "multi",
                                alg_effect = "algae",
                                fish_effect = "fish",
                                invert_effect = "invert"))

metric_labels <- c(
  "alg_sync" = "algae",
  "fish_herb_sync" = "fish",
  "invert_herb_sync" = "invert",
  "multi_sync" = "multi",
  "alg_effect" = "algae",
  "fish_effect" = "fish",
  "invert_effect" = "invert")

sync_together<-sync_df_plot %>% 
  ggplot(aes(x = system, y = value, fill = metric)) +
  geom_hline(yintercept = 0) +
  #geom_sina(size = 2, pch =21, alpha = 0.4) +
  geom_sina(size = 2, pch =21, alpha = 0.4,position = position_dodge(width = .75)) +
  geom_boxplot(alpha = 0.5,outlier.shape = NA, notch = T) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_fill_manual(
    values = c(algae="#009E73",
               fish = "#CC70A9",
               invert = "#56B4E9",
               multi = "darkgrey"
    )) +
  labs(x = "Organizational scale", y = "Variance-ratio") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    axis.text.x = element_blank()) +
  scale_y_continuous(limits = c(NA, NA),breaks = c(-3,-1,0,1,3,5)) +
  facet_wrap(~period, nrow = 1, scales = "free") +
  removeGrid()

sync_together<-tag_facet_Wtop(sync_together, tag_pool = c("A","B"))


stab_df_plot <- stab_df %>% 
  mutate(metric = recode_factor(metric,
                                alg_stability = "algae",
                                fish_herb_stability = "fish",
                                invert_herb_stability = "invert",
                                multi_stability = "multi",
                                alg_CV_eff = "algae",
                                fish_CV_eff = "fish",
                                invert_CV_eff = "invert"))
stab_together<-stab_df_plot %>% 
  ggplot(aes(x = system, y = value, fill = metric)) +
  geom_hline(yintercept = 0) +
  geom_sina(size = 2, pch =21, alpha = 0.4,position = position_dodge(width = .75)) +
  geom_boxplot(alpha = 0.5,outlier.shape = NA, notch = T) +
  scale_fill_manual(
    values = c(algae="#009E73",
               fish = "#CC70A9",
               invert = "#56B4E9",
               multi = "darkgrey")) +
  labs(x = "Organizational scale", y = "Stability") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()) +
  scale_y_continuous(limits = c(NA, NA),breaks = c(-4,-2,0,2,5,10,15,20,25)) +
  facet_wrap(~period, scales = "free") +
  removeGrid()

stab_together<-tag_facet_Wtop(stab_together, tag_pool = c("C","D"))


system_plot<-sync_together + stab_together + plot_layout(nrow = 2,guides = "collect") & theme(legend.position = "bottom")



ggsave("./Figures/Figure_3.pdf",
       plot=system_plot,
       width = 8,
       height = 7)

ggsave("./Figures/Figure_3.TIFF",
       plot=system_plot,
       width = 8,
       height = 7)

length(unique(stab_df_plot$site)) # 266 sites with greater than 9 yrs data





## ============ Figure S5 - richness ===============

algae_cont<-com_sync_stab_long %>% 
  filter(metric %in% c("alg_sync","alg_rich","alg_stability"))  %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  mutate(mode = "algae") %>%
  rename(rich = alg_rich,
         sync = alg_sync,
         stab = alg_stability) %>%
  select(site,dataset,system,mode,rich,sync,stab)

fish_cont<-com_sync_stab_long %>% 
  filter(metric %in% c("fish_herb_sync","fish_herb_stability","fish_herb_rich"))  %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  mutate(mode = "fish_herbivore") %>%
  rename(rich = fish_herb_rich,
         sync = fish_herb_sync,
         stab = fish_herb_stability) %>% 
  select(site,dataset,system,mode,rich,sync,stab)

invert_cont<-com_sync_stab_long %>% 
  filter(metric %in% c("invert_herb_sync","invert_herb_stability","invert_herb_rich"))  %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  mutate(mode = "invert_herbivore") %>%
  rename(rich = invert_herb_rich,
         sync = invert_herb_sync,
         stab = invert_herb_stability) %>% 
  select(site,dataset,system,mode,rich,sync,stab)

rich_list<-list(algae_cont,fish_cont,invert_cont)
richness_df <- do.call(rbind, rich_list)
#richness_df<-rbind(algae_cont,herb_cont)

richness_df <- richness_df %>% mutate(mode = recode(mode, fish_herbivore = "fish herbivore",
                                                    invert_herbivore = "invert. herbivore"))


richness_df_upd<-richness_df %>% mutate(mode_color = case_when(system == "tropical" & str_detect(mode, "algae") ~ "green",
                                                               system == "temperate" & str_detect(mode, "algae")~ "darkgreen",
                                                               system == "tropical" &  str_detect(mode, "fish")~ "red",
                                                               system == "temperate" & str_detect(mode, "fish")~ "darkred",
                                                               system == "tropical" &  str_detect(mode, "invert")~ "blue",
                                                               system == "temperate" & str_detect(mode, "invert")~ "darkblue"))


richness_df_upd$sys_mode<-paste(richness_df_upd$system,richness_df_upd$mode, sep = "-")

  
# Z-score richness - homogeneity of variances 
library(DHARMa)
richness_df_upd_scaled<-richness_df_upd %>% 
  dplyr::group_by(mode,system) %>% 
  mutate(rich_scaled = scale(rich))



# try with LMER
sync_rich_lmer<-lmer(sync~rich_scaled *
                  system *
                  mode * 
                  (1|dataset),
                  richness_df_upd_scaled,
                REML = T)

summary(sync_rich_lmer)
ranova(sync_rich_lmer)
anova(sync_rich_lmer)

sync_rich_predict<-ggpredict(sync_rich_lmer, terms = c("rich_scaled [all]","mode","system"))
sync_rich_predict_df <- as.data.frame(sync_rich_predict)
plot(sync_rich_predict, add.data = T, limit.range = T)
names(sync_rich_predict_df) <- c("rich_scaled","sync","std.error","conf.low","conf.high","mode","system")
sync_rich_predict_df$sys_mode<-paste(sync_rich_predict_df$system,sync_rich_predict_df$mode, sep = "-")



# filter by range
rich_sync_range <- richness_df_upd_scaled %>%
  group_by(mode,system) %>%
  summarize(min_rich = min(rich_scaled,na.rm = T), max_rich = max(rich_scaled,na.rm = T))

sync_rich_predict_df_filt <- sync_rich_predict_df %>%
  inner_join(rich_sync_range, by = c("mode","system")) %>%
  filter(rich_scaled >= min_rich, rich_scaled <= max_rich)

dim(sync_rich_predict_df_filt)

# Extract results for table
sync_rich_table<-data.frame(anova(sync_rich_lmer))
sync_rich_table$variable <-rownames(sync_rich_table)
rownames(sync_rich_table) <-c("richness", "system", "trophic group (TG)", "rich.:system","rich.:TG","system:TG","rich:system:TG")
sync_rich_table_export<-sync_rich_table

# format p's correctly
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


# Add a column with subscript text using the variable
DenDF_values <- sapply(round(sync_rich_table$DenDF,2), function(x) HTML(paste0("<i>F</i><sub>", x, "</sub>")))
sync_rich_table$rounded_p_value <- sapply(sync_rich_table$Pr..F., round_p_value)
sync_rich_table$p_star <- makeStars(sync_rich_table$Pr..F.)
sync_rich_table$label <- paste0("</p>", DenDF_values, "  =  ",round(sync_rich_table$F.value,2),"<sup>", sync_rich_table$p_star, "</sup>")
sync_rich_table$variable <-c("richness", "system", "trophic group (TG)", "rich.:system","rich.:TG","system:TG","rich:system:TG")

# Calculate CIs and slopes
sync_rich_CI<-emtrends(sync_rich_lmer,pairwise~system|mode, var = "rich_scaled")$emtrends %>% 
  as.data.frame() %>%
  group_by(system,mode) %>%
  dplyr::mutate(facet = system,group_col = mode)


sync_rich_CI$y <-c(6.5,6.5,6.0,6.0,5.5,5.5)


# For stability and richness
stab_rich_lmer<-lmer(stab~rich_scaled *
                       system *
                       mode * 
                       (1|dataset),
                     richness_df_upd_scaled,
                     REML = T)

summary(stab_rich_lmer)
ranova(stab_rich_lmer)
anova(stab_rich_lmer)

stab_rich_predict<-ggpredict(stab_rich_lmer, terms = c("rich_scaled [all]","mode","system"))
plot(stab_rich_predict, add.data = T, limit.range = T)
stab_rich_predict_df <- as.data.frame(stab_rich_predict)
names(stab_rich_predict_df) <- c("rich_scaled","stab","std.error","conf.low","conf.high","mode","system")
stab_rich_predict_df$sys_mode<-paste(stab_rich_predict_df$system,stab_rich_predict_df$mode, sep = "-")
#stab_rich_predict$x

# manually extract coefficient estimates from ggpredict
# ignore warnings about perfect fit....they SHOULD be! These are the slopes

# filter by range
rich_stab_range <- richness_df_upd_scaled %>%
  group_by(mode,system) %>%
  summarize(min_rich = min(rich_scaled,na.rm = T), max_rich = max(rich_scaled,na.rm = T))

stab_rich_predict_df_filt <- stab_rich_predict_df %>%
  inner_join(rich_stab_range, by = c("mode","system")) %>%
  filter(rich_scaled >= min_rich, rich_scaled <= max_rich)


# Calculate CIs and slopes
stab_rich_CI<-emtrends(stab_rich_lmer,pairwise~system|mode, var = "rich_scaled")$emtrends %>% 
  as.data.frame() %>%
  group_by(system,mode) %>%
  dplyr::mutate(facet = system,group_col = mode)

stab_rich_CI$y <-c(24,24,22,22,20,20)

# Extract results for table
stab_rich_table<-data.frame(anova(stab_rich_lmer))
stab_rich_table$variable <-rownames(stab_rich_table)
rownames(stab_rich_table) <-c("richness", "system", "trophic group (TG)", "rich.:system","rich.:TG","system:TG","rich:system:TG")
stab_rich_table$variable <-c("richness", "system", "trophic group (TG)", "rich.:system","rich.:TG","system:TG","rich:system:TG")
stab_rich_table_export<-stab_rich_table

# Add a column with subscript text using the variable
DenDF_values_stab <- sapply(round(stab_rich_table$DenDF,2), function(x) HTML(paste0("<i>F</i><sub>", x, "</sub>")))
stab_rich_table$rounded_p_value <- sapply(stab_rich_table$Pr..F., round_p_value)
stab_rich_table$p_star <- makeStars(stab_rich_table$Pr..F.)
stab_rich_table$label <- paste0("</p>", DenDF_values_stab, "= ",round(stab_rich_table$F.value,2),"<sup>", stab_rich_table$p_star, "</sup>")


#sync_rich_table$estimate <-"synchrony"
#stab_rich_table$estimate <-"stability"


sync_rich_table_merge<-sync_rich_table %>% mutate(synchrony = label) %>% dplyr::select(variable,synchrony)
stab_rich_table_merge<-stab_rich_table %>% mutate(stability = label) %>% dplyr::select(variable,stability)

stab_rich_table_merge$id  <- 1:nrow(stab_rich_table_merge)

# finished table
rich_combined_table<-merge(sync_rich_table_merge,stab_rich_table_merge)
rich_combined_table<-rich_combined_table[order(rich_combined_table$id), ]

# remove rownames
row.names(rich_combined_table)<-NULL
rich_combined_table$id<-NULL

# complete table
#rich_table_done<-rich_combined_table %>% 
#  kbl(format = "html", escape = F, col.names = c("Fixed effects", "Synchrony","Stability")) %>%
#  kable_classic(html_font = "Times") %>%
#  footnote(general = "Results from linear mixed effects modeling (LME) with system, trophic group (TG), and richness as fixed effects. Dataset was treated as a random effect. Richness is scaled to satisfy assumptions of homogeneity of variances. Subscripts are denominator degrees of freedom from Satterhwaite approximation. Superscripts represent significance based on p-value: ns is p > 0.05; * p < 0.05; ** p < 0.01; *** p < 0.001; **** p < 0.0001",
#           escape = T) %>%
#  kable_styling() %>%
#  kable_save()



rich_combined_table_complete<-rich_combined_table %>% 
  #select(variable,label) %>%
  gt() %>% 
  fmt_markdown(columns = c(synchrony,stability)) %>%
  tab_source_note(source_note = md("Results from linear mixed effects modeling (LME) with system, trophic group (TG), and richness as fixed effects. Dataset was treated as a random effect. Richness is scaled to satisfy assumptions of homogeneity of variances. Subscripts are denominator degrees of freedom from Satterhwaite approximation. Superscripts represent significance based on p-value: ns is p > 0.05; * p < 0.05; ** p < 0.01; *** p < 0.001; **** p < 0.0001")) %>%
  gt::tab_options(table.font.names = "Times New Roman") %>%
  cols_label(
    variable = md("**Variable**"),
    synchrony = md("**Synchrony**"),
    stability ~ md("**Stability**")) #%>% 
  huxtable::to_rtf()
  #as_rtf()



# Export table as csvs and build manually

sync_rich_table_exp_merge<-sync_rich_table_export %>% mutate(group = "synchrony")
stab_rich_table_exp_merge<-stab_rich_table_export %>% mutate(group = "stability")
rich_table_exp<-rbind(sync_rich_table_exp_merge,stab_rich_table_exp_merge)
write.csv(rich_table_exp, "./Tables/rich_table.csv", row.names = F)



richness_point_sync<-richness_df_upd_scaled %>% mutate(x = rich_scaled,
                                                       predicted = sync,
                                                       facet = system,
                                                       group = mode)

rich_A<-plot(sync_rich_predict, limit.range = T) +
  theme_bw() +
  geom_point(data = richness_point_sync, aes(x = rich_scaled, y= predicted, color = group, fill = group), color = "black", pch =21, size =3, alpha = 0.6) +  #geom_point(color = "black",size = 3, pch =21, alpha = 0.6) +
  geom_line(linewidth = 1) +
  scale_fill_manual(values = c("#009E73","#CC70A9","#56B4E9")) +
  scale_color_manual(values = c("#009E73","#CC70A9","#56B4E9")) +
  removeGrid() +
  guides(color = "none", fill = "none") +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = 1,color = "black", linetype = "dashed") +
  #geom_vline(xintercept = 0,color = "black", linetype = "solid") +
  #geom_hline(yintercept = 0,color = "black", linetype = "solid") +
  labs(x = "Species richness (scaled)", y = "Monotrophic synchrony") +
  geom_text(data = sync_rich_CI, aes(label = paste("B = ",round(rich_scaled.trend   ,2), "; CI = [",round(lower.CL,2),", ",round(upper.CL,2),"]", sep = ""), 
                                         x = 1.2, y = y), hjust = 0, size = 2.5)

richness_point_stab<-richness_df_upd_scaled %>% mutate(x = rich_scaled,
                                                                   predicted = stab,
                                                                   facet = system,
                                                                   group = mode)

rich_B<-plot(stab_rich_predict,limit.range = T) +
  theme_bw() +
  geom_point(data = richness_point_stab, aes(x = rich_scaled, y= predicted, color = group, fill = group), color = "black", pch =21, size =3, alpha = 0.6) +  #geom_point(color = "black",size = 3, pch =21, alpha = 0.6) +
  geom_line(linewidth = 1) +
  scale_fill_manual(values = c("#009E73","#CC70A9","#56B4E9")) +
  scale_color_manual(values = c("#009E73","#CC70A9","#56B4E9")) +
  removeGrid() +
  guides(fill = "none") +
  theme(plot.title = element_blank(),
        strip.text.x = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank()) +
  #geom_hline(yintercept = 1,color = "black", linetype = "dashed") +
  #geom_vline(xintercept = 0,color = "black", linetype = "solid") +
  #geom_hline(yintercept = 0,color = "black", linetype = "solid") +
  labs(x = "Species richness (scaled)", y = "Monotrophic stability", color = "Trophic group") +
  geom_text(data = stab_rich_CI, aes(label = paste("B = ",round(rich_scaled.trend   ,2), "; CI = [",round(lower.CL,2),", ",round(upper.CL,2),"]", sep = ""), 
                                     x = 1.2, y = y), hjust = 0, size = 2.5)
# ADD MULTITROPHIC SYNCHRONY AND RICHNESS TO THIS!!!

richness_regression_lmer<-rich_A + rich_B + plot_layout(ncol = 1, guides = "collect")


ggsave("./Figures/Figure_S5_lmer.pdf",
       plot=richness_regression_lmer,
       width = 8.5,
       height = 6.5)

ggsave("./Figures/Figure_S5_lmer.jpeg",
       plot=richness_regression_lmer,
       width = 8.5,
       height = 6.5)

ggsave("./Figures/Figure_S5_lmer.TIFF",
       plot=richness_regression_lmer,
       width = 8.5,
       height = 6.5)

rich_lm_test<-richness_df_upd %>% 
  filter(mode == "fish herbivore", system == "tropical") %>%
  lm(sync ~ rich, data = .)




## ============ Figure 4 ===============
y_seq<-c(0.7,0.6,0.5,0.7,0.6,0.5) # for plotting text evenly across plots

## First row - panels A and B ###

stab_est_df<-stab_df_plot %>% filter(period == "Estimate") %>% mutate(stability = value) %>% select(-value)
sync_est_df<-sync_df_plot %>% filter(period == "Estimate") %>% mutate(sync = value) %>% select(-value)

sync_stab_est_df<-merge(stab_est_df,sync_est_df)
sync_stab_est_df$sys_mode<-paste(sync_stab_est_df$system,sync_stab_est_df$metric, sep = "-")

dim(sync_stab_est_df)

stab_effect_df<-stab_df_plot %>% filter(period == "Effect") %>% mutate(stability = value) %>% select(-value)
sync_effect_df<-sync_df_plot %>% filter(period == "Effect") %>% mutate(sync = value) %>% select(-value)

dim(stab_effect_df)
dim(sync_effect_df)

sync_stab_reg_df<-merge(stab_effect_df,sync_effect_df)
sync_stab_reg_df$sys_mode<-paste(sync_stab_reg_df$system,sync_stab_reg_df$metric, sep = "-")

dim(sync_stab_reg_df)


## Second row - panels C and D ###
estimates_merge<-sync_stab_est_df %>% 
  mutate(stability_est= stability,
         synchrony_est= sync) %>%
  select(-c(stability,sync,period))

effects_merge<-sync_stab_reg_df %>% 
  mutate(stability_eff= stability,
         synchrony_eff= sync) %>%
  select(-c(stability,sync,period))

est_eff_stab<-merge(estimates_merge,effects_merge)


# try regressing trophic level sync to multi
multitroph_stab<-sync_stab_est_df %>% 
  filter(metric == "multi") %>%
  mutate(stability_multi= stability,
         synchrony_multi= sync) %>%
  select(-c(stability,sync,period,metric,sys_mode,color,system,dataset))

mono_to_multi_df<-merge(estimates_merge, multitroph_stab, by = "site", all = TRUE)



# Figure 4 modified - LMER ####

# try with LMER
mono_to_multi_df<-merge(estimates_merge, multitroph_stab, by = "site", all = TRUE)

# Monotrophic synchrony v. monotrophic stability
sync_stab_est_df_nomulti<-sync_stab_est_df %>% 
  filter(!metric == "multi")

mono_v_mono_lmer <- lmer(stability~sync *
                           system *
                           metric *
                           (1|dataset),
                         sync_stab_est_df_nomulti,
                         REML = F)

summary(mono_v_mono_lmer)
ranova(mono_v_mono_lmer)
anova(mono_v_mono_lmer)


mono_v_mono_predict<-ggpredict(mono_v_mono_lmer, terms = c("sync [all]","metric","system"))
mono_v_mono_predict_df <- as.data.frame(mono_v_mono_predict)
test_plot<-plot(mono_v_mono_predict, add.data = T, limit.range = T)
names(mono_v_mono_predict_df) <- c("sync","stability","std.error","conf.low","conf.high","metric","system")
mono_v_mono_predict_df$sys_mode<-paste(mono_v_mono_predict_df$system,mono_v_mono_predict_df$metric, sep = "-")




# Calculate CIs and slopes
mono_v_mono_CI<-emtrends(mono_v_mono_lmer,pairwise~system|metric, var = "sync")$emtrends %>% 
  as.data.frame() %>%
  group_by(system,metric) %>%
  dplyr::mutate(facet = system,group_col = metric)

mono_v_mono_CI$y <-c(10,10,9.3,9.3,8.6,8.6)


# better limits
mono_v_mono<-plot(mono_v_mono_predict, limit.range = T) +
  theme_bw() +
  geom_line(linewidth = 1) +
  scale_fill_manual(values = c("#009E73","#CC70A9","#56B4E9")) +
  scale_color_manual(values = c("#009E73","#CC70A9","#56B4E9")) +
  removeGrid() +
  theme(plot.title = element_blank()) +
  geom_vline(xintercept = 1,color = "black", linetype = "dashed") +
  geom_vline(xintercept = 0,color = "black", linetype = "solid") +
  geom_hline(yintercept = 0,color = "black", linetype = "solid") +
  labs(x = "Monotrophic synchrony", y = "Monotrophic stability", color = "Trophic group") +
  geom_text(data = mono_v_mono_CI, aes(label = paste("B = ",round(sync.trend,2), "; CI = [",round(lower.CL,2),", ",round(upper.CL,2),"]", sep = ""), 
                                       x = 3, y = y), hjust = 0, size = 2.5) +
  coord_cartesian(ylim = c(0,10), xlim = c(0,5))



# Monotrophic synchrony v. multitrophic stability
mono_to_multi_df_nomulti<-mono_to_multi_df %>% 
  filter(!metric == "multi")

mono_v_multi_lmer <-lmer(stability_multi~synchrony_est *
                           system *
                           metric *
                           (1|dataset),
                         mono_to_multi_df_nomulti,
                         REML = F)

summary(mono_v_multi_lmer)
ranova(mono_v_multi_lmer)
anova(mono_v_multi_lmer)

mono_v_multi_predict<-ggpredict(mono_v_multi_lmer, terms = c("synchrony_est [all]","metric","system"))
mono_v_multi_predict_df <- as.data.frame(mono_v_multi_predict)
plot(mono_v_multi_predict, add.data = T, limit.range = T)
names(mono_v_multi_predict_df) <- c("synchrony_est","stability_multi","std.error","conf.low","conf.high","metric","system")
mono_v_multi_predict_df$sys_mode<-paste(mono_v_multi_predict_df$system,mono_v_multi_predict_df$metric, sep = "-")


# Calculate CIs and slopes
mono_v_multi_CI<-emtrends(mono_v_multi_lmer,pairwise~system|metric, var = "synchrony_est")$emtrends %>% 
  as.data.frame() %>%
  group_by(system,metric) %>%
  dplyr::mutate(facet = system,group_col = metric)

mono_v_multi_CI$y <-c(10,10,9.3,9.3,8.6,8.6)



# better limits
mono_v_multi<-plot(mono_v_multi_predict, limit.range = T) +
  theme_bw() +
  geom_line(linewidth = 1) +
  scale_fill_manual(values = c("#009E73","#CC70A9","#56B4E9")) +
  scale_color_manual(values = c("#009E73","#CC70A9","#56B4E9")) +
  removeGrid() +
  theme(plot.title = element_blank()) +
  geom_vline(xintercept = 1,color = "black", linetype = "dashed") +
  geom_vline(xintercept = 0,color = "black", linetype = "solid") +
  geom_hline(yintercept = 0,color = "black", linetype = "solid") +
  labs(x = "Monotrophic synchrony", y = "Multitrophic stability", color = "Trophic group") +
  geom_text(data = mono_v_multi_CI, aes(label = paste("B = ",round(synchrony_est.trend,2), "; CI = [",round(lower.CL,2),", ",round(upper.CL,2),"]", sep = ""), 
                                        x = 3, y = y), hjust = 0, size = 2.5) +
  coord_cartesian(ylim = c(0,10), xlim = c(0,5))


# Monotrophic synchrony v. multitrophic stability

multi_v_multi_lmer<-lmer(stability~sync *
                           system *
                           metric *
                           (1|dataset),
                         sync_stab_reg_df,
                         REML = F)

summary(multi_v_multi_lmer)
ranova(multi_v_multi_lmer)
anova(multi_v_multi_lmer)

multi_v_multi_predict<-ggpredict(multi_v_multi_lmer, terms = c("sync [all]","metric","system"))
multi_v_multi_predict_df <- as.data.frame(multi_v_multi_predict)
plot(multi_v_multi_predict, add.data = T, limit.range = T)
names(multi_v_multi_predict_df) <- c("sync","stability","std.error","conf.low","conf.high","metric","system")
multi_v_multi_predict_df$sys_mode<-paste(multi_v_multi_predict_df$system,multi_v_multi_predict_df$metric, sep = "-")

# Calculate CIs and slopes
multi_v_multi_CI<-emtrends(multi_v_multi_lmer,pairwise~system|metric, var = "sync")$emtrends %>% 
  as.data.frame() %>%
  group_by(system,metric) %>%
  dplyr::mutate(facet = system,group_col = metric)


multi_v_multi_CI$y <-c(10,10,9,9,8,8)




# better limits
multi_v_multi<-plot(multi_v_multi_predict, limit.range = T) +
  theme_bw() +
  geom_line(linewidth = 1) +
  scale_fill_manual(values = c("#009E73","#CC70A9","#56B4E9")) +
  scale_color_manual(values = c("#009E73","#CC70A9","#56B4E9")) +
  removeGrid() +
  theme(plot.title = element_blank()) +
  geom_vline(xintercept = 1,color = "black", linetype = "dashed") +
  geom_vline(xintercept = 0,color = "black", linetype = "solid") +
  geom_hline(yintercept = 0,color = "black", linetype = "solid") +
  labs(x = "Synchrony contribution", y = "Multitrophic stability", color = "Trophic group") +
  geom_text(data = multi_v_multi_CI, aes(label = paste("B = ",round(sync.trend,2), "; CI = [",round(lower.CL,2),", ",round(upper.CL,2),"]", sep = ""), 
                                         x = 2, y = y), hjust = 0, size = 2.5) +
  coord_cartesian(ylim = c(-5,10), xlim = c(-4,6))


mono_v_mono_tag<-tag_facet_Wtop(mono_v_mono, tag_pool = c("A","B"))

mono_v_multi_tag<-tag_facet_NP(mono_v_multi, tag_pool = c("C","D"))

multi_v_multi_tag<-tag_facet_NP(multi_v_multi, tag_pool = c("E","F"))


figure_4_plot_lmer<-mono_v_mono_tag + mono_v_multi_tag + multi_v_multi_tag + plot_layout(ncol = 1, guides = "collect") & theme(legend.position = "bottom")

ggsave("./Figures/Figure_4_lmer.pdf",
       plot=figure_4_plot_lmer,
       width = 8,
       height = 10)

ggsave("./Figures/Figure_4_lmer.TIFF",
       plot=figure_4_plot_lmer,
       width = 8,
       height = 10)

# Extract results for table
# (1) Mono v mono
mono_v_mono_table<-data.frame(anova(mono_v_mono_lmer))
mono_v_mono_table$variable <-rownames(mono_v_mono_table)
rownames(mono_v_mono_table) <-c("sync", "system", "trophic group (TG)", "sync:system","sync:TG","system:TG","sync:system:TG")
mono_v_mono_table$variable <-c("sync", "system", "trophic group (TG)", "sync:system","sync:TG","system:TG","sync:system:TG")
mono_v_mono_table_export<-mono_v_mono_table

# Add a column with subscript text using the variable
DenDF_values_stab_1 <- sapply(round(mono_v_mono_table$DenDF,2), function(x) HTML(paste0("<i>F</i><sub>", x, "</sub>")))
mono_v_mono_table$rounded_p_value <- sapply(mono_v_mono_table$Pr..F., round_p_value)
mono_v_mono_table$p_star <- makeStars(mono_v_mono_table$Pr..F.)
#stab_rich_table$label <- paste0("</p>", DenDF_values, "</p>", "  =  ",round(stab_rich_table$F.value,2),stab_rich_table$rounded_p_value)
mono_v_mono_table$label <- paste0("</p>", DenDF_values_stab_1, "</p>", "  =  ",round(mono_v_mono_table$F.value,2),"<sup>", mono_v_mono_table$p_star, "</sup>")
mono_v_mono_table$mono_V_mono <-mono_v_mono_table$label


# (2) mono v. multi
mono_v_multi_table<-data.frame(anova(mono_v_multi_lmer))
mono_v_multi_table$variable <-rownames(mono_v_multi_table)
rownames(mono_v_multi_table) <-c("sync", "system", "trophic group (TG)", "sync:system","sync:TG","system:TG","sync:system:TG")
mono_v_multi_table$variable <-c("sync", "system", "trophic group (TG)", "sync:system","sync:TG","system:TG","sync:system:TG")
mono_v_multi_table_export<-mono_v_multi_table

# Add a column with subscript text using the variable
DenDF_values_stab_2 <- sapply(round(mono_v_multi_table$DenDF,2), function(x) HTML(paste0("<i>F</i><sub>", x, "</sub>")))
mono_v_multi_table$rounded_p_value <- sapply(mono_v_multi_table$Pr..F., round_p_value)
mono_v_multi_table$p_star <- makeStars(mono_v_multi_table$Pr..F.)
#stab_rich_table$label <- paste0("</p>", DenDF_values, "</p>", "  =  ",round(stab_rich_table$F.value,2),stab_rich_table$rounded_p_value)
mono_v_multi_table$label <- paste0("</p>", DenDF_values_stab_2, "</p>", "  =  ",round(mono_v_multi_table$F.value,2),"<sup>", mono_v_multi_table$p_star, "</sup>")
mono_v_multi_table$mono_V_multi <-mono_v_multi_table$label


# (3) multi v. multi
multi_v_multi_table<-data.frame(anova(multi_v_multi_lmer))
multi_v_multi_table$variable <-rownames(multi_v_multi_table)
rownames(multi_v_multi_table) <-c("sync", "system", "trophic group (TG)", "sync:system","sync:TG","system:TG","sync:system:TG")
multi_v_multi_table$variable <-c("sync", "system", "trophic group (TG)", "sync:system","sync:TG","system:TG","sync:system:TG")
multi_v_multi_table_export<-multi_v_multi_table

# Add a column with subscript text using the variable
DenDF_values_stab_3 <- sapply(round(multi_v_multi_table$DenDF,2), function(x) HTML(paste0("<i>F</i><sub>", x, "</sub>")))
multi_v_multi_table$rounded_p_value <- sapply(multi_v_multi_table$Pr..F., round_p_value)
multi_v_multi_table$p_star <- makeStars(multi_v_multi_table$Pr..F.)
#stab_rich_table$label <- paste0("</p>", DenDF_values, "</p>", "  =  ",round(stab_rich_table$F.value,2),stab_rich_table$rounded_p_value)
multi_v_multi_table$label <- paste0("</p>", DenDF_values_stab_3, "</p>", "  =  ",round(multi_v_multi_table$F.value,2),"<sup>", multi_v_multi_table$p_star, "</sup>")
multi_v_multi_table$multi_V_multi <-multi_v_multi_table$label

# bring together

mono_v_mono_merge<-mono_v_mono_table %>% dplyr::select(variable,mono_V_mono)
mono_v_multi_merge<-mono_v_multi_table %>% dplyr::select(variable,mono_V_multi)
multi_v_multi_merge<-multi_v_multi_table %>% dplyr::select(variable,multi_V_multi)

multi_v_multi_merge$id  <- 1:nrow(multi_v_multi_merge)

# finished table
sync_table_list <- list(mono_v_mono_merge,mono_v_multi_merge,multi_v_multi_merge)
sync_combined_table<-Reduce(function(x, y) merge(x, y, all=TRUE), sync_table_list, accumulate=FALSE)

sync_combined_table<-sync_combined_table[order(sync_combined_table$id), ]

# remove rownames
row.names(sync_combined_table)<-NULL
sync_combined_table$id<-NULL

# complete table
#sync_table_done<-
sync_combined_table %>% 
  kbl(format = "html", escape = F, col.names = c("Fixed effects", "Mono. to mono.","Mono. to multi.","Multi. to multi.")) %>%
  kable_classic(html_font = "Times") %>%
  footnote(general = "Results from linear mixed effects modeling (LME) with system, trophic group (TG), and synchrony as fixed effects. Dataset was treated as a random effect. Subscripts are denominator degrees of freedom from Satterhwaite approximation. Superscripts represent significance based on p-value: ns is p > 0.05; * p < 0.05; ** p < 0.01; *** p < 0.001; **** p < 0.0001",
           escape = T) %>%
  kable_styling() 





## ============ Visual abstract - part 2 ===============

stab_effect_df<-stab_df_plot %>% filter(period == "Effect") %>% mutate(stability = value) %>% select(-value)
sync_effect_df<-sync_df_plot %>% filter(period == "Effect") %>% mutate(sync = value) %>% select(-value)

dim(stab_effect_df)
dim(sync_effect_df)

sync_stab_reg_df<-merge(stab_effect_df,sync_effect_df)
dim(sync_stab_reg_df)




# Use LMER results
# extract Betas and plot those with CIs




# plot it -- original
lm_dotplot<-mono_v_mono_CI %>% 
  #filter(grepl("sync", coef)) %>%
  mutate(system = recode(system, temperate = "on temperate multitrophic stability", tropical = "on tropical multitrophic stability")) %>%
  ggplot(aes(x = sync.trend, y = "", fill = metric), color = metric) +
  geom_vline(xintercept = 0) +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL, color = metric), width = 0.01,position=position_dodge(width = -1)) +
  geom_point(size = 6, pch =21, color = "black",position=position_dodge(width = -1)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.15,0.8),
        axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA)) +
  facet_wrap(~system, ncol = 1) +
  removeGrid() +
  labs(y = "", x= "Estimate") +
  scale_color_manual(
    values = c(invert="#56B4E9",
                fish = "#CC70A9",
                algae = "#009E73"),
    labels = c("algae", "herb. fish", "herb. invert")) +
  scale_fill_manual(
    values = c(algae = "#009E73",
               fish = "#CC70A9",
               invert="#56B4E9"),na.value="transparent",
    labels = c("algae", "herb. fish", "herb. invert")) +
  guides(color = "none")

ggsave("./Figures/vis_abstract/dot_plot.pdf",
       plot=lm_dotplot,
       width = 4,
       height = 2.8)



# Two panel dotplot
# dotplot part A

#mono_dotplot$metric <- factor(mono_dotplot$metric,levels = c("invert","fish","algae"))

group_A<-mono_v_mono_CI %>% 
  #filter(grepl("sync", coef)) %>%
  #mutate(system = recode(system, temperate = "on temperate multitrophic stability", tropical = "on tropical multitrophic stability")) %>%
  #ggplot(aes(x = Estimate, y = "", fill = if_else(p_val< 0.05, metric,NA), group = metric, color = metric)) +
  ggplot(aes(x = sync.trend, y = metric,  fill = metric, color = metric)) +
  geom_vline(xintercept = 0) +
  #geom_vline(xintercept = 1, linetype = "dashed") +
  #geom_hline(yintercept = 0) +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0.01,position=position_dodge(width = -1)) +
  geom_point(size = 4, pch =21,position=position_dodge(width = -1)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        #legend.position = c(0.15,0.8),
        legend.position = "none",
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        #legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_blank(), 
        strip.background = element_blank()) +
  removeGrid() +
  labs(y = "", x= "Estimate") +
  scale_color_manual(
    values = c(invert="#56B4E9",
               fish = "#CC70A9",
               algae = "#009E73"),
    labels = c("herb. invert","herb. fish","algae")) +
  scale_fill_manual(
    values = c(algae = "#009E73",
               fish = "#CC70A9",
               invert="#56B4E9"),
    na.value="white",labels = c("herb. invert","herb. fish","algae","ns")) +
  guides(color = "none",
         fill = "none") +
  xlim(-2.5,1.2) +
  facet_wrap(~system, ncol = 1) #+
  guides(y = ggh4x::guide_axis_nested(delim = "&", angle = 90))


combined_results$metric <- factor(combined_results$metric,levels = c("invert","fish","algae"))

group_B<-multi_v_multi_CI %>% 
  #filter(grepl("sync", coef)) %>%
  ggplot(aes(x = sync.trend, y = metric,  fill = metric, color = metric)) +
  geom_vline(xintercept = 0) +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0.01,position=position_dodge(width = -1)) +
  geom_point(size = 4, pch =21,position=position_dodge(width = -1)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA_character_), # necessary to avoid drawing plot outline
        axis.text.y = element_blank(),
        strip.text = element_blank(), 
        strip.background = element_blank()) +
  facet_wrap(~system, ncol = 1) +
  removeGrid() +
  labs(y = "", x= "Estimate") +
  scale_color_manual(
    values = c(algae = "#009E73",
               fish = "#CC70A9",
               invert="#56B4E9"),
    labels = c("algae","herb. fish","herb. invert")) +
  scale_fill_manual(
    values = c(algae = "#009E73",
               fish = "#CC70A9",
               invert="#56B4E9"),na.value="white",
    labels = c("algae","herb. fish","herb. invert")) +
  guides(color = "none",
         fill = guide_legend(override.aes = list(fill = c("#009E73", "#CC70A9", "#56B4E9")),
                              #title = "Category",
                              values = c("#009E73", "#CC70A9", "#56B4E9", "white"),
                              labels = c("herb. invert","herb. fish","algae", "ns"))) +
  xlim(-2.5,1.2)
  
group_A_titled<-group_A + ggtitle("monotrophic synchrony on multitrophic stability") + 
  theme(plot.title = element_text(size = 12))

group_B_titled<-group_B + ggtitle("synchrony contribution on multitrophic stability") +
  theme(plot.title = element_text(size = 12))


modified_dotplot<-group_A_titled / group_B_titled + plot_layout(ncol = 1, guides = "collect") & theme(legend.position = "bottom")

ggsave("./Figures/vis_abstract/modified_dotplot_A.pdf",
       plot=group_A,
       width = 4,
       height = 1.5)

ggsave("./Figures/vis_abstract/modified_dotplot_B.pdf",
       plot=group_B,
       width = 4,
       height = 1.5)


ggsave("./Figures/vis_abstract/modified_dotplot.pdf",
       plot=modified_dotplot,
       width = 4,
       height = 4)







## Supplementary figures ####

stab_df_plot %>% 
  group_by(site,system) %>% 
  summarise_all(mean) %>%
  group_by(system) %>%
  summarize(site_count = n())



length(unique(sync_df_plot$site))

### Figure S4 - Comparison across datasets ####
# Rename the IMOS/RLS/ATRC dataset
stab_df_plot<-mutate(stab_df_plot, dataset = recode(dataset, IMOS = "RLS/ATRC"))
sync_df_plot<-mutate(sync_df_plot, dataset = recode(dataset, IMOS = "RLS/ATRC"))

sync_df_plot$dataset <- factor(sync_df_plot$dataset, levels = c("RLS/ATRC", "PISCO", "MCR", "GBR-LTMP", "TCRMP"))

sync_together_dataset<-sync_df_plot %>% 
  ggplot(aes(x = dataset, y = value, fill = metric)) +
  geom_hline(yintercept = 0) +
  #geom_sina(size = 2, pch =21, alpha = 0.4) +
  geom_sina(size = 2, pch =21, alpha = 0.4,position = position_dodge(width = .75)) +
  geom_boxplot(alpha = 0.5,outlier.shape = NA, notch = T) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_fill_manual(
    values = c(algae="#009E73",
               fish = "#CC70A9",
               invert = "#56B4E9",
               multi = "darkgrey")) +
  labs(x = "Organizational scale", y = "Variance-ratio") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    axis.text.x = element_blank()) +
  scale_y_continuous(limits = c(NA, NA),breaks = c(-3,-1,0,1,3,5)) +
  facet_wrap(~period, nrow = 1, scales = "free") +
  removeGrid()

sync_together_dataset<-tag_facet_Wtop(sync_together_dataset, tag_pool = c("A","B"))

stab_df_plot$dataset <- factor(stab_df_plot$dataset, levels = c("RLS/ATRC", "PISCO", "MCR", "GBR-LTMP", "TCRMP"))

unique(stab_df_plot$dataset)
stab_together_dataset<-stab_df_plot %>% 
  mutate(dataset = recode(dataset, "GBR-LTMP" = "AIMS-LTMP")) %>%
  ggplot(aes(interaction(dataset, system, sep = "&"), y = value, fill = metric)) +
  geom_hline(yintercept = 0) +
  geom_sina(size = 2, pch =21, alpha = 0.4,position = position_dodge(width = .75)) +
  geom_boxplot(alpha = 0.5,outlier.shape = NA, notch = T) +
  scale_fill_manual(
    values = c(algae="#009E73",
               fish = "#CC70A9",
               invert = "#56B4E9",
               multi = "darkgrey"
    )) +
  labs(x = "Dataset", y = "Stability") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()) +
  scale_y_continuous(limits = c(NA, NA),breaks = c(-4,-2,0,2,5,10,15,20,25)) +
  facet_wrap(~period, scales = "free") +
  removeGrid() +
  guides(x = ggh4x::guide_axis_nested(delim = "&"))

stab_together_dataset<-tag_facet_Wtop(stab_together_dataset, tag_pool = c("C","D"))


system_plot_dataset<-sync_together_dataset + stab_together_dataset + plot_layout(nrow = 2,guides = "collect") & theme(legend.position = "bottom")

ggsave("./Figures/Figure_S4.pdf",
       plot=system_plot_dataset,
       width = 9,
       height = 7)

ggsave("./Figures/Figure_S4.TIFF",
       plot=system_plot_dataset,
       width = 9,
       height = 7)




# LMER ####

# LMER for synchrony
# by system
# by dataset
sync_estimates<-sync_df_plot %>% filter(period == "Estimate")
sync_effect<-sync_df_plot %>% filter(period == "Effect")

# Estimates
sync_est_lmer<-lmer(value ~
                      system * 
                      metric *
                      (1|dataset),
                    data = sync_estimates,
                    REML = F)


summary(sync_est_lmer)
ranova(sync_est_lmer)
anova(sync_est_lmer)
plot(ggpredict(sync_est_lmer, terms = c("metric","system")))

# Effect
sync_eff_lmer<-lmer(value ~
                      system * 
                      metric *
                      (1|dataset),
                    data = sync_effect,
                    REML = F)


summary(sync_eff_lmer)
ranova(sync_eff_lmer)
anova(sync_eff_lmer)
plot(ggpredict(sync_eff_lmer, terms = c("metric","system")))

#sync_eff_temp_aov<-aov(value ~ system* metric * dataset,sync_effect)
#plot(ggpredict(sync_eff_temp_aov, terms = c("metric","system","dataset")))


# LMER for stability
# by system
# by dataset

stab_estimates<-stab_df_plot %>% filter(period == "Estimate")
stab_effect<-stab_df_plot %>% filter(period == "Effect")

# Estimates
stab_est_lmer<-lmer(value ~
                      system * 
                      metric *
                      (1|dataset),
                    data = stab_estimates,
                    REML = F)


summary(stab_est_lmer)
ranova(stab_est_lmer)
anova(stab_est_lmer)
ranova(stab_est_lmer)

plot(ggpredict(stab_est_lmer, terms = c("metric","system")))

# Effect
stab_eff_lmer<-lmer(value ~
                      system * 
                      metric *
                      (1|dataset),
                    data = stab_effect,
                    REML = F)


summary(stab_eff_lmer)
ranova(stab_eff_lmer)
anova(stab_eff_lmer)
plot(ggpredict(stab_eff_lmer, terms = c("metric","system")))


stab_eff_aov<-aov(value ~ system * metric, stab_effect)
summary(stab_eff_aov)
plot(ggpredict(stab_eff_lmer, terms = c("metric","system")))





# test within system
# tropical
sync_effect_trop<-sync_effect %>% filter(system == "tropical")
sync_effect_trop_aov<-aov(value ~ metric * dataset,sync_effect_trop)
summary(sync_effect_trop_aov) # differences across datasets -- but not interaction between metrics
plot(ggpredict(sync_effect_trop_aov, terms = c("metric","dataset")))

stab_effect_trop<-stab_effect %>% filter(system == "tropical")
stab_eff_trop_aov<-aov(value ~ metric * dataset,stab_effect_trop)
summary(stab_eff_trop_aov) # differences across datasets -- but not interaction between metrics
plot(ggpredict(stab_eff_trop_aov, terms = c("metric","dataset")))





#temperate
sync_effect_temp<-sync_effect %>% filter(system == "temperate")
sync_eff_temp_aov<-aov(value ~ metric * dataset,sync_effect_temp)
summary(sync_eff_temp_aov) # all significant
plot(ggpredict(sync_eff_temp_aov, terms = c("metric","dataset")))

stab_effect_temp<-stab_effect %>% filter(system == "temperate")
stab_eff_temp_aov<-aov(value ~ metric * dataset,stab_effect_temp)
summary(stab_eff_temp_aov) # all significant
plot(ggpredict(stab_eff_temp_aov, terms = c("metric","dataset")))



# Summary tables ####
sync_estimates %>% group_by(system,metric) %>% summarize(mean_sync = mean(value,na.rm = T),
                                                         sd_sync = sd(value, na.rm = T))

sync_effect %>% group_by(system,metric) %>% summarize(mean_sync = mean(value,na.rm = T),
                                                      sd_sync = sd(value, na.rm = T))

sync_estimates %>% 
  filter(metric == "multi") %>%
  group_by(system) %>% 
  summarize(mean_sync = mean(value,na.rm = T),
            sd_sync = sd(value, na.rm = T))



stab_estimates %>% group_by(system,metric) %>% summarize(mean_stab = mean(value,na.rm = T),
                                                         sd_stab = sd(value, na.rm = T))


stab_estimates %>% group_by(dataset) %>% summarize(site_count = n_distinct(site))



# Extract results for table
# (1) Synchrony estimates
sync_est_lmer_table<-data.frame(anova(sync_est_lmer))
sync_est_lmer_table$variable <-rownames(sync_est_lmer_table)
rownames(sync_est_lmer_table) <-c("system", "trophic group (TG)","system:TG")
sync_est_lmer_table$variable <-c("system", "trophic group (TG)","system:TG")
sync_est_table_export<-sync_est_lmer_table

# Add a column with subscript text using the variable
DenDF_sync_est <- sapply(round(sync_est_lmer_table$DenDF,2), function(x) HTML(paste0("<i>F</i><sub>", x, "</sub>")))
sync_est_lmer_table$rounded_p_value <- sapply(sync_est_lmer_table$Pr..F., round_p_value)
sync_est_lmer_table$p_star <- makeStars(sync_est_lmer_table$Pr..F.)
#stab_rich_table$label <- paste0("</p>", DenDF_values, "</p>", "  =  ",round(stab_rich_table$F.value,2),stab_rich_table$rounded_p_value)
sync_est_lmer_table$label <- paste0("</p>", DenDF_sync_est, "</p>", "  =  ",round(sync_est_lmer_table$F.value,2),"<sup>", sync_est_lmer_table$p_star, "</sup>")
sync_est_lmer_table$sync_est <-sync_est_lmer_table$label


# (2) Synchrony effect
sync_eff_lmer_table<-data.frame(anova(sync_eff_lmer))
sync_eff_lmer_table$variable <-rownames(sync_eff_lmer_table)
rownames(sync_eff_lmer_table) <-c("system", "trophic group (TG)","system:TG")
sync_eff_lmer_table$variable <-c("system", "trophic group (TG)","system:TG")
sync_eff_lmer_table_export<-sync_eff_lmer_table

# Add a column with subscript text using the variable
DenDF_sync_eff <- sapply(round(sync_eff_lmer_table$DenDF,2), function(x) HTML(paste0("<i>F</i><sub>", x, "</sub>")))
sync_eff_lmer_table$rounded_p_value <- sapply(sync_eff_lmer_table$Pr..F., round_p_value)
sync_eff_lmer_table$p_star <- makeStars(sync_eff_lmer_table$Pr..F.)
#stab_rich_table$label <- paste0("</p>", DenDF_values, "</p>", "  =  ",round(stab_rich_table$F.value,2),stab_rich_table$rounded_p_value)
sync_eff_lmer_table$label <- paste0("</p>", DenDF_sync_eff, "</p>", "  =  ",round(sync_eff_lmer_table$F.value,2),"<sup>", sync_eff_lmer_table$p_star, "</sup>")
sync_eff_lmer_table$sync_eff <-sync_eff_lmer_table$label


# (3) Stability estimate
stab_est_lmer_table<-data.frame(anova(stab_est_lmer))
stab_est_lmer_table$variable <-rownames(stab_est_lmer_table)
rownames(stab_est_lmer_table) <-c("system", "trophic group (TG)","system:TG")
stab_est_lmer_table$variable <-c("system", "trophic group (TG)","system:TG")
stab_est_lmer_table_export<-stab_est_lmer_table

# Add a column with subscript text using the variable
DenDF_stab_est <- sapply(round(stab_est_lmer_table$DenDF,2), function(x) HTML(paste0("<i>F</i><sub>", x, "</sub>")))
stab_est_lmer_table$rounded_p_value <- sapply(stab_est_lmer_table$Pr..F., round_p_value)
stab_est_lmer_table$p_star <- makeStars(stab_est_lmer_table$Pr..F.)
#stab_rich_table$label <- paste0("</p>", DenDF_values, "</p>", "  =  ",round(stab_rich_table$F.value,2),stab_rich_table$rounded_p_value)
stab_est_lmer_table$label <- paste0("</p>", DenDF_stab_est, "</p>", "  =  ",round(stab_est_lmer_table$F.value,2),"<sup>", stab_est_lmer_table$p_star, "</sup>")
stab_est_lmer_table$stab_est <-stab_est_lmer_table$label


# (4) Stability effect
stab_eff_lmer_table<-data.frame(anova(stab_eff_lmer))
stab_eff_lmer_table$variable <-rownames(stab_eff_lmer_table)
rownames(stab_eff_lmer_table) <-c("system", "trophic group (TG)","system:TG")
stab_eff_lmer_table$variable <-c("system", "trophic group (TG)","system:TG")
stab_eff_lmer_table_export<-stab_eff_lmer_table

# Add a column with subscript text using the variable
DenDF_stab_eff <- sapply(round(stab_eff_lmer_table$DenDF,2), function(x) HTML(paste0("<i>F</i><sub>", x, "</sub>")))
stab_eff_lmer_table$rounded_p_value <- sapply(stab_eff_lmer_table$Pr..F., round_p_value)
stab_eff_lmer_table$p_star <- makeStars(stab_eff_lmer_table$Pr..F.)
#stab_rich_table$label <- paste0("</p>", DenDF_values, "</p>", "  =  ",round(stab_rich_table$F.value,2),stab_rich_table$rounded_p_value)
stab_eff_lmer_table$label <- paste0("</p>", DenDF_stab_eff, "</p>", "  =  ",round(stab_eff_lmer_table$F.value,2),"<sup>", stab_eff_lmer_table$p_star, "</sup>")
stab_eff_lmer_table$stab_eff <-stab_eff_lmer_table$label

# bring together

sync_est_merge<-sync_est_lmer_table %>% dplyr::select(variable,sync_est)
sync_eff_merge<-sync_eff_lmer_table %>% dplyr::select(variable,sync_eff)
stab_est_merge<-stab_est_lmer_table %>% dplyr::select(variable,stab_est)
stab_eff_merge<-stab_eff_lmer_table %>% dplyr::select(variable,stab_eff)

sync_est_merge$id  <- 1:nrow(sync_est_merge)

# finished table
eff_est_table_list <- list(sync_est_merge,sync_eff_merge,stab_est_merge,stab_eff_merge)
eff_est_table<-Reduce(function(x, y) merge(x, y, all=TRUE), eff_est_table_list, accumulate=FALSE)

eff_est_table<-eff_est_table[order(eff_est_table$id), ]

# remove rownames
row.names(eff_est_table)<-NULL
eff_est_table$id<-NULL

# complete table
#sync_table_done<-
eff_est_table %>% 
  kbl(format = "html", escape = F, col.names = c("Fixed effects", "Synchrony estimate", "Synchrony contr.", "Stability estimate", "Stability contr.")) %>%
  kable_classic(html_font = "Times") %>%
  footnote(general = "Results from linear mixed effects modeling (LME) with system and trophic group (TG) as fixed effects. Dataset was treated as a random effect. Subscripts are denominator degrees of freedom from Satterhwaite approximation. Superscripts represent significance based on p-value: ns is p > 0.05; * p < 0.05; ** p < 0.01; *** p < 0.001; **** p < 0.0001",
           escape = T) %>%
  kable_styling() 




# Markdown


rmarkdown::render("./Tables/tables.Rmd")
#rmarkdown::render("asynchrony_synthesis_code.Rmd")

#rmarkdown::render("asynchrony_synthesis_code_synth_curation.Rmd")
#knit('asynchrony_synthesis_code_synth_curation.Rmd')

# END # 

