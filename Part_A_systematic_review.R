# Code in support of "Understanding diversity-synchrony-stability relationships in multitrophic communities"
# in Nature Ecology & Evolution 2024
# Griffin Srednick and Stephen Swearer
# DOI: XXX


# ====== Part A - Systematic Review ======


# Packages
library(tidyverse)
library(patchwork)
library(ggpmisc)
library(ggExtra)
library(ggrepel)





#  Load Data ####
trophic_lit<-read.csv("./Data/curated_literature_data_anon.csv")

names(trophic_lit)

  
# PART A - Generate figures for systematic review ####


## Figure 1 ####
# Elements built here -- assembled in Adobe Illustrator

trophic_lit$pop_com_both <- tolower(trophic_lit$pop_com_both)
unique(trophic_lit$pop_com_both)

trophic_lit$pop_com_both <- factor(trophic_lit$pop_com_both, levels=c("community","population","predator-prey","plant-pollinator","host-parasite"))


# Figure 1A - Monotrophic vs. Multitrophic studies
P1_bar<-trophic_lit %>%
  group_by(trophic_study) %>%
  count() %>%
  mutate(ypos = cumsum(n)- 0.5*n) %>%
  ggplot(aes(x="", y=n, fill=trophic_study)) +
  geom_bar(stat="identity", width=1, color="white") +
  theme_void() +
  geom_text(aes(y = ypos, 
                label = n), size=5) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#F8766D","#6BAED6")) 



ggsave("./Figures/Figure_1/P1_bar.pdf",
       plot=P1_bar,
       width = 1.5,
       height = 4)


# Figure 1B - Monotrophic organizational scale  

P2_A_df<-trophic_lit %>%
  filter(trophic_study == "monotrophic") %>%
  mutate(system_type = recode_factor(system_type, 
                              AT = "MAT",
                              MAT = "MAT")) %>%
  group_by(pop_com_both,system_type) %>%
  count() 

P2_A <- P2_A_df %>%
  mutate(ypos = cumsum(n)- 0.5*n) %>%
  ggplot(aes(x=as.character(system_type), y=n, fill=pop_com_both, group = pop_com_both)) +
  geom_bar(stat="identity", width=1, position = 'stack',color="white") +
  scale_fill_manual(values = c("#440154FF","#3B528BFF")) +
  geom_text(aes(x=system_type,label=n),
            position = position_stack(vjust = .5), color = "white") +
  theme_bw() +
  removeGrid() +
  labs(y = "no. papers") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")


ggsave("./Figures/Figure_1/P2_A.pdf",
       plot=P2_A,
       width = 6,
       height = 4)

# Figure 1B - Multitrophic organizational scale  
P2_B_df<-trophic_lit %>%
  filter(!trophic_study == "monotrophic") %>%
  mutate(system_type = recode_factor(system_type, 
                              AT = "MAT",
                              MAT = "MAT")) %>%
  group_by(pop_com_both,system_type) %>%
  count() 

P2_B<-P2_B_df %>%
  #mutate(ypos = cumsum(n)- 0.5*n) %>%
  ggplot(aes(x=as.character(system_type), y=n, fill=pop_com_both, group = pop_com_both)) +
  geom_bar(stat="identity", width=1, position = 'stack',color="white") +
  scale_fill_viridis_d() +
  labs(y = "no. papers",
       x = "system type") +
  geom_text(aes(x=system_type,label=n),
            position = position_stack(vjust = .5), color = "white") +
  theme_bw() +
  #scale_fill_manual(values = c("#440154FF","#3B528BFF")) +
  theme(legend.position=c(0.4,0.85),
        legend.title = element_blank()) +
  removeGrid() +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))


g <- ggplot_build(P2_B)
unique(g$data[[1]]["fill"])


ggsave("./Figures/Figure_1/P2_B.pdf",
       plot=P2_B,
       width = 6,
       height = 4)

# Combine
P2<-P2_A + P2_B + plot_layout(nrow = 2) 

ggsave("./Figures/Figure_1/P2.pdf",
       plot=P2,
       width = 5,
       height = 8)

pie_colors<-colorRampPalette(brewer.pal(9,"Blues"))(7)


# Format column
trophic_lit$asynch_is_stab_at_com<-tolower(trophic_lit$asynch_is_stab_at_com)

trophic_lit_upd<-trophic_lit %>% 
  mutate(asynch_is_stab_at_com = case_when(
    str_detect(asynch_is_stab_at_com, "yes") ~ "yes",
    str_detect(asynch_is_stab_at_com, "both") ~ "mixed",
    TRUE~asynch_is_stab_at_com))

# Figure 1C
# diversity bar

mono_div_bar<-trophic_lit %>%
  mutate(richness_stab = recode_factor(richness_stab, # standardize labels
                                Yes_and_No = "Mixed",
                                "Mixed_one season no" = "Mixed",
                                Both = "Mixed",
                                no_effect = "No effect",
                                Weak = "Yes")) %>% 
  filter(!pop_com_both == "population") %>%
  drop_na(richness_stab) %>%
  filter(trophic_study == "monotrophic") %>%
  group_by(richness_stab) %>%
  count() %>%
  ungroup() %>%
  mutate(Prop = n/sum(n)) %>%
  mutate(richness_stab = reorder(richness_stab,-1* Prop)) %>%
  ggplot(aes(x="", y=Prop, fill=richness_stab)) +
  
  geom_bar(aes(alpha =richness_stab), stat = "identity", fill = "#F8766D",color = "black", position = position_fill(reverse = T)) +
  geom_text(aes(label = richness_stab),position = position_stack(reverse = T, vjust = 0.5)) +
  theme_bw() +
  labs(fill = "Eff. of diversity") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        #legend.position = "bottom",
        legend.background = element_blank(), 
        legend.box.background = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +  
  removeGrid() 

ggsave("./Figures/Figure_1/mono_div_bar.pdf",
       plot=mono_div_bar,
       width = 1.5,
       height = 4)


multi_div_bar<-trophic_lit %>%
  mutate(richness_stab = recode_factor(richness_stab, # standardize labels
                                Yes_and_No = "Mixed",
                                "Mixed_one season no" = "Mixed",
                                Both = "Mixed",
                                no_effect = "No effect",
                                Weak = "Yes")) %>% 
  filter(!pop_com_both == "population") %>%
  drop_na(richness_stab) %>%
  filter(trophic_study == "multitrophic") %>%
  group_by(richness_stab) %>%
  count() %>%
  ungroup() %>%
  mutate(Prop = n/sum(n)) %>%
  mutate(richness_stab = reorder(richness_stab,-1* Prop)) %>%
  ggplot(aes(x="", y=Prop, fill=richness_stab)) +
  geom_bar(aes(alpha =richness_stab), stat = "identity", fill = "#6BAED6",color = "black", position = position_fill(reverse = F)) +
  geom_text(aes(label = richness_stab),position = position_stack(reverse = F, vjust = 0.5)) +
  theme_bw() +
  labs(fill = "Eff. of diversity") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.background = element_blank(), 
        legend.box.background = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +  
  removeGrid() 

ggsave("./Figures/Figure_1/multi_div_bar.pdf",
       plot=multi_div_bar,
       width = 1.5,
       height = 4)


# the order here needs to be corrected
multi_colors <-data.frame(levels = c("system stable","prey stablize","predators stablize","system unstable","predators destablize","no effect","other"),
                          colors = c("#56B4E9","blue","#0072b2","red","#FF6666","grey","brown"))


# Figure 1D - Vertical bar for asynchrony effects - Monotrophic
mono_bar_vert<-trophic_lit_upd %>%
  filter(trophic_study == "monotrophic",
         !asynch_is_stab_at_com == "NA") %>%    
  mutate(asynch_is_stab_at_com = recode_factor(asynch_is_stab_at_com, 
                                        yes = "destabilize", # this has been modified to fit the results for visual abstract
                                        "yes and no_pred-prey" = "mixed",
                                        both = "destabilize",  # this has been modified to fit the results for visual abstract
                                        no = "stablize",  # this has been modified to fit the results for visual abstract
                                        Mixed = "mixed")) %>%
  group_by(asynch_is_stab_at_com) %>%
  count() %>%
  ungroup() %>%
  mutate(Prop = n/sum(n)) %>%
  ggplot(aes(x="", y=Prop, fill=asynch_is_stab_at_com)) +
  #geom_bar(stat="identity",color="black",position = position_fill(reverse = TRUE))  +
  geom_bar(aes(alpha = asynch_is_stab_at_com),stat="identity", width=1, color="black",fill ="#F8766D",position = position_fill(reverse = TRUE)) +
  #coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c("red","white","#56B4E9"))+
  labs(fill = "Eff. of asynchrony") +
  geom_text(aes(x="",label = ifelse(Prop >= 0.2,str_wrap(asynch_is_stab_at_com, width = 10), "")),
            position = position_stack(vjust = 0.5, reverse = T), color = "black") +
  geom_text_repel(aes(label=ifelse(Prop <= 0.2, str_wrap(asynch_is_stab_at_com, width = 10), "")),
                  position = position_stacknudge(vjust = 0.5, x = -1, reverse = T),
                  size = 4,
                  direction = "x",
                  color = "black",
                  min.segment.length = unit(0, 'lines'),
                  force = 1) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        #legend.position = "bottom",
        legend.background = element_blank(), 
        legend.box.background = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +  
  removeGrid() 


ggsave("./Figures/Figure_1/mono_bar_vert.pdf",
       plot=mono_bar_vert,
       width = 1.5,
       height = 4)


# Figure 1D - Vertical bar for asynchrony effects - Multitrophic
multi_bar_vert_df<-trophic_lit %>%
  filter(trophic_study == "multitrophic",
         !pred_effect_red == "NA") %>%
  group_by(pred_effect_red) %>%
  mutate(pred_effect_red = gsub("_", " ", pred_effect_red)) %>%
  count() %>%
  ungroup() %>%
  mutate(Prop = n/sum(n)) %>%
  mutate(pred_effect_red = recode_factor(pred_effect_red, levels = c("system stable","prey stablize","predators stablize","system unstable","predators destablize","no effect","other")))

multi_bar_vert<-multi_bar_vert_df %>% 
  mutate(pred_effect_red = reorder(pred_effect_red,-1* Prop)) %>%
  ggplot(aes(x="", y=Prop, fill=pred_effect_red, label =pred_effect_red)) +
  geom_bar(stat="identity",color="black")  +
  theme_bw() +
  scale_fill_manual(values = pie_colors) +
  labs(fill = "Eff. of asynchrony") +
  guides(fill = guide_legend(nrow = 3)) +
  geom_text(aes(x="",label = ifelse(Prop >= 0.2,str_wrap(pred_effect_red, width = 10), "")),
            position = position_stack(vjust = 0.5), color = "black") +
  geom_text_repel(aes(label=ifelse(Prop <= 0.2, str_wrap(pred_effect_red, width = 10), "")),
                  position = position_stacknudge(vjust = 0.5, x = -1),
                  size = 4,
                  direction = "x",
                  color = "black",
                  min.segment.length = unit(0, 'lines'),
                  force = 1) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.background = element_blank(), 
        legend.box.background = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +
  removeGrid() 



ggsave("./Figures/Figure_1/multi_bar_vert.pdf",
       plot=multi_bar_vert,
       width = 1.5,
       height = 4)











## Figure S2 ###

# Figure S2A - temporal variation in monotrophic vs. multitrophic studies 
year_pub<-trophic_lit %>% 
  group_by(Year.of.publication,trophic_study) %>% 
  filter(!Year.of.publication %in% c("","na","MAYBE")) %>%
  dplyr::summarize(count = length(Year.of.publication)) %>% 
  ggplot(aes(Year.of.publication,count, fill = trophic_study)) +
  geom_bar(stat = 'identity', position = 'stack', size = 0.25) +
  theme_bw() +
  labs(x = "Year published", y = "no. papers") +
  scale_fill_manual(values = c("#F8766D","#6BAED6")) +
  theme(legend.position = c(0.3,0.8),
        legend.title = element_blank(),
        text = element_text(size=12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_x_continuous(breaks = seq(1992,2022,2)) +
  removeGrid()

# Figure S2B - type of organisms
organism_type<-trophic_lit %>% group_by(type_of_org,trophic_study) %>% 
  filter(!type_of_org %in% c("","na",NA,"MAYBE","MAYBE NOT INCLUDE","DON'T INCLUDE")) %>%
  mutate(type_of_org = tolower(type_of_org)) %>%
  dplyr::summarize(count = length(type_of_org)) %>% 
  ggplot(aes(reorder(type_of_org, -count),count, fill = trophic_study)) +
  geom_bar(stat = 'identity', position = 'stack', size = 0.25) +
  theme_bw() +
  xlab("Type of organism") +
  theme(legend.title = element_blank(),
        legend.position = "none",
        text = element_text(size=12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6)
        )+
  scale_fill_manual(values = c("#F8766D","#6BAED6")) +
  labs(x = "Type of organism", y = "no. papers") +
  scale_x_discrete(label=function(x) abbreviate(x, minlength=5)) +
  removeGrid()



# Figure 2C - effect of diversity across monotrophic vs. multitrophic studies
div_stab_plot<-trophic_lit %>%
  mutate(richness_stab = recode(richness_stab, # standardize labels
                                Yes_and_No = "Mixed",
                                "Mixed_one season no" = "Mixed",
                                Both = "Mixed",
                                no_effect = "No effect",
                                Weak = "Yes")) %>% 
  filter(!pop_com_both == "population") %>%
  group_by(trophic_study,richness_stab) %>%
  dplyr::summarize(count = length(richness_stab)) %>%
  ggplot(aes(x = reorder(richness_stab,-count),y = count, fill = trophic_study)) +
  geom_bar(stat = 'identity', position = 'stack', size = 0.25) +
  scale_fill_manual(values = c("#F8766D","#6BAED6")) +
  labs(x = "Stabilizing effect of diversity", y = "no. papers") +
  theme_bw() +
  removeGrid() +
  guides(fill = "none")

# Combine
exp_plot_V2 <-year_pub / organism_type / div_stab_plot + plot_layout() & theme(text = element_text(size = 12), legend.background = element_blank()) & plot_annotation(tag_levels = "A")

ggsave("./Figures/Figure_S2.pdf",
       plot=exp_plot,
       width = 4,
       height = 8)


ggsave("./Figures/MS/Figure_S2.pdf",
       plot=exp_plot_V2,
       width = 4,
       height = 8)


## Figure 2 ####

trophic_plot<-trophic_lit %>% group_by(trophic_study) %>% 
  summarize(count = length(trophic_study)) %>% 
  ggplot(aes(trophic_study,count, fill = trophic_study)) +
  geom_bar(stat = 'identity', position = 'dodge', size = 0.25) +
  theme_bw() +
  labs(x = "Study type", y = "no. papers") +
  theme(legend.title = element_blank(),
        legend.position = "none",
        #axis.title.y = element_blank(),
        text = element_text(size=12)) +
  scale_x_discrete(labels = c("monotrophic","multitrophic")) +
  geom_text(aes(x=trophic_study,label=count), vjust=-1) +
  ylim(0,350) +
  removeGrid()



trophic_lit$pred_effect_red <- factor(trophic_lit$pred_effect_red, 
                                       levels=c("predators_stablize", "predators_destablize", "prey_stablize", "prey_destabilize","system_stable","system_unstable","other","no_effect"))

multi_studies<-trophic_lit %>%
  filter(trophic_study == "multitrophic")

# number of predators 
trophic_stability<-trophic_lit %>%
  filter(trophic_study == "multitrophic") %>%
  mutate(study_type = if_else(study_type == "model","theoretical",study_type),
         pred_effect_red = replace_na(pred_effect_red,"no_effect")) %>% 
  group_by(pred_effect_red,study_type) %>% 
  summarize(count = length(pred_effect_red)) %>% 
  #ggplot(aes(reorder(pred_effect_red, -count),count)) +
  ggplot(aes(x = pred_effect_red,y = count, group = study_type)) +
  #ggplot(aes(x = pred_effect_red,y = count)) +
  #geom_bar(stat = 'identity', position = 'dodge', fill = "#00BFC4", size = 0.25) +
  geom_bar(aes(alpha = study_type), stat = 'identity', position = 'stack', size = 0.25,fill = "#00BFC4",colour="black") +
  theme_bw() +
  labs(x = "Multitrophic effect", y = "no. papers") +
  theme(legend.title = element_blank(),
        legend.position = c(0.2,0.88),
        legend.background = element_blank(),
        text = element_text(size=12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_text(aes(x=pred_effect_red,label=count),
            position = position_stack(vjust = .5)) +
  scale_x_discrete(labels = c("pred. stabilize", "pred. destab.", "prey stabilize", "system stable", "system unstable", "other", "no effect")) +
  ylim(0,30) +
  removeGrid()


trophic_stability_df<-trophic_lit %>%
  filter(trophic_study == "multitrophic",
         !pred_effect_red == "NA") %>%
  mutate(study_type = if_else(study_type == "model","theoretical",study_type))


year_pub_troph<-trophic_lit %>% 
  group_by(Year.of.publication,trophic_study) %>% 
  dplyr::summarize(count = length(Year.of.publication)) %>% 
  ggplot(aes(x = as.numeric(Year.of.publication), y = count, color = trophic_study)) +
  geom_point()+
  theme_bw() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) +
  stat_poly_eq(use_label(c("R2","P")),
               label.y = c(0.62,0.44)) +
  labs(x = "Year published", y = "no. papers", color = "Trophic study") +
  theme(legend.position = c(0.2,0.82),
        legend.background = element_blank(),
        text = element_text(size=12),
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous(breaks = seq(1992,2022,2)) +
  removeGrid()



Figure_2<-trophic_plot +  year_pub_troph + trophic_stability + plot_layout(ncol = 1, guides = "keep") & theme(axis.title = element_text(size = 12)) & plot_annotation(tag_levels = "A")


ggsave("./Figures/Figure_2.pdf",
       plot=Figure_2,
       width = 4,
       height = 10)


  
## Figure S3 ####
  
# Figure S3A - synchrony effects across biomes 
  
  mono_response_system<-trophic_lit_upd %>%
    filter(trophic_study == "monotrophic",
           !asynch_is_stab_at_com == "NA") %>%    
    mutate(asynch_is_stab_at_com = recode(asynch_is_stab_at_com, 
                                          yes = "stabilize",
                                          "yes and no_pred-prey" = "mixed",
                                          both = "stabilize",
                                          no = "destablize",
                                          Mixed = "mixed"),
           system_type = if_else(system_type == "AT", "MAT", system_type)) %>%
    group_by(asynch_is_stab_at_com,system_type) %>%
    count() %>%
    ggplot(aes(x=system_type, y=n, fill=asynch_is_stab_at_com, group = asynch_is_stab_at_com)) +
    geom_bar(stat="identity", width=1, position = 'stack',color="white") +
    geom_text(aes(x=system_type,label=ifelse(n > 2, n, "")),position = position_stack(vjust = 0.5), color = "black") +
    geom_text_repel(aes(label=ifelse(n <= 2, n, "")),
                    position = position_stacknudge(vjust = 0.5, x = 0.4),
                    size = 4,
                    color = "black",
                    min.segment.length = unit(0, 'lines')) +
    theme_bw() +
    removeGrid() +
    labs(y = "no. papers", fill = "Effect of asynchrony") +
    theme(axis.title.x = element_blank(),
          legend.background = element_blank(),
          legend.position=c(.3,.8)) +
    annotate("text", x = -Inf, y = Inf, label = "A",
             hjust = -.5, vjust = 1.5, size = 6, color = "black")
  
  # Figure S3B - synchrony effects across organizational scales
  
  mono_response_org<-trophic_lit_upd %>%
    filter(trophic_study == "monotrophic",
           !asynch_is_stab_at_com == "NA") %>%    
    mutate(asynch_is_stab_at_com = recode(asynch_is_stab_at_com, 
                                          yes = "stabilize",
                                          "yes and no_pred-prey" = "mixed",
                                          both = "stabilize",
                                          no = "destablize",
                                          Mixed = "mixed")) %>%
    group_by(asynch_is_stab_at_com,pop_com_both) %>%
    count() %>%
    ggplot(aes(x=pop_com_both, y=n, fill=asynch_is_stab_at_com, group = asynch_is_stab_at_com)) +
    geom_bar(stat="identity", width=1, position = 'stack',color="white") +
    geom_text(aes(x=pop_com_both,label=ifelse(n > 3, n, "")),position = position_stack(vjust = 0.5), color = "black") +
    geom_text_repel(aes(label=ifelse(n <= 3, n, "")),
                    position = position_stacknudge(vjust = 0.5, x = 0.4),
                    size = 4,
                    color = "black",
                    min.segment.length = unit(0, 'lines')) +
    theme_bw() +
    removeGrid() +
    labs(y = "no. papers", fill = "Effect of asynchrony") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none") +
    annotate("text", x = -Inf, y = Inf, label = "B",
             hjust = -.5, vjust = 1.5, size = 6, color = "black") +
    ylim(0,220)
  
  # Figure S3A & B combined
  
  mono_comparison<-mono_response_system + mono_response_org + plot_layout(ncol = 2)
  
  # Figure S3C - synchrony effects across biomes - multitrophic
  
  multi_response_system<-trophic_lit %>%
    filter(trophic_study == "multitrophic",
           !pred_effect_red == "NA") %>%
    group_by(pred_effect_red,system_type) %>%
    mutate(pred_effect_red = str_replace_all(pred_effect_red, "_", " "),
           system_type = if_else(system_type == "AT", "MAT", system_type)) %>%
    count() %>%
    ggplot(aes(x=system_type, y=n, fill=pred_effect_red, group = pred_effect_red)) +
    geom_bar(stat="identity", width=1, position = 'stack',color="white") +
    geom_text(aes(x=system_type,label=n),position = position_stack(vjust = 0.5), color = "black") +
    theme_bw() +
    removeGrid() +
    labs(y = "no. papers", fill = "Multitrophic effect") +
    theme(axis.title.x = element_blank(),
          legend.position="none") +
    annotate("text", x = -Inf, y = Inf, label = "C",
             hjust = -.5, vjust = 1.5, size = 6, color = "black")
  
  trophic_lit$pop_com_both <- factor(trophic_lit$pop_com_both, levels=c("community","population","predator-prey","plant-pollinator","host-parasite"))
  
  
  # Figure S3D - synchrony effects across organizational scales - multitrophic
  
  multi_response_org<-trophic_lit %>%
    filter(trophic_study == "multitrophic",
           !pred_effect_red == "NA") %>%
    group_by(pred_effect_red,pop_com_both) %>%
    mutate(pred_effect_red = str_replace_all(pred_effect_red, "_", " ")) %>%
    count() %>%
    ggplot(aes(x=pop_com_both, y=n, fill=pred_effect_red, group = pred_effect_red)) +
    geom_bar(stat="identity", width=1, position = 'stack',color="white") +
    geom_text(aes(x=pop_com_both,label=ifelse(n > 2, n, "")),position = position_stack(vjust = 0.5), color = "black") +
    geom_text_repel(aes(label=ifelse(n <= 2, n, "")),
                    position = position_stacknudge(vjust = 0.5, x = 0.4),
                    size = 4,
                    color = "black",
                    min.segment.length = unit(0, 'lines')) +
    theme_bw() +
    removeGrid() +
    labs(y = "no. papers", fill = "Multitrophic effect") +
    theme(axis.title.x = element_blank(),
          legend.position= c(0.8,0.6),
          legend.background = element_blank()) +
    annotate("text", x = -Inf, y = Inf, label = "D",
             hjust = -.5, vjust = 1.5, size = 6, color = "black") +
    ylim(c(0,50))
  
  
  
  
  
multi_comparison<-multi_response_system + multi_response_org + plot_layout(ncol = 2) #+ plot_annotation(tag_levels = "C")
  
mono_multi_comparison<-mono_comparison / multi_comparison + plot_layout(ncol = 1)

ggsave("./Figures/Figure_S3.pdf",
       plot=mono_multi_comparison,
       width = 10,
       height = 8)

ggsave("./Figures/Figure_S3.png",
       plot=mono_multi_comparison,
       width = 10,
       height = 8)




# Part B - Formal analysis for systematic review ####


 ## 1. Mono vs. multitrophic - study count ####
  T1_table<-trophic_lit %>%
    select(trophic_study) %>%
    mutate_all(~ if(is.character(.)) tolower(.) else .) %>%
    table()
  
  prop.table(T1_table)
  chisq.test(T1_table)
  
  T1_res<-data.frame(test = "T1", X = unclass(chisq.test(T1_table, simulate.p.value = TRUE))$statistic, p = unclass(chisq.test(T1_table, simulate.p.value = TRUE))$p.value)
  
 ## 2. Study design v. Realm (this is already done) - count - looking for differences in the evenness of distribution ####
  T2_table<-trophic_lit %>%
    select(trophic_study,system_type) %>%
    mutate_all(~ if(is.character(.)) tolower(.) else .) %>%
    table()
  
  chisq.test(T2_table)
  chisq.test(T2_table, simulate.p.value = TRUE)
  
  T2_res<-data.frame(test = "T2", X = unclass(chisq.test(T2_table, simulate.p.value = TRUE))$statistic, p = unclass(chisq.test(T2_table, simulate.p.value = TRUE))$p.value)
  
 ### 2a. Within monotrophic (not done) ####
  T2A_table<-trophic_lit %>%
    filter(trophic_study == "monotrophic") %>%
    select(system_type) %>%
    mutate_all(~ if(is.character(.)) tolower(.) else .) %>%
    table()
  
  chisq.test(T2A_table)
  chisq.test(T2A_table, simulate.p.value = TRUE)
  
  T2A_res<-data.frame(test = "T2A", X = unclass(chisq.test(T2A_table, simulate.p.value = TRUE))$statistic, p = unclass(chisq.test(T2A_table, simulate.p.value = TRUE))$p.value)
  
  ### 2b. Within multitrophic (not done) ####
  T2B_table<-trophic_lit %>%
    filter(trophic_study == "multitrophic") %>%
    select(system_type) %>%
    mutate_all(~ if(is.character(.)) tolower(.) else .) %>%
    table()
  
  chisq.test(T2B_table)
  chisq.test(T2B_table, simulate.p.value = TRUE)
  
  T2B_res<-data.frame(test = "T2B", X = unclass(chisq.test(T2B_table, simulate.p.value = TRUE))$statistic, p = unclass(chisq.test(T2B_table, simulate.p.value = TRUE))$p.value)
  
 ## 3. Study design v. Organizational scale (this is already done) - count - looking for differences in the evenness of distribution ####
  T3_table<-trophic_lit %>%
    select(trophic_study,pop_com_both) %>%
    mutate_all(~ if(is.character(.)) tolower(.) else .) %>%
    table()
  
  chisq.test(T3_table)
  chisq.test(T3_table, simulate.p.value = TRUE)
  
  T3_res<-data.frame(test = "T3", X = unclass(chisq.test(T3_table, simulate.p.value = TRUE))$statistic, p = unclass(chisq.test(T3_table, simulate.p.value = TRUE))$p.value)
  
 ### 3a. Within monotrophic (not done) ####
  T3A_table<-trophic_lit %>%
    filter(trophic_study == "monotrophic") %>%
    select(trophic_study,pop_com_both) %>%
    mutate_all(~ if(is.character(.)) tolower(.) else .) %>%
    table()
  chisq.test(T3A_table)
  chisq.test(T3A_table, simulate.p.value = TRUE)
  
  T3A_res<-data.frame(test = "T3A", X = unclass(chisq.test(T3A_table, simulate.p.value = TRUE))$statistic, p = unclass(chisq.test(T3A_table, simulate.p.value = TRUE))$p.value)
  
  
  T3A_extra_table<-trophic_lit %>%
    filter(trophic_study == "monotrophic") %>%
    select(pop_com_both,system) %>%
    mutate_all(~ if(is.character(.)) tolower(.) else .) %>%
    table()
  
  chisq.test(T3A_extra_table)
  chisq.test(T3A_extra_table, simulate.p.value = TRUE)
  
  
 ### 3b. Within multitrophic (not done) ####
  T3B_table<-trophic_lit %>%
    filter(trophic_study == "multitrophic") %>%
    select(trophic_study,pop_com_both) %>%
    mutate_all(~ if(is.character(.)) tolower(.) else .) %>%
    table()
  chisq.test(T3B_table)
  chisq.test(T3B_table, simulate.p.value = TRUE)
  
  T3B_res<-data.frame(test = "T3B", X = unclass(chisq.test(T3B_table, simulate.p.value = TRUE))$statistic, p = unclass(chisq.test(T3B_table, simulate.p.value = TRUE))$p.value)
  
 ## 4. Effect of synchrony on stability in monotrophic ####
  unique(trophic_lit$asynch_is_stab_at_com)
  T4_table<-trophic_lit %>%
    filter(trophic_study == "monotrophic") %>%
    filter(!asynch_is_stab_at_com %in% c("","na",NA,"MAYBE NOT INCLUDE","DON'T INCLUDE")) %>%
    mutate_all(~ if(is.character(.)) tolower(.) else .) %>%
    mutate(asynch_is_stab_at_com = recode(asynch_is_stab_at_com, 
                                          yes = "stabilize",
                                          "yes and no_pred-prey" = "mixed",
                                          both = "stabilize",
                                          no = "destablize",
                                          Mixed = "mixed")) %>%
    select(asynch_is_stab_at_com) %>%
    table()
  
  prop.table(T4_table)
  chisq.test(T4_table)
  chisq.test(T4_table, simulate.p.value = TRUE)
  
  T4_res<-data.frame(test = "T4", X = unclass(chisq.test(T4_table, simulate.p.value = TRUE))$statistic, p = unclass(chisq.test(T4_table, simulate.p.value = TRUE))$p.value)
  
 ### 4A. Within realm ####
  T4A_table<-trophic_lit %>%
    filter(trophic_study == "monotrophic") %>%
    filter(!asynch_is_stab_at_com %in% c("","na",NA,"MAYBE NOT INCLUDE","DON'T INCLUDE")) %>%
    mutate_all(~ if(is.character(.)) tolower(.) else .) %>%
    mutate(asynch_is_stab_at_com = recode(asynch_is_stab_at_com, 
                                          yes = "stabilize",
                                          "yes and no_pred-prey" = "mixed",
                                          both = "stabilize",
                                          no = "destablize",
                                          Mixed = "mixed")) %>%
    select(asynch_is_stab_at_com,system_type) %>%
    table()
  
  prop.table(T4A_table)
  chisq.test(T4A_table)
  chisq.test(T4A_table, simulate.p.value = TRUE)
  
  T4A_res<-data.frame(test = "T4A", X = unclass(chisq.test(T4A_table, simulate.p.value = TRUE))$statistic, p = unclass(chisq.test(T4A_table, simulate.p.value = TRUE))$p.value)
  
 ### 4B. Across organizational scale ####
  T4B_table<-trophic_lit %>%
    filter(trophic_study == "monotrophic") %>%
    filter(!asynch_is_stab_at_com %in% c("","na",NA,"MAYBE NOT INCLUDE","DON'T INCLUDE")) %>%
    mutate_all(~ if(is.character(.)) tolower(.) else .) %>%
    mutate(asynch_is_stab_at_com = recode(asynch_is_stab_at_com, 
                                          yes = "stabilize",
                                          "yes and no_pred-prey" = "mixed",
                                          both = "stabilize",
                                          no = "destablize",
                                          Mixed = "mixed")) %>%
    select(asynch_is_stab_at_com,pop_com_both) %>%
    filter(pop_com_both %in% c("community","population")) %>% 
    droplevels %>%
    table()
    

  prop.table(T4B_table)
  chisq.test(T4B_table)
  chisq.test(T4B_table, simulate.p.value = TRUE)
  
  T4B_res<-data.frame(test = "T4B", X = unclass(chisq.test(T4B_table, simulate.p.value = TRUE))$statistic, p = unclass(chisq.test(T4B_table, simulate.p.value = TRUE))$p.value)
  
 ## 5. Effect of synchrony on stability in multitrophic ####
  T5_table<- trophic_lit %>%
    filter(trophic_study == "multitrophic",
           !pred_effect_red == "NA") %>%
    mutate_all(~ if(is.character(.)) tolower(.) else .) %>%
    select(pred_effect_red) %>%
    table()
  
  prop.table(T5_table)
  chisq.test(T5_table)
  chisq.test(T5_table, simulate.p.value = TRUE)
  
  T5_res<-data.frame(test = "T5", X = unclass(chisq.test(T5_table, simulate.p.value = TRUE))$statistic, p = unclass(chisq.test(T5_table, simulate.p.value = TRUE))$p.value)
  
 ### 5A. Within realm  ####
  T5A_table<- trophic_lit %>%
    filter(trophic_study == "multitrophic",
           !pred_effect_red == "NA") %>%
    mutate_all(~ if(is.character(.)) tolower(.) else .) %>%
    mutate(pred_effect_red = factor(pred_effect_red)) %>%
    select(pred_effect_red,system_type) %>% 
    table()
  
  
  prop.table(T5A_table)
  chisq.test(T5A_table)
  chisq.test(T5A_table, simulate.p.value = TRUE)
 
  T5A_res<-data.frame(test = "T5A", X = unclass(chisq.test(T5A_table, simulate.p.value = TRUE))$statistic, p = unclass(chisq.test(T5A_table, simulate.p.value = TRUE))$p.value)
  
  ### 5B. Across organizational scale ####
  T5B_table<- trophic_lit %>%
    filter(trophic_study == "multitrophic",
           !pred_effect_red == "NA") %>%
    mutate_all(~ if(is.character(.)) tolower(.) else .) %>%
    mutate(pred_effect_red = factor(pred_effect_red)) %>%
    select(pred_effect_red,pop_com_both) %>%
    table()
  
  prop.table(T5B_table)
  chisq.test(T5B_table)
  T5B_res<-data.frame(test = "T5B", X = unclass(chisq.test(T5B_table, simulate.p.value = TRUE))$statistic, p = unclass(chisq.test(T5B_table, simulate.p.value = TRUE))$p.value)
  
  
  
  ## 6. Effect of diversity on stability ####
  div_stab_table<-trophic_lit %>%
    mutate(richness_stab = recode(richness_stab, # standardize labels
                                  Yes_and_No = "Mixed",
                                  "Mixed_one season no" = "Mixed",
                                  Both = "Mixed",
                                  no_effect = "No effect",
                                  Weak = "Yes"),
           richness_stab = if_else(is.na(richness_stab),"not_assessed",richness_stab))%>% 
    filter(!pop_com_both == "population") %>%
    #group_by(trophic_study,richness_stab) %>%
    #dplyr::summarize(count = length(richness_stab)) %>%
    select(richness_stab,richness_stab) %>%
    table()
  
  
  prop.table(div_stab_table)
  chisq.test(div_stab_table)
  T6_res<-data.frame(test = "T6", X = unclass(chisq.test(div_stab_table, simulate.p.value = TRUE))$statistic, p = unclass(chisq.test(div_stab_table, simulate.p.value = TRUE))$p.value)
  
  
  
  # Add these all to a table
  chi_results<-rbind(T1_res,T2_res,T2A_res,T2B_res,T3_res,T3A_res,T3B_res,T4_res,T4A_res,T4B_res,T5_res,T5A_res,T5B_res,T6_res)
  
  

  organ_test_df<-trophic_lit %>%
    select(trophic_study,system_type)
  
  organ_test_3_df<-trophic_lit %>%
    select(trophic_study,pop_com_both,system_type) %>%
    mutate_all(as.factor)

  organ_test_3_table <- table(organ_test_3_df$trophic_study,
                              organ_test_3_df$pop_com_both,
                              organ_test_3_df$system_type)
  
  organ_test_3_table <- table(organ_test_3_df)
  
  is.table(organ_test_3_table)

  
  mantelhaen.test(organ_test_3_table)
  chisq.test(organ_test_3_table)
  
  # Two tests
  # (1) Differences in frequency of system across monotrophic vs. multitrophic study
  T1_chi_table<-trophic_lit %>%
    select(trophic_study,system_type) %>%
    mutate_all(~ if(is.character(.)) tolower(.) else .) %>%
    table()
  chisq.test(T1_chi_table)
  chisq.test(T1_chi_table, simulate.p.value = TRUE)
  
  # (2) Differences in frequency of organizational scale across monotrophic vs. multitrophic study
  T2_chi_table<-trophic_lit %>%
    select(trophic_study,pop_com_both) %>%
    mutate_all(~ if(is.character(.)) tolower(.) else .) %>%
    table()
  chisq.test(T2_chi_table)
  chisq.test(T2_chi_table, simulate.p.value = TRUE)
  
  
### END ###
