########## CLUSTER ANALYSIS OF CCFPS STRIPED BASS MOVEMENTS ###############
# Setup ###################################################################
source(file.path("2. Code","0_Setup.R"))
dat <- readRDS(file.path("1. Data", "Outputs" , "Cluster_Analysis_Data.rds"))

library(vegan)
library(factoextra)
library(cluster)
library(flexclust)
set.seed(123)

## Formatting ==============================================================
datT1 <- decostand(dat, method = "total")
datT1 <- decostand(datT1, method = "max")
datT2 <- datT1 %>%
  select(-AvgRes_INSIDE,-AvgRes_OUTSIDE,
         q25_distance, avg_time_btwn_detects, lifetime_outside, - TAL)
datBray <- vegdist(datT2, na.rm = TRUE)


# Exploration ==============================================================

## Correlation plot ========================================================
corr = cor(datT2)
testRes = corrplot::cor.mtest(datT2, conf.level = 0.95)
corrplot::corrplot(corr, p.mat = testRes$p, method = 'circle', type = 'lower', 
                   insig='blank',
                   #addCoef.col ='black', 
                   number.cex = 0.2, order = 'AOE', diag=FALSE)

## Principal Components ====================================================
fit <- princomp(datT2, cor = TRUE)
summary(fit) # print variance accounted for
loadings(fit) # pc loadings
plot(fit,type="lines") # scree plot
fit$scores # the principal components scores for each individual
biplot(fit, cex = c(0.1,.75), scale = 1) #PCA Biplot

pdf(file.path("1. Data","Figures","PCA_biplot.pdf"), width = 10, height = 10)
biplot(fit, cex = c(0.1,.75), scale = 1) #PCA Biplot
dev.off()

## Dendogram ===============================================================
datBray <- vegdist(datT2)
bray_agnes <- cluster::agnes(datBray, method = "ward")

pdf(file.path("1. Data","Figures","Dendrogram.pdf"), 
    width = 8, height = 6)
plot(bray_agnes, which.plots = 2, main="Movement Groups", cex=0.1)
abline(h=3.5, col="blue", lwd = 2)
abline(h=4.5, col="green", lwd = 2)
dev.off()

## Identify Potential Clusters =============================================

### WSS Plot ---------------------------------------------------------------
# Want to look for the bend in the curve
pdf(file.path("1. Data","Figures","WSS_Plot.pdf"), 
    width = 8, height = 6)
fviz_nbclust(datT2, kmeans, method = "wss", nboot = 1000)
dev.off ()

### Silhouette Plot --------------------------------------------------------
# Want to get the highest mean sihouette 
pdf(file.path("1. Data","Figures","Sil_Plot.pdf"), 
    width = 8, height = 6)
fviz_nbclust(datT2, kmeans, method = "silhouette", nboot = 1000)
dev.off()

### Gap Statistic ----------------------------------------------------------
gap_stat <- cluster::clusGap(datT2, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)
pdf(file.path("1. Data","Figures","GapStat_Plot.pdf"), 
    width = 8, height = 6)
fviz_gap_stat(gap_stat)
dev.off()

# Perform Clustering #######################################################
# 3 to 5 clusters has relatively good support:
# Dendrogram: 3 or 4
# WSS Plot: 3 to 5
# Silhouette Plot: 5
# Gap Statistic: 10 (lots of overlaps)

# K-means Cluster with k=3 
kmeans_3 <- eclust(datT2, "kmeans", k = 3, nstart = 25, nboot = 1000, 
                   verbose = TRUE)

pdf(file.path("1. Data","Figures","ClusterOutput.pdf"), 
    width = 8, height = 6)
fviz_cluster(kmeans_3, axes = c(1,3))+
  theme_classic() 
dev.off()

# Assign the resulting cluster ID to each individual
cluster_assign <- data.frame("cluster" = factor(kmeans_3$cluster,
                                                levels = c(1,2,3))) %>%
  rownames_to_column("TagAge")

assigned_dat <- dat %>%
  bind_cols(cluster_assign)

# Interpretation ##########################################################
## ANOVA ==================================================================
# Run anova on all the variables, all vars should be significant
anova <- aov(cbind(InsideSiteVisits, OutsideSiteVisits, numSites, TAL,                   
                   lifetime_inside, lifetime_outside, lifetime_all, AgeBin, 
                   mean_WkRcvrs, max_WkRcvrs, CaptureAgeBin, 
                   max_time_btwn_detects, avg_time_btwn_detects,
                   mean_time_btwn_mvmts, max_time_btwn_mvmts, mean_distance,       
                   q25_distance, q50_distance, q75_distance, q99_distance, 
                   emmigrant, FirstExitAgeBin, mean_exit_woy, mean_entry_woy,     
                   Entry, Exit, All_Moves, Exit_freq, Entry_freq, All_Moves_freq,       
                   AvgRes_INSIDE, AvgRes_OUTSIDE, TotRes_INSIDE,TotRes_OUTSIDE
                   ) ~ cluster, 
             data = assigned_dat)
# get a summary of anova results
summary(anova) # All except q25 are significant

# View the means of each variable for each cluster
assigned_dat %>%
  group_by(cluster) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE))) %>%
  select(-q25_distance) %>% #remove insignificant value
  t()

### Feature Importance ====================================================
# Identify those features which have the most weight in the clustering process

#devtools::install_github("o1iv3r/FeatureImpCluster")

# FeatureImpCluster requires input from flexclust so convert
res <- flexclust::as.kcca(kmeans_3, data = datT2, k = 3)

# Get the Feature Importance
set.seed(916)
nr_seeds <- 100
seeds_vec <- sample(1:1000, nr_seeds)

savedImp <- data.frame(matrix(0,nr_seeds,dim(datT2)[2]))
count <- 1
for (s in seeds_vec) {
  set.seed(s)
  res <- flexclust::kcca(datT2, k=3)
  set.seed(s)
  FeatureImp_res <- FeatureImpCluster::FeatureImpCluster(res,as.data.table(datT2),
                                                         sub = 1,biter = 1)
  savedImp[count,] <- FeatureImp_res$featureImp[sort(names(FeatureImp_res$featureImp))]
  count <- count + 1
}
names(savedImp) <- sort(names(FeatureImp_res$featureImp))

savedImp <- savedImp  %>%
  pivot_longer(cols = AgeBin:TotRes_OUTSIDE, 
               names_to = "variable",
               values_to = "importance")

savedImp %>%  
  group_by(variable) %>%
  summarise(mean_imp = mean(importance)) %>%
  arrange(desc(mean_imp)) %>%
  pull(variable) -> var_orders

ggplot(savedImp %>%
         mutate(importance = ifelse(importance == 0, 0.0001,
                                    importance)) %>%
         group_by(variable) #%>%
         #mutate(mean_importance = mean(importance))
       , 
       aes(x = factor(variable,levels = var_orders), y = importance))+
  geom_point()+
  geom_boxplot()+
  geom_smooth()+
  theme_classic()+
  scale_y_log10()+
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0),
    axis.text.y = element_blank()
  )+
  labs(x = "", y = "Importance")
ggsave(file.path("1. Data","Figures","FeatureImportance.pdf"), device = "pdf",
       width = 8, height = 6)

# q50_distance, max_WkRcvrs, avg_time_btwn_detects, lifetime_outside, mean_distance,
# q25_distance, q75_distance, are the lowest scoring.



## Examine the variables by cluster --------------------------------------

# Residency: Total Residency is the #1 Strongest predictor of Cluster
residency <- ggplot(assigned_dat %>%
                      select(cluster, TotRes_INSIDE, TotRes_OUTSIDE) %>%
                      pivot_longer(cols = TotRes_INSIDE:TotRes_OUTSIDE,
                                   names_to = c("stat","location"),
                                   values_to = "value",
                                   names_sep = "_"))+
  geom_boxplot(aes(x = location, y = value, fill = location))+
  facet_grid(cluster ~ .)+
  scale_fill_discrete(name = "Residency Location",
                      breaks = c("INSIDE","OUTSIDE"),
                      labels = c("Inside","Outside"))+
  scale_y_continuous(name = "Total Weeks of Residency")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        legend.position = "top",
        axis.text.y = element_blank(),
        axis.title.y = element_blank())+
  coord_flip()
residency
ggsave(file.path("1. Data","Figures","Cluster_01_Residency.pdf"), device = "pdf",
       width = 8, height = 6)

# Emmigrant Status was Second
emmigrant_plot <- ggplot(assigned_dat %>%
                           group_by(cluster) %>%
                           mutate(total = n()) %>%
                           group_by(emmigrant, cluster) %>%
                           summarise(prop = n()) %>%
                           distinct())+
  geom_bar(aes(x = emmigrant, y = prop, fill = emmigrant), stat = "identity")+
  scale_fill_discrete(name = "",
                      breaks = c(TRUE,FALSE),
                      labels = c("Emigrant","Non-Emigrant"))+
  facet_grid(cluster ~ .)+
  coord_flip()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.text.y = element_blank(),
        legend.position = "top")+
  labs(x = "", y = "Number of Individuals")
emmigrant_plot
ggsave(file.path("1. Data","Figures","Cluster_02_Emmigration.pdf"), device = "pdf",
       width = 8, height = 6)

# 3rd was Total Number of Inside Sites Visited 
sites_plot1 <- ggplot(assigned_dat %>%
                        select(cluster,TagAge,"Inside" = InsideSiteVisits,
                               "Outside" = OutsideSiteVisits) %>%
                        pivot_longer(cols = Inside:Outside,
                                     names_to = "Site",
                                     values_to = "Total"))+
  geom_boxplot(aes(y = Total, x = Site, fill = Site))+
  facet_grid(cluster ~ ., scales = 'free_y')+
  coord_flip()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0))+
  labs(y = "Total Number of Site Visits at Each Location")
sites_plot1
ggsave(file.path("1. Data","Figures","Cluster_03_TotalSites.pdf"), device = "pdf",
       width = 8, height = 6)

# 4th and 5th are the Average Weekly Site Visits
sites_plot2 <- ggplot(assigned_dat %>%
                        select(cluster,TagAge,"Inside" = lifetime_inside,
                               "Outside" = lifetime_outside,
                               "All" = lifetime_all) %>%
                        pivot_longer(cols = Inside:All,
                                     names_to = "Site",
                                     values_to = "Total"))+
  geom_boxplot(aes(y = Total, x = Site, fill = Site))+
  facet_grid(cluster ~ ., scales = 'free_y')+
  coord_flip()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0))+
  labs(y = "Average Number of Sites Visits per Week at Large")
sites_plot2
ggsave(file.path("1. Data","Figures","Cluster_04_AvgWeeklySites.pdf"), device = "pdf",
       width = 8, height = 6)

# 6th was 99-%tile Weekly Distance
distance_plot <- ggplot(assigned_dat %>% 
                          select(cluster, q25_distance:q99_distance) %>%
                          pivot_longer(cols = q25_distance:q99_distance, 
                                       names_to = "metric",
                                       values_to = "distance"),
                        aes(x = metric, y = distance, group = cluster))+
  stat_summary(fun.data=mean_se,geom="linerange",color="black", size = 1)+
  stat_summary(fun.y=mean,geom="line", size = 01)+
  geom_hline(aes(yintercept = 1.6), linetype = "dashed", linewidth = 0.5,
             color = "blue",alpha = 0.65)+
  geom_text(aes(x = 0.5, y = 1.9, label = "1/2 Distance From Intake Canal to Radial Gates"),
            hjust = 0)+
  geom_hline(aes(yintercept = 3.8), linetype = "dashed", linewidth = 0.5,
             color = "blue",alpha = 0.65)+
  geom_text(aes(x = 0.5, y = 4.1, label = "Distance From Intake Canal to ORS1"),
            hjust = 0)+
  geom_hline(aes(yintercept = 5.1), linetype = "dashed", linewidth = 0.5,
             color = "blue",alpha = 0.65)+
  geom_text(aes(x = 0.5, y = 5.4, label = "Distance From Intake Canal to Closest Array Boundary"),
            hjust = 0)+
  facet_grid(cluster ~ .)+
  scale_x_discrete(name = "Quantile",
                   breaks = c("q25_distance","q50_distance","q75_distance",
                              "q99_distance"),
                   labels = c("25%","50%","70%","99%"))+
  scale_y_continuous(name = "Mean Weekly Distance Travelled")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        panel.grid.major.y = element_line(color = "grey90"),
        panel.background = element_rect(color = "grey50"))
distance_plot
ggsave(file.path("1. Data","Figures","Cluster_05_WeeklyDistance.pdf"), device = "pdf",
       width = 8, height = 6)

# 7th was Time between Movements
time_move_plot <- ggplot(assigned_dat %>%
                           select(cluster,TagAge,mean_time_btwn_mvmts,
                                  max_time_btwn_mvmts) %>%
                           pivot_longer(cols = mean_time_btwn_mvmts:max_time_btwn_mvmts,
                                        names_to = "Statistic",
                                        values_to = "value"))+
  geom_boxplot(aes(y = value, x = 1, fill = Statistic))+
  facet_grid(cluster ~ ., scales = 'free_y')+
  coord_flip()+
  theme_classic()+
  theme()+
  scale_x_discrete(breaks = c("mean_time_btwn_mvmts"),
                   labels = c("Mean"))+
  labs(x = "", y = "Weeks between Transits")
time_move_plot
ggsave(file.path("1. Data","Figures","Cluster_06_TimeBtwnMovements.pdf"), device = "pdf",
       width = 8, height = 6)

wk_rcvrs_plot <- ggplot(assigned_dat %>%
                          select(cluster,TagAge,mean_WkRcvrs,max_WkRcvrs) %>%
                          pivot_longer(cols = mean_WkRcvrs:max_WkRcvrs,
                                       names_to = "stat",
                                       values_to = "value"))+
  geom_violin(aes(y = value, x = stat, fill = stat))+
  facet_grid(cluster ~ ., scales = 'free_y')+
  coord_flip()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0))+
  scale_x_discrete(breaks = c("mean_WkRcvrs","max_WkRcvrs"),
                   labels = c("Mean", "Max"))+
  labs(x = "Statistic", y = "Number of Receivers per Week")+
  theme(legend.position = "none")
wk_rcvrs_plot
ggsave(file.path("1. Data","Figures","Cluster_07_WeeklyReceiverVisits.pdf"), device = "pdf",
       width = 8, height = 6)

# 8th is the Age of Various Events
age_plot <- ggplot(assigned_dat %>%
                     select(cluster, FirstExitAgeBin, CaptureAgeBin,
                            AgeBin) %>%
                     pivot_longer(cols = FirstExitAgeBin:AgeBin,
                                  names_to = "event",values_to = "age_bin") %>%
                     mutate(event = case_when(
                       event == "FirstExitAgeBin" ~ "First Exit",
                       event == "CaptureAgeBin" ~ "Capture",
                       event == "AgeBin" ~ "Current"
                     )))+
  geom_violin(aes(x = event, fill = event, y = age_bin))+
  facet_grid(cluster ~ .)+
  coord_flip()+
  scale_fill_discrete(name = "Event")+
  scale_y_continuous(name = "Age of Event",
                     breaks = c(1,2,3,4),
                     labels = c("0","1-2","3-5","6+"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())+
  theme(legend.position = "top")
age_plot
ggsave(file.path("1. Data","Figures","Cluster_08_Age.pdf"), device = "pdf",
       width = 8, height = 6)

# 9th is Movements (Specifically Exits)
movements <- ggplot(assigned_dat %>%
                      select(cluster,Exit,Entry, 
                             All_Moves) %>%
                      pivot_longer(cols = Exit:All_Moves,
                                   names_to = "movement",
                                   values_to = "count"))+
  geom_boxplot(aes(y = count,
                   x = 1,
                   fill = movement),
  )+
  facet_grid(cluster ~ .)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top"
  )+
  scale_y_continuous(name = "Number of Transits",
                     breaks = seq(0,10,2))+
  scale_x_continuous(name = "", limits= c(0.3,1.6))+
  scale_fill_discrete(name = "")+
  coord_flip()
movements
ggsave(file.path("1. Data","Figures","Cluster_08_Age.pdf"), device = "pdf",
       width = 8, height = 6)

mean_exit_woy <- ggplot(assigned_dat%>%
                          separate(TagID,into = c("TagID","AgeBin"), sep = "-",
                                   extra = "merge") %>%
                          filter(TagID %in% c(emmigrant_tags$TagID)))+
  geom_violin(aes(y = lubridate::ymd("2017-12-31")+(mean_exit_woy*7), x = 1), 
              fill = "grey50",
              trim = FALSE, bw = 7)+
  facet_grid(cluster ~ .,drop = FALSE)+
  coord_flip()+
  theme_classic()+
  scale_y_date(name = "Date", date_labels = "%b-%d")+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  labs(title = "Mean Week of Exit")
mean_exit_woy

mean_entry_woy <- ggplot(assigned_dat%>%
                           separate(TagID,into = c("TagID","AgeBin"), sep = "-",
                                    extra = "merge") %>%
                           filter(mean_entry_woy > 0))+
  geom_violin(aes(y = lubridate::ymd("2017-12-31")+(mean_entry_woy*7), x = 1), fill = "grey50",
              scale = "width",
              trim = FALSE, bw = 7)+
  facet_grid(cluster ~ .,drop = FALSE)+
  coord_flip()+
  theme_classic()+
  scale_y_date(name = "Date", date_labels = "%b-%d")+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  labs(title = "Mean Week of Entry")
mean_entry_woy 

for(i in unique(assigned_dat$TagID)){
  tmp <- TagBKM_Bin %>%
    filter(TagID == i) %>%
    filter(SiteCode != "Tag Failure") %>%
    mutate(TagID = as.character(TagID)) %>%
    left_join(cluster_assign)

  if(nrow(tmp)>1){
    ggplot()+
      geom_point(data = tmp,
                 aes(x = studyweek_startdate(Week), y = SiteCode, group = TagID))+
      geom_step(data = tmp,
                aes(x = studyweek_startdate(Week), y = SiteCode, group = TagID))+
      scale_x_date(name = "Date",
                   breaks = scales::breaks_pretty(n = 8),
                   date_labels = "%b-%d")+
      scale_y_discrete()+
      theme_classic()+
      theme(panel.spacing = unit(0, "lines"),
            axis.text.x = element_text(angle = -45, hjust = 0))+
      facet_grid(TagID ~.)
    ggsave(file.path("1. Data","Figures","DetectionPlots",paste(tmp$TagID[1],tmp$cluster[1],"ex.pdf",sep = "-")),
           device = "pdf", width = 5, height = 2.5)
  } else {
    print("Tag not in cluster")
  }
}
