########## CLUSTER ANALYSIS OF CCFPS STRIPED BASS MOVEMENTS ###############
# Setup ###################################################################
source(file.path("2. Code","0_Setup.R"))
dat <- readRDS(file.path("1. Data", "Outputs" , "Cluster_Analysis_Data.rds"))

library(vegan)
library(factoextra)
library(cluster)
library(flexclust)
library(dendextend)
set.seed(123)

## Formatting ==============================================================
datT1_mig <- dat %>% 
  filter(emmigrant ==TRUE) %>% 
  decostand(method = "standardize") %>% 
  select(-emmigrant)
datT1_res <- dat %>% 
  filter(emmigrant ==FALSE) %>% 
  decostand(method = "standardize") %>%
  select(-emmigrant, -outside_sites_total, -outside_site_freq, -exit_total, 
         -exit_freq, -entry_total, -entry_freq, -all_total, -all_freq,
         -residence_max_OUTSIDE, -residence_mean_OUTSIDE, -residence_total_OUTSIDE,
         -mean_exit_woy, -mean_entry_woy)
# Exploration ==============================================================

## Correlation plot ========================================================
corr_mig = cor(datT1_mig)
corr_res = cor(datT1_res)

testRes_mig = corrplot::cor.mtest(datT1_mig, conf.level = 0.95)
testRes_res = corrplot::cor.mtest(datT1_res, conf.level = 0.95)

corrplot::corrplot(corr_mig, p.mat = testRes_mig$p, method = 'circle', type = 'lower', 
                   insig='blank',
                   #addCoef.col ='black', 
                   number.cex = 0.2, order = 'AOE', diag=FALSE)

corrplot::corrplot(corr_res, p.mat = testRes_res$p, method = 'circle', type = 'lower', 
                   insig='blank',
                   #addCoef.col ='black', 
                   number.cex = 0.2, order = 'AOE', diag=FALSE)

corr_mig <- as_tibble(corr_mig)
corr_res <- as_tibble(corr_res)
(corr_mig <- corr_mig %>%
  mutate(V1 = names(corr_mig)) %>%
  pivot_longer(cols = AgeBin:gate_perc_open,
               names_to = "V2",values_to = "corr") %>%
  filter(V1 != V2) %>%
  arrange(desc(corr)))
(corr_res <- corr_res %>%
  mutate(V1 = names(corr_res)) %>%
  pivot_longer(cols = AgeBin:gate_perc_open,
               names_to = "V2",values_to = "corr") %>%
  filter(V1 != V2) %>%
  arrange(desc(corr)))

# Remove Variables with >95% correlation (same variable only)
datT2_mig <- datT1_mig %>%
  select(-time_btwn_mvmts_max, -residence_max_OUTSIDE,
         -CaptureAgeBin)
datT2_res <- datT1_res %>%
  select(-time_btwn_detects_max, -time_btwn_mvmts_max,
         -CaptureAgeBin, all_site_freq, all_sites_total, )

## Principal Components ====================================================
fit_mig <- prcomp(datT2_mig)
fit_res <- prcomp(datT2_res)

summary(fit_mig) # print variance accounted for
summary(fit_res) # print variance accounted for


plot(fit_mig, type="lines") # scree plot
plot(fit_res,type="lines") # scree plot

biplot(fit_mig, cex = c(0.1,.75), scale = 1) #PCA Biplot
biplot(fit_res, cex = c(0.1,.75), scale = 1) #PCA Biplot

png(file.path("1. Data","Figures","Migrant_PCA_biplot.png"), width = 10, height = 10,
    units = "in", res = 600)
biplot(fit_mig, cex = c(0.1,.75), scale = 1) #PCA Biplot
dev.off()

png(file.path("1. Data","Figures","Resident_PCA_biplot.png"), width = 10, height = 10,
    units = "in", res = 600)
biplot(fit_res, cex = c(0.1,.75), scale = 1) #PCA Biplot
dev.off()

## Dendogram ===============================================================
datEuc_mig <- vegdist(datT2_mig, method = "euclidean", na.rm = TRUE)
datEuc_res <- vegdist(datT2_res, method = "euclidean", na.rm = TRUE)

euc_agnes_mig <- cluster::agnes(datEuc_mig, method = "ward")
euc_agnes_res <- cluster::agnes(datEuc_res, method = "ward")


labelCol_season <- function(x) {
  if (is.leaf(x)) {
    ## fetch label
    label <- attr(x, "label") 
    ## set label color to red for A and B, to blue otherwise
    attr(x, "nodePar") <- list(lab.col=case_when(grepl("Spawn",label) ~ "orange",
                                                 grepl("Immigration",label) ~ "forestgreen",
                                                 grepl("Residence",label) ~ "plum",
                                                 grepl("Emmigration", label) ~ "steelblue"),
                               lab.cex = 0.2,
                               cex = NA, pch = NA)
  }
  return(x)
}

labelCol_age <- function(x) {
  if (is.leaf(x)) {
    ## fetch label
    label <- attr(x, "label") 
    ## set label color to red for A and B, to blue otherwise
    attr(x, "nodePar") <- list(lab.col=case_when(grepl("1-2",label) ~ "orange",
                                                 grepl("3-5",label) ~ "forestgreen",
                                                 grepl("6+", label) ~ "steelblue"),
                               lab.cex = 0.2,
                               cex = NA, pch = NA)
  }
  return(x)
}

dg_mig <- as.dendrogram(euc_agnes_mig)
dg_mig <- dendextend::remove_nodes_nodePar(dg_mig)
dg_mig_season <- dendrapply(dg_mig, labelCol_season)
dg_mig_age <- dendrapply(dg_mig, labelCol_age)

dg_res <- as.dendrogram(euc_agnes_res)
dg_res <- dendextend::remove_nodes_nodePar(dg_res)
dg_res_season <- dendrapply(dg_res, labelCol_season)
dg_res_age <- dendrapply(dg_res, labelCol_age)

plot(dg_mig_age)
plot(dg_res_age)

png(file.path("1. Data","Figures","Migrant_Dendrogram.png"), 
    width = 12, height = 6,
    units = "in", res = 600)
par(mar = c(0.1,3,0.1,0.11))
plot(dg_mig_age)
#dendextend::rect.dendrogram(dg, k = 5)
#abline(h=5, col="blue", lwd = 2)
#abline(h=4, col="green", lwd = 2)
dev.off()

png(file.path("1. Data","Figures","Resident_Dendrogram.png"), 
    width = 12, height = 6,
    units = "in", res = 600)
par(mar = c(0.1,3,0.1,0.11))
plot(dg_res_age)
#dendextend::rect.dendrogram(dg, k = 5)
#abline(h=5, col="blue", lwd = 2)
#abline(h=4, col="green", lwd = 2)
dev.off()

## Identify Potential Clusters =============================================

### WSS Plot ---------------------------------------------------------------
# Want to look for the bend in the curve
fviz_nbclust(datT2_mig, kmeans, method = "wss", nboot = 1000)
fviz_nbclust(datT2_res, kmeans, method = "wss", nboot = 1000)

png(file.path("1. Data","Figures","Migrant_WSS_Plot.png"), 
    width = 8, height = 6,
    units = "in", res = 600)
fviz_nbclust(datT2_mig, kmeans, method = "wss", nboot = 1000)
dev.off ()

png(file.path("1. Data","Figures","Resident_WSS_Plot.png"), 
    width = 8, height = 6,
    units = "in", res = 600)
fviz_nbclust(datT2_res, kmeans, method = "wss", nboot = 1000)
dev.off ()

### Silhouette Plot --------------------------------------------------------
# Want to get the highest mean sihouette 
png(file.path("1. Data","Figures","Migrant_Sil_Plot.png"), 
    width = 8, height = 6,
    units = "in", res = 600)
fviz_nbclust(datT2_mig, kmeans, method = "silhouette", nboot = 1000)
dev.off()

png(file.path("1. Data","Figures","Resident_Sil_Plot.png"), 
    width = 8, height = 6,
    units = "in", res = 600)
fviz_nbclust(datT2_res, kmeans, method = "silhouette", nboot = 1000)
dev.off()


### Gap Statistic ----------------------------------------------------------
gap_stat_mig <- cluster::clusGap(datT2_mig, FUN = kmeans, nstart = 25,
                             K.max = 10, B = 300)
gap_stat_res <- cluster::clusGap(datT2_res, FUN = kmeans, nstart = 25,
                                 K.max = 10, B = 300)

png(file.path("1. Data","Figures","Migrant_GapStat_Plot.png"), 
    width = 8, height = 6, units = "in", res = 300)
fviz_gap_stat(gap_stat_mig)
dev.off()

png(file.path("1. Data","Figures","Resident_GapStat_Plot.png"), 
    width = 8, height = 6, units = "in", res = 300)
fviz_gap_stat(gap_stat_res)
dev.off()

# Perform Clustering #######################################################
# 3 to 5 clusters has relatively good support:
# Dendrogram: 3 or 4
# WSS Plot: 3 to 5
# Silhouette Plot: 5
# Gap Statistic: 10 (lots of overlaps)

# K-means Cluster with k=3 
kmeans_mig <- eclust(datT2_mig, "kmeans", k = 2, nstart = 25, nboot = 10000, 
                 verbose = TRUE)
kmeans_res <- eclust(datT2_res, "kmeans", k = 3, nstart = 25, nboot = 10000, 
                     verbose = TRUE)

png(file.path("1. Data","Figures","Migrant_ClusterOutput.png"), 
    width = 8, height = 6,
    units = "in", res = 600)
fviz_cluster(kmeans_mig, stand = FALSE, geom = "point", axes = c(1,2))+
  theme_classic()
dev.off()

png(file.path("1. Data","Figures","Resident_ClusterOutput.png"), 
    width = 8, height = 6,
    units = "in", res = 600)
fviz_cluster(kmeans_res, stand = FALSE, geom = "point", axes = c(1,2))+
  theme_classic()
dev.off()

fviz_cluster(kmeans_mig, stand = FALSE, geom = "point", axes = c(1,3))+
  theme_classic() 
fviz_cluster(kmeans_mig, stand = FALSE, geom = "point", axes = c(1,4))+
  theme_classic() 

fviz_cluster(kmeans_res, stand = FALSE, geom = "point", axes = c(1,3))+
  theme_classic() 
fviz_cluster(kmeans_res, stand = FALSE, geom = "point", axes = c(1,4))+
  theme_classic() 


# Assign the resulting cluster ID to each individual
cluster_assign_mig <- data.frame("cluster" = factor(kmeans_mig$cluster,
                                                levels = sort(unique(kmeans_mig$cluster)))) %>%
  rownames_to_column("TagAge") %>%
  mutate(cluster = factor(paste("M",cluster, sep = ""), 
                          levels = paste("M",sort(unique(kmeans_mig$cluster)),sep ="")))
cluster_assign_res <- data.frame("cluster" = factor(kmeans_res$cluster,
                                                    levels = sort(unique(kmeans_res$cluster)))) %>%
  rownames_to_column("TagAge") %>%
  mutate(cluster = factor(paste("R",cluster, sep = ""), 
                          levels = paste("R",sort(unique(kmeans_res$cluster)),sep ="")))

cluster_assign <- bind_rows(cluster_assign_mig,cluster_assign_res)

assigned_dat <- dat %>%
  rownames_to_column("TagAge") %>%
  left_join(cluster_assign)

# library(plotly)
# p <- plot_ly(assigned_dat, x=~Comp.3, y=~Comp.2, 
#              z=~Comp.1, color=~cluster) %>%
#   add_markers(size=1.5)
# print(p)

# Interpretation ############################################################
## PERMANOVA ================================================================
datT2_assigned <- dat %>%
  decostand(method = "standardize") %>%
  rownames_to_column("TagAge") %>%
  left_join(cluster_assign) %>%
  select(-TagAge)

permanova <- vegan::adonis2(datT2_assigned[,1:(length(datT2_assigned)-2)] ~ cluster, 
                            data = datT2_assigned,
                            method = "euc",
                            na.action = na.omit)
permanova

## Pairwise PERMANOVA =======================================================

pairwise <- pairwise.adonis2(resp = vegdist(datT2_assigned[,1:(length(datT2_assigned)-2)]), 
                             fact = factor(datT2_assigned$cluster))

pairwise$p.value

## ANOVA ==================================================================
# Run anova on all the variables, all vars should be significant

b.sig = 0.05/length(dat)
anova<-aov(cbind(AgeBin, inside_sites_total, outside_sites_total, 
                 all_sites_total, inside_site_freq, outside_site_freq, 
                 all_site_freq, weekly_receivers_mean, weekly_receivers_max, 
                 CaptureAgeBin, time_btwn_detects_max, 
                 time_btwn_detects_mean, time_btwn_mvmts_max, 
                 time_btwn_mvmts_mean, mean_distance, q05_distance, q25_distance, 
                 q50_distance, q75_distance, q99_distance, emmigrant, 
                 mean_exit_woy, mean_entry_woy, entry_total, exit_total, all_total, 
                 entry_freq, exit_freq, all_freq, 
                 residence_total_INSIDE, residence_total_OUTSIDE, 
                 residence_mean_INSIDE, residence_mean_OUTSIDE, 
                 residence_max_INSIDE, residence_max_OUTSIDE, wyt, 
                 gate_perc_open
) ~ cluster, 
data = assigned_dat)

# get a summary of anova results
sum <- summary(anova) # All except q25 are significant

resp_df <- data.frame()
for(i in 1:length(sum)){
  response = stringr::str_split(names(sum[i])," ")[[1]][3]
  tmp <- sum[[i]][1,]
  tmp$response <- response
  resp_df[i,1:6] <- tmp[1,1:6]
}
(resp_df <- resp_df %>% select(response, Df:`Pr(>F)`) %>% arrange(`Pr(>F)`) %>%
    mutate(`Pr(>F)` = round(`Pr(>F)`,4)))
write.csv(resp_df, file.path("1. Data","Outputs","ANOVA_Response.csv"))

# View the means of each variable for each cluster
assigned_dat %>%
  group_by(cluster) %>%
  summarise(count = n(),
            emmigrant = sum(emmigrant),
            across(where(~is.numeric(.x)), ~mean(.x, na.rm = TRUE))) %>%
  t()

### Feature Importance ====================================================
# Identify those features which have the most weight in the clustering process

#devtools::install_github("o1iv3r/FeatureImpCluster")

# FeatureImpCluster requires input from flexclust so convert
res_mig <- flexclust::as.kcca(kmeans_mig, data = datT2_mig, k = 2)
res_res <- flexclust::as.kcca(kmeans_res, data = datT2_res, k = 3)

# Get the Feature Importance
nr_seeds <- 1000
seeds_vec <- sample(1:10000, nr_seeds)

savedImp_mig <- data.frame(matrix(0,nr_seeds,dim(datT2_mig)[2]))
count <- 1
for (s in seeds_vec) {
  set.seed(s)
  res <- flexclust::kcca(datT2_mig, k=2)
  set.seed(s)
  FeatureImp_res <- FeatureImpCluster::FeatureImpCluster(res,as.data.table(datT2_mig),
                                                         sub = 1,biter = 1)
  savedImp_mig[count,] <- FeatureImp_res$featureImp[sort(names(FeatureImp_res$featureImp))]
  count <- count + 1
}
names(savedImp_mig) <- sort(names(FeatureImp_res$featureImp))

savedImp_res <- data.frame(matrix(0,nr_seeds,dim(datT2_res)[2]))
count <- 1
for (s in seeds_vec) {
  set.seed(s)
  res <- flexclust::kcca(datT2_res, k=2)
  set.seed(s)
  FeatureImp_res <- FeatureImpCluster::FeatureImpCluster(res,as.data.table(datT2_res),
                                                         sub = 1,biter = 1)
  savedImp_res[count,] <- FeatureImp_res$featureImp[sort(names(FeatureImp_res$featureImp))]
  count <- count + 1
}
names(savedImp_res) <- sort(names(FeatureImp_res$featureImp))

savedImp_mig <- savedImp_mig  %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "misclassification") %>%
  mutate(group = "migrant")

savedImp_res <- savedImp_res  %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "misclassification") %>%
  mutate(group = "resident")

savedImp <- bind_rows(savedImp_mig, savedImp_res)

mig_imp_plot <- ggplot(savedImp_mig %>%
                         mutate(misclassification = ifelse(misclassification == 0, 0.00001,
                                                           misclassification)) %>%
                         group_by(group, variable) #%>%
                       #mutate(mean_importance = mean(importance))
                       , 
                       aes(x = reorder(variable, misclassification, decreasing = TRUE), 
                           y = misclassification, group = interaction(group, variable)))+
  geom_boxplot(outlier.shape = NA)+
  geom_smooth()+
  facet_grid(group ~ ., scales = "free")+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0)
  )+
  labs(x = "", y = "Rate of Misclassification")

res_imp_plot <- ggplot(savedImp_res %>%
         mutate(misclassification = ifelse(misclassification == 0, 0.00001,
                                           misclassification)) %>%
         group_by(group, variable) #%>%
       #mutate(mean_importance = mean(importance))
       , 
       aes(x = reorder(variable, misclassification, decreasing = TRUE), 
           y = misclassification, group = interaction(group, variable)))+
  geom_boxplot(outlier.shape = NA)+
  geom_smooth()+
  facet_grid(group ~ ., scales = "free")+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0)
  )+
  labs(x = "", y = "Rate of Misclassification")

png(file.path("1. Data","Figures","FeatureImportance.png"),
    width = 8, height = 6,  units = "in", res = 600)
print(gridExtra::grid.arrange(mig_imp_plot, res_imp_plot, nrow = 2))
dev.off()

## Examine the variables by cluster --------------------------------------

# Residency
residency <- ggplot(assigned_dat %>%
                      select(cluster, residence_total_INSIDE:residence_max_OUTSIDE) %>%
                      pivot_longer(cols = residence_total_INSIDE:residence_max_OUTSIDE,
                                   names_to = c("ignore","stat","location"),
                                   values_to = "value",
                                   names_sep = "_") %>%
                      mutate(Statistic = case_when(
                        grepl("mean", stat) ~ "Average",
                        grepl("max", stat) ~ "Maximum",
                        TRUE ~ "Total")))+
  geom_boxplot(aes(x = Statistic, y = value, color = location),
               linewidth = 0.75, fill = "grey75")+
  facet_grid(cluster ~ .)+
  scale_color_discrete(name = "Residency Location",
                       breaks = c("INSIDE","OUTSIDE"),
                       labels = c("Inside","Outside"))+
  scale_y_continuous(name = "Average Weeks of Residency per Stay")+
  
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        legend.position = "top",
        panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5))+
  coord_flip()
residency
ggsave(file.path("1. Data","Figures","Residency.png"), device = "png",
       width = 8, height = 6, units = "in",dpi = 600)

# Site Visits
sites_plot1 <- ggplot(assigned_dat %>%
                        select(cluster,TagAge, inside_site_freq:all_site_freq) %>%
                        pivot_longer(cols = inside_site_freq:all_site_freq,
                                     names_to = c("Site","ignore","stat"),
                                     values_to = "value",
                                     names_sep = "_"))+
  geom_boxplot(aes(y = value, x = stat, color = Site), fill = "grey75", linewidth = .75)+
  facet_grid(cluster ~ ., scales = 'free_y')+
  coord_flip()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5))+
  labs(y = "Frequency of Site Visits at Each Location", x = "")
sites_plot1
ggsave(file.path("1. Data","Figures","SiteVisitFrequency.png"), device = "png",
       width = 8, height = 6,        units = "in",dpi = 600)

# Distance
dist_df <- crossing(data.frame("x" = rep(c(0,100), each = 3),
                    "intercept" = rep(c(1.6,3.8,5.1),times = 2), 
                    "color" = c("green","orange","blue"),
                    "label" = c("1/2 Distance From Intake Canal to Radial Gates",
                                "Distance From Intake Canal to ORS1",
                                "Distance From Intake Canal to Closest Array Boundary")),
         "group" = c("M1","M2","R1","R2","R3"))

distance_plot <- ggplot(assigned_dat %>% 
                          select(cluster, q05_distance:q99_distance) %>%
                          pivot_longer(cols = q05_distance:q99_distance, 
                                       names_to = "metric",
                                       values_to = "distance") %>%
                          mutate(metric = as.numeric(stringr::str_sub(metric,2,3))),
                        aes(x = metric, y = distance, group = cluster))+
  stat_summary(fun.data=mean_se,geom="linerange",color="black", linewidth = 0.75)+
  stat_summary(fun = "mean", geom="line", linewidth = 0.75)+
  geom_line(data = dist_df,
            aes(x = x, y = intercept, color = color),
             linetype = "dashed", linewidth = 0.75, inherit.aes = FALSE)+
  facet_grid(cluster ~ .)+
  scale_x_continuous(name = "Quantile",
                     breaks = c(05, 25, 50, 75, 99),
                     labels = c("5%","25%","50%","75%","99%"),
                     limits = c(0,100))+
  scale_y_continuous(name = "Mean Weekly Distance Travelled (km)")+
  scale_color_discrete(name = "Distance",
                       breaks = c("green","orange","blue"),
                       labels= c("1/2 Distance From Intake Canal to Radial Gates",
                                   "Distance From Intake Canal to ORS1",
                                   "Distance From Intake Canal to Closest Array Boundary"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        panel.grid.major.y = element_line(color = "grey90"),
        panel.background = element_rect(color = "grey50"),
        legend.position = "bottom",
        legend.key.width = unit(2,"lines"),
        legend.title = element_blank(),
        legend.justification = "center")+
  guides(color = guide_legend(nrow = 2,
                              byrow = TRUE))
distance_plot
ggsave(file.path("1. Data","Figures","WeeklyDistance.png"), device = "png",
       width = 8, height = 6,        units = "in",dpi = 600)

# Average Weekly Site Visits
sites_plot2 <- ggplot(assigned_dat %>%
                        select(cluster,weekly_receivers_mean,weekly_receivers_max) %>%
                        pivot_longer(cols = weekly_receivers_mean:weekly_receivers_max,
                                     names_to = c("ignore_1","ignore_2","stat"),
                                     names_sep = "_",
                                     values_to = "value"))+
  geom_boxplot(aes(y = value, x = stat, group = stat), linewidth = 0.75, 
               fill = "grey75")+
  facet_grid(cluster ~ ., scales = 'free_y')+
  coord_flip()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5))+
  labs(y = "Average Number of Sites Visits per Week",
       x = "")
sites_plot2
ggsave(file.path("1. Data","Figures","AvgWeeklySites.png"), device = "png",
       width = 8, height = 6,        units = "in",dpi = 600)

# Number and Frequency of Movements
movements <- ggplot(assigned_dat %>%
                      select(cluster,entry_freq:all_freq) %>%
                      pivot_longer(cols = entry_freq:all_freq,
                                   names_to = c("Direction","stat"),
                                   names_sep = "_",
                                   values_to = "value") %>%
                      mutate(value = ifelse(stat == "freq", value*10, value)))+
  geom_boxplot(aes(x = stat, y = value,
                   color = Direction), linewidth = 0.75, fill = "grey75"
  )+
  facet_grid(cluster ~ .)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        legend.position = "top",
        panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5)
  )+
  scale_y_continuous(name = "Frequency of Transits (transits/week)")+
  scale_fill_discrete(name = "")+
  coord_flip()
movements

ggsave(file.path("1. Data","Figures","Movements.png"), device = "png",
       width = 8, height = 6,        units = "in",dpi = 600)

#Exit and Entry timing
seasons <- data.frame(
  season = c("Immigration_a","Spawn","Emmigration","Residence","Immigration_b"),
  start =  c(ymd("2018-01-01"),ymd("2018-03-15"),ymd("2018-05-15"),ymd("2018-08-01"), ymd("2018-11-01")),
  end = c(ymd("2018_03-15"),ymd("2018-05-15"),ymd("2018-08-01"),ymd("2018-11-01"),ymd("2019-01-01"))) %>%
  mutate(start_doy = yday(start),
         end_doy = yday(end),)

seasons_timing <- seasons %>%
  crossing(x = c(-2,2))

mean_exit_woy <- ggplot()+
  geom_ribbon(data = seasons_timing, aes(x = x, 
                                  ymin = start, 
                                  ymax = end,
                                  fill = season),
              alpha = 0.5)+
  geom_violin(data= assigned_dat%>%
                separate(TagAge,into = c("TagID","AgeBin"), sep = "-",
                         extra = "merge") %>%
                filter(emmigrant == TRUE & exit_total >= 1),
              aes(y = lubridate::ymd("2017-12-31")+(mean_exit_woy*7), x = 1), 
              fill = "grey50",
              trim = TRUE, bw = 7)+
  geom_text(data = assigned_dat%>%
              separate(TagAge,into = c("TagID","AgeBin"), sep = "-",
                       extra = "merge") %>%
              filter(emmigrant == TRUE & exit_total >= 1) %>%
              select(TagID, mean_exit_woy, cluster) %>%
              distinct() %>%
              group_by(cluster) %>%
              summarise(mean = mean(mean_exit_woy),
                        count = n()),
            aes(y = lubridate::ymd("2017-12-31")+(mean*7), x = 1.25, 
                label = paste("N =",count, sep = " ")),
            inherit.aes = FALSE)+
  scale_fill_manual(name = "Season",
                    breaks = c("Immigration_a","Spawn","Emmigration",
                    "Residence","Immigration_b"),
                    values = c("olivedrab1","steelblue2","tomato2",
                              "goldenrod","olivedrab1"),
                    labels = c("Immigration","Spawn","Emmigration",
                               "Residence","Immigration"))+
  facet_grid(cluster ~ .,drop = FALSE)+
  coord_flip(xlim =c(0.5,1.5), expand = FALSE)+
  theme_classic()+
  scale_y_date(name = "Date", date_labels = "%b",
               breaks = "1 month")+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5))+
  labs(title = "Mean Week of Exit", x = "")
mean_exit_woy

ggsave(file.path("1. Data","Figures","ExitWeek.png"), device = "png",
       width = 8, height = 6,        units = "in",dpi = 600)

mean_entry_woy <- ggplot()+
  geom_ribbon(data = seasons_timing, aes(x = x, 
                                  ymin = start, 
                                  ymax = end,
                                  fill = season),
              alpha = 0.5)+
  geom_violin(data= assigned_dat%>%
                separate(TagAge,into = c("TagID","AgeBin"), sep = "-",
                         extra = "merge") %>%
                filter(emmigrant == TRUE & entry_total >= 1),
              aes(y = lubridate::ymd("2017-12-31")+(mean_entry_woy*7), x = 1), 
              fill = "grey50",
              trim = TRUE, bw = 7)+
  geom_text(data = assigned_dat%>%
              separate(TagAge,into = c("TagID","AgeBin"), sep = "-",
                       extra = "merge") %>%
              filter(emmigrant == TRUE & entry_total >= 1) %>%
              select(TagID, mean_entry_woy, cluster) %>%
              distinct() %>%
              group_by(cluster) %>%
              summarise(mean = mean(mean_entry_woy),
                        count = n()),
            aes(y = lubridate::ymd("2017-12-31")+(mean*7), x = 1.25, 
                label = paste("N =",count, sep = " ")),
            inherit.aes = FALSE)+
  facet_grid(cluster ~ .,drop = FALSE)+
  coord_flip(xlim =c(0.5,1.5), expand = FALSE)+
  theme_classic()+
  scale_fill_manual(name = "Season",
                    breaks = c("Immigration_a","Spawn","Emmigration",
                               "Residence","Immigration_b"),
                    values = c("olivedrab1","steelblue2","tomato2",
                               "goldenrod","olivedrab1"),
                    labels = c("Immigration","Spawn","Emmigration",
                               "Residence","Immigration"))+
  scale_y_date(name = "Date", date_labels = "%b-%d")+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5))+
  labs(title = "Mean Week of Entry")
mean_entry_woy
ggsave(file.path("1. Data","Figures","EntryWeek.png"), device = "png",
       width = 8, height = 6,        units = "in",dpi = 600)

### Time between Detects Plot ---------------------------------------------
time_detect_plot <- ggplot(assigned_dat %>%
                             select(cluster,TagAge,time_btwn_detects_mean:time_btwn_detects_max) %>%
                             pivot_longer(cols = time_btwn_detects_mean:time_btwn_detects_max,
                                          names_to = "Statistic",
                                          values_to = "value"))+
  geom_boxplot(aes(y = value, x = 1, color = Statistic), linewidth = 0.75,
               fill = "grey75")+
  geom_text(
    aes(label = paste(after_stat(count/2),"Detections"),
        x = 2,
        y = after_stat(100)),
    stat = "count",
    #position = position_nudge(x = -7, y = )
  )+
  facet_grid(cluster ~ ., scales = 'free_y')+
  coord_flip()+
  theme_classic()+
  theme(panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5))+
  scale_x_discrete(breaks = c("time_btwn_detects_mean","time_btwn_detects_max"),
                   labels = c("Mean","Max"))+
  scale_color_discrete(breaks = c("time_btwn_detects_mean","time_btwn_detects_max"),
                       labels = c("Mean","Max"))+
  labs(x = "", y = "Weeks between Detections")
time_detect_plot
ggsave(file.path("1. Data","Figures","TimeBtwnDetects.png"), device = "png",
       width = 8, height = 6, units = "in",dpi = 600)

### Time between Movements Plot --------------------------------------------
time_move_plot <- ggplot(assigned_dat %>%
                           select(cluster,TagAge,time_btwn_mvmts_mean,
                                  time_btwn_mvmts_max) %>%
                           pivot_longer(cols = c("time_btwn_mvmts_mean",
                                                 "time_btwn_mvmts_max"),
                                        names_to = "Statistic",
                                        values_to = "value") %>%
                           filter(value != 100))+
  geom_boxplot(aes(y = value, x = factor(Statistic), color = Statistic), linewidth = 0.75,
               fill = "grey75")+ 
  geom_text(
    aes(label = paste(after_stat(count/2),"Movements"),
        x = 2,
        y = after_stat(70)),
    stat = "count",
    #position = position_nudge(x = -7, y = )
  )+
  facet_grid(cluster ~ ., scales = 'free_y', drop = FALSE)+
  coord_flip()+
  theme_classic()+
  theme(panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5))+
  scale_x_discrete(breaks = c("time_btwn_mvmts_mean",
                              "time_btwn_mvmts_max"),
                   labels = c("Mean", "Max"))+
  scale_color_discrete(breaks = c("time_btwn_mvmts_mean",
                                  "time_btwn_mvmts_max"),
                       labels = c("Mean", "Max"))+
  labs(x = "", y = "Weeks between Transits")
time_move_plot
ggsave(file.path("1. Data","Figures","TimeBtwnMovements.png"), device = "png",
       width = 8, height = 6,        units = "in",dpi = 600)

### Age Plot ---------------------------------------------------------------
age_plot <- ggplot(assigned_dat %>%
                     select(cluster, AgeBin) %>%
                     pivot_longer(cols = AgeBin,
                                  names_to = "event",values_to = "age_bin") %>%
                     mutate(event = case_when(
                       event == "FirstExitAgeBin" ~ "First Exit",
                       event == "CaptureAgeBin" ~ "Capture",
                       event == "AgeBin" ~ "Current"
                     )) %>%
                     filter(event == "Current"))+
  geom_bar(aes(x = "",fill = factor(age_bin,
                                    levels = c(1:4)), group = age_bin), 
           stat = "count", position = "fill")+
  facet_grid(cluster ~ .)+
  coord_flip()+
  scale_fill_discrete(name = "Age Bin", breaks = c(1:4),
                      labels = c("0","1-2","3-5","6+"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5))+
  theme(legend.position = "top")
age_plot
ggsave(file.path("1. Data","Figures","Cluster_08_Age.png"), device = "png",
       width = 8, height = 6,        units = "in",dpi = 600)

### Water Year Type Plot ---------------------------------------------------
wyt_weights <- data.frame(year = c(2013,2014,2015,2016,2017,2018)) %>%
  mutate(wyt = get_water_year_type(year),
         duration = 52) %>%
  group_by(wyt) %>%
  summarise(duration = sum(duration)) %>%
  mutate(wyt = factor(wyt, levels =  c("Critical","Dry","Below Normal","Above Normal",
                                       "Wet")))
  


wyt_plot <- ggplot(assigned_dat %>%
                     select(cluster, wyt) %>%
                     mutate(wyt = factor(wyt, levels = c(1:5),
                                         labels =  c("Critical","Dry","Below Normal","Above Normal",
                                                     "Wet"))) %>%
                     left_join(wyt_weights) %>%
                     group_by(cluster, wyt) %>%
                     summarise(count = n(),
                               duration = mean(duration),
                               adj_count = count/duration))+
  geom_bar(aes(x = "", y = adj_count, fill = reorder(as.character(wyt), adj_count), group = wyt), 
           stat = "identity", position = "fill")+
  facet_grid(cluster ~ .)+
  coord_flip()+
  scale_fill_discrete(name = "Water Year Type", breaks = c("Critical","Dry","Below Normal","Above Normal",
                                 "Wet"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5))+
  theme(legend.position = "top")
wyt_plot

### Season Plot-----------------------------------------------------------------
season_weights <- seasons %>%
  distinct(season, start,end) %>%
  mutate(duration = (end-start)/7) %>%
  group_by(season) %>%
  summarise(weight = round(sum(duration),0))

season_plot <- ggplot(assigned_dat %>%
                        group_by(cluster,season) %>%
                        tally(name = "count") %>%
                        mutate(season = factor(season,
                                               levels = c(1:4),
                                               labels = c("Immigration","Spawn",
                                                          "Emmigration",
                                                          "Residence"))) %>%
                        left_join(season_weights) %>%
                        mutate(adj_count = count/weight)
)+
  geom_bar(aes(x = "", y = adj_count, fill = reorder(season,adj_count)), 
           stat = "identity", position = "fill")+
  facet_grid(cluster ~ .)+
  coord_flip()+
  scale_fill_discrete(name = "season", 
                      breaks = c("Immigration","Spawn","Emmigration","Residence"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5),
        legend.position = "top")
season_plot
ggsave(file.path("1. Data","Figures","Season.png"), device = "png",
       width = 8, height = 6,        units = "in",dpi = 600)

### Detection Histories ----------------------------------------------------
TagBKM_Bin <- readRDS(file.path("1. Data", "Outputs", "TagBKM_Bin.rds"))
seasons %>%
  select(season,start_doy, end_doy) %>%
  crossing(year = c(2012:2017)) %>%
  mutate(start_date = ymd(paste(year,"12","31"))+days(start_doy),
         end_date = as_date(ifelse(start_doy > 300,
                           ymd(paste(year+2,"1","1")),
                           ymd(paste(year,"12","31"))+days(end_doy)))) %>%
  mutate(season_year = paste(season,year,sep = "-")) %>%
  crossing(y = c(0,16)) %>%
  arrange(start_date)-> all_seasons

assigned_dat%>%
  separate(TagAge,into = c("TagID","AgeBin"), sep = "-",
           extra = "merge") %>%
  select(TagID, AgeBin, cluster) %>%
  separate("AgeBin", into = c("AgeBin","extra"), sep = "-",
           extra = "merge") %>%
  group_by(TagID,AgeBin) %>%
  select(-extra) %>%
  mutate(AgeBin = case_when(
    AgeBin == "1" ~ "1-2",
    AgeBin == "3" ~ "3-5",
    TRUE ~ "6+"
  )) %>%
  summarise(cluster_list = list(as.character(cluster)))%>%
  mutate(primary_cluster = getmode(unlist(cluster_list)),
         secondary_cluster = getmode(unlist(cluster_list)[which(unlist(cluster_list) != primary_cluster)])) %>%
  mutate(primary_cluster = as.character(primary_cluster),
         secondary_cluster = as.character(secondary_cluster)) ->tag_clusters

for(i in (tag_clusters %>% pull(var = "TagID"))){
  tmp <- TagBKM_Bin %>%
    filter(TagID == i) %>%
    filter(SiteCode != "Tag Failure") %>%
    mutate(TagID = as.character(TagID)) %>%
    left_join(tag_clusters)
  
  if(nrow(tmp)>1){
    ggplot()+
      geom_ribbon(data = all_seasons %>%
                    mutate(y = ifelse(y == 16, length(unique(tmp$SiteCode))+1, y)),
                  aes(y = y,
                  xmin = start_date,
                  xmax = end_date,
                  fill = season,
                  group = interaction(season,year)),
                  alpha = 0.35)+
      geom_point(data = tmp,
                 aes(x = studyweek_startdate(Week), y = SiteCode, group = TagID))+
      geom_step(data = tmp,
                aes(x = studyweek_startdate(Week), y = SiteCode, group = TagID))+
      geom_hline(yintercept = as.numeric(factor("RGD1",
                                                levels = levels(tmp$SiteCode))) + 0.5,
                                                color = "blue", 
                 linewidth = 2, linetype = "dashed")+
      scale_x_date(name = "Date",
                   breaks = scales::breaks_pretty(n = 8),
                   date_labels = "%b-%d")+
      scale_y_discrete()+
      scale_fill_manual(name = "Season",
                        breaks = c("Immigration_a","Spawn","Emmigration",
                                   "Residence","Immigration_b"),
                        values = c("olivedrab1","steelblue2","tomato2",
                                   "goldenrod","olivedrab1"),
                        labels = c("Immigration","Spawn","Emmigration",
                                   "Residence","Immigration"))+
      theme_classic()+
      theme(panel.spacing = unit(0, "lines"),
            axis.text.x = element_text(angle = -45, hjust = 0),
            axis.title.y = element_blank())+
      labs(subtitle = paste("Primary Cluster:",tmp$primary_cluster[1]))+
      coord_cartesian(xlim = c(studyweek_startdate(min(tmp$Week)-1),
                               studyweek_startdate(max(tmp$Week)+1)))
    ggsave(file.path("1. Data","Figures","DetectionPlots",paste(tmp$TagID[1],tmp$primary_cluster[1],"ex.png",sep = "-")),
           device = "png", width = 10, height = 5, dpi = 600)
  } else {
    print("Tag not in cluster")
  }
}
