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
datT1 <- decostand(dat, method = "standardize")

# Exploration ==============================================================

## Correlation plot ========================================================
corr = cor(datT1)
testRes = corrplot::cor.mtest(datT1, conf.level = 0.95)
corrplot::corrplot(corr, p.mat = testRes$p, method = 'circle', type = 'lower', 
                   insig='blank',
                   #addCoef.col ='black', 
                   number.cex = 0.2, order = 'AOE', diag=FALSE)

corr <- as_tibble(corr)
corr <- corr %>%
  mutate(V1 = names(corr)) %>%
  pivot_longer(cols = AgeBin:gate_perc_open,
               names_to = "V2",values_to = "corr") %>%
  filter(V1 != V2) %>%
  arrange(desc(corr))

# Remove Variables with >95% correlation (same variable only)
datT2 <- dat %>%
  select(-time_btwn_detects_max, -time_btwn_mvmts_max, -residence_max_OUTSIDE,
  -CaptureAgeBin) %>%
  decostand(method = "standardize")

## Principal Components ====================================================
fit <- prcomp(datT2)
summary(fit) # print variance accounted for
plot(fit,type="lines") # scree plot
biplot(fit, cex = c(0.1,.75), scale = 1) #PCA Biplot

png(file.path("1. Data","Figures","PCA_biplot.png"), width = 10, height = 10,
    units = "in", res = 600)
biplot(fit, cex = c(0.1,.75), scale = 1) #PCA Biplot
dev.off()

## Dendogram ===============================================================
datEuc <- vegdist(datT2, method = "euclidean", na.rm = TRUE)
euc_agnes <- cluster::agnes(datEuc, method = "ward")

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

dg <- as.dendrogram(euc_agnes)
dg <- dendextend::remove_nodes_nodePar(dg)
dg_season <- dendrapply(dg, labelCol_season)
dg_age <- dendrapply(dg, labelCol_age)

png(file.path("1. Data","Figures","Dendrogram.png"), 
    width = 12, height = 6,
    units = "in", res = 600)
par(mar = c(0.1,3,0.1,0.11))
plot(dg_age)
dendextend::rect.dendrogram(dg, k = 5)
#abline(h=5, col="blue", lwd = 2)
#abline(h=4, col="green", lwd = 2)
dev.off()

## Identify Potential Clusters =============================================

### WSS Plot ---------------------------------------------------------------
# Want to look for the bend in the curve
png(file.path("1. Data","Figures","WSS_Plot.png"), 
    width = 8, height = 6,
    units = "in", res = 600)
fviz_nbclust(datT2, kmeans, method = "wss", nboot = 1000)
dev.off ()

### Silhouette Plot --------------------------------------------------------
# Want to get the highest mean sihouette 
png(file.path("1. Data","Figures","Sil_Plot.png"), 
    width = 8, height = 6,
    units = "in", res = 600)
fviz_nbclust(datT2, kmeans, method = "silhouette", nboot = 1000)
dev.off()

### Gap Statistic ----------------------------------------------------------
gap_stat <- cluster::clusGap(datT2, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 300)
png(file.path("1. Data","Figures","GapStat_Plot.png"), 
    width = 8, height = 6, units = "in", res = 300)
fviz_gap_stat(gap_stat)
dev.off()

# Perform Clustering #######################################################
# 3 to 5 clusters has relatively good support:
# Dendrogram: 3 or 4
# WSS Plot: 3 to 5
# Silhouette Plot: 5
# Gap Statistic: 10 (lots of overlaps)

# K-means Cluster with k=3 
kmeans <- eclust(datT2, "kmeans", k = 3, nstart = 25, nboot = 10000, 
                   verbose = TRUE)

png(file.path("1. Data","Figures","ClusterOutput.png"), 
    width = 8, height = 6,
    units = "in", res = 600)
fviz_cluster(kmeans, stand = FALSE, geom = "point", axes = c(1,2))+
  theme_classic()
dev.off()
fviz_cluster(kmeans, stand = FALSE, geom = "point", axes = c(1,3))+
  theme_classic() 
fviz_cluster(kmeans, stand = FALSE, geom = "point", axes = c(1,4))+
  theme_classic() 

# Assign the resulting cluster ID to each individual
cluster_assign <- data.frame("cluster" = factor(kmeans$cluster,
                                                levels = sort(unique(kmeans$cluster)))) %>%
  rownames_to_column("TagAge")

assigned_dat <- dat %>%
  bind_cols(cluster_assign) %>%
  left_join(as.data.frame(fit$scores) %>% rownames_to_column(var = "TagAge"))

# library(plotly)
# p <- plot_ly(assigned_dat, x=~Comp.3, y=~Comp.2, 
#              z=~Comp.1, color=~cluster) %>%
#   add_markers(size=1.5)
# print(p)

# Interpretation ############################################################
## PERMANOVA ================================================================
datT2_assigned <- datT2 %>%
  bind_cols(cluster_assign)

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

b.sig = 0.05/length(datT2)
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
res <- flexclust::as.kcca(kmeans, data = datT2, k = 4)

# Get the Feature Importance
nr_seeds <- 1000
seeds_vec <- sample(1:10000, nr_seeds)

savedImp <- data.frame(matrix(0,nr_seeds,dim(datT2)[2]))
count <- 1
for (s in seeds_vec) {
  set.seed(s)
  res <- flexclust::kcca(datT2, k=4)
  set.seed(s)
  FeatureImp_res <- FeatureImpCluster::FeatureImpCluster(res,as.data.table(datT2),
                                                         sub = 1,biter = 1)
  savedImp[count,] <- FeatureImp_res$featureImp[sort(names(FeatureImp_res$featureImp))]
  count <- count + 1
}
names(savedImp) <- sort(names(FeatureImp_res$featureImp))

savedImp <- savedImp  %>%
  pivot_longer(cols = everything(), 
               names_to = "variable",
               values_to = "misclassification")

savedImp %>%  
  group_by(variable) %>%
  summarise(mean_imp = mean(misclassification)) %>%
  arrange(desc(mean_imp)) %>%
  pull(variable) -> var_orders

ggplot(savedImp %>%
         mutate(misclassification = ifelse(misclassification == 0, 0.0001,
                                    misclassification)) %>%
         group_by(variable) #%>%
         #mutate(mean_importance = mean(importance))
       , 
       aes(x = factor(variable,levels = var_orders), y = misclassification))+
  geom_boxplot(outlier.shape = NA)+
  geom_smooth()+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0)
  )+
  labs(x = "", y = "Rate of Misclassification")
ggsave(file.path("1. Data","Figures","FeatureImportance.png"), device = "png",
       width = 8, height = 6,
       units = "in", dpi = 600)

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
                        select(cluster,TagAge, inside_sites_total:all_site_freq) %>%
                        pivot_longer(cols = inside_sites_total:all_site_freq,
                                     names_to = c("Site","ignore","stat"),
                                     values_to = "value",
                                     names_sep = "_"))+
  geom_boxplot(aes(y = value, x = stat, color = Site), fill = "grey75", linewidth = .75)+
  facet_grid(cluster ~ ., scales = 'free_y')+
  coord_flip()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5))+
  labs(y = "Total Number of Site Visits at Each Location", x = "")
sites_plot1
ggsave(file.path("1. Data","Figures","TotalSiteVisits.png"), device = "png",
       width = 8, height = 6,        units = "in",dpi = 600)

# Distance
distance_plot <- ggplot(assigned_dat %>% 
                          select(cluster, q05_distance:q99_distance) %>%
                          pivot_longer(cols = q05_distance:q99_distance, 
                                       names_to = "metric",
                                       values_to = "distance") %>%
                          mutate(metric = as.numeric(stringr::str_sub(metric,2,3))),
                        aes(x = metric, y = distance, group = cluster)
                        )+
  stat_summary(fun.data=mean_se,geom="linerange",color="black", linewidth = 1)+
  stat_summary(fun = "mean", geom="line", linewidth = 01)+
  geom_hline(aes(yintercept = 1.6), linetype = "dashed", linewidth = 0.5,
             color = "blue",alpha = 0.5)+
  geom_hline(aes(yintercept = 3.8), linetype = "dashed", linewidth = 0.5,
             color = "blue",alpha = 0.5)+
  
  geom_text(aes(x = 25, y = 1.9, label = "1/2 Distance From Intake Canal to Radial Gates"),
            hjust = 0)+
  geom_hline(aes(yintercept = 3.8), linetype = "dashed", linewidth = 0.5,
             color = "blue",alpha = 0.5)+
  geom_text(aes(x = 25, y = 4.1, label = "Distance From Intake Canal to ORS1"),
            hjust = 0)+
  geom_hline(aes(yintercept = 5.1), linetype = "dashed", linewidth = 0.5,
             color = "blue",alpha = 0.5)+
  geom_text(aes(x = 25, y = 5.4, label = "Distance From Intake Canal to Closest Array Boundary"),
            hjust = 0)+
  facet_grid(cluster ~ .)+
  scale_x_continuous(name = "Quantile",
                   breaks = c(05, 25, 50, 75, 99),
                   labels = c("5%","25%","50%","75%","99%"),
                   limits = c(5,100))+
  scale_y_continuous(name = "Mean Weekly Distance Travelled")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        panel.grid.major.y = element_line(color = "grey90"),
        panel.background = element_rect(color = "grey50"))
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
                      select(cluster,entry_total:all_freq) %>%
                      pivot_longer(cols = entry_total:all_freq,
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
  scale_y_continuous(name = "Number of Transits")+
  scale_fill_discrete(name = "")+
  coord_flip()
movements

ggsave(file.path("1. Data","Figures","Movements.png"), device = "png",
       width = 8, height = 6,        units = "in",dpi = 600)

#Exit and Entry timing
mean_exit_woy <- ggplot()+
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
  facet_grid(cluster ~ .,drop = FALSE)+
  coord_flip()+
  theme_classic()+
  scale_y_date(name = "Date", date_labels = "%b-%d")+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5))+
  labs(title = "Mean Week of Exit")
mean_exit_woy

ggsave(file.path("1. Data","Figures","ExitWeek.png"), device = "png",
       width = 8, height = 6,        units = "in",dpi = 600)

mean_entry_woy <- ggplot()+
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
  coord_flip()+
  theme_classic()+
  scale_y_date(name = "Date", date_labels = "%b-%d")+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5))+
  labs(title = "Mean Week of Entry")
mean_entry_woy
ggsave(file.path("1. Data","Figures","EntryWeek.png"), device = "png",
       width = 8, height = 6,        units = "in",dpi = 600)

# Time between Detects
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
       width = 8, height = 6,        units = "in",dpi = 600)

# Emmigrant Status
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
        legend.position = "top",
        panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5))+
  labs(x = "", y = "Number of Individuals")
emmigrant_plot
ggsave(file.path("1. Data","Figures","Emmigration.png"), device = "png",
       width = 8, height = 6,        units = "in",dpi = 600)

# Time between Movements
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
  facet_grid(cluster ~ ., scales = 'free_y')+
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

# Age
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

for(i in (assigned_dat%>%
                separate(TagAge,into = c("TagID","AgeBin"), sep = "-",
                         extra = "merge") %>%
                select(TagID) %>% pull())){
  tmp <- TagBKM_Bin %>%
    filter(TagID == i) %>%
    filter(SiteCode != "Tag Failure") %>%
    mutate(TagID = as.character(TagID)) %>%
    left_join(assigned_dat%>%
                separate(TagAge,into = c("TagID","AgeBin"), sep = "-",
                         extra = "merge"))

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
      facet_grid(TagID ~.)+
      labs(subtitle = paste("Cluster:",tmp$cluster[1]))
    ggsave(file.path("1. Data","Figures","DetectionPlots",paste(tmp$TagID[1],tmp$cluster[1],"ex.png",sep = "-")),
           device = "png", width = 5, height = 2.5)
  } else {
    print("Tag not in cluster")
  }
}
