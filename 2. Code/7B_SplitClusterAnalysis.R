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
dat_trns <- dat %>% 
  filter(residence_total_OUTSIDE > 0 |
           outside_sites_total > 0)

datT1_trns <- dat_trns %>%
  decostand(method = "standardize") %>%
  select(-ccf_cross, - wyt, -gate_perc_open)

dat_res <- dat %>% 
  filter(residence_total_OUTSIDE == 0 &
           outside_sites_total == 0) 

datT1_res <- dat_res %>%
  decostand(method = "standardize") %>%
  select(-outside_sites_total, -all_sites_total, -exit_total, -entry_total, 
         -all_total, -residence_mean_OUTSIDE,-residence_mean_INSIDE, -residence_total_OUTSIDE, 
         -mean_exit_woy, -mean_entry_woy, -time_btwn_mvmts_mean, 
         -ccf_cross, -wyt, -gate_perc_open)

datT1_all <-dat %>%
  decostand(method = "standardize") %>%
  select(-ccf_cross, -wyt, -gate_perc_open)

# Exploration ==============================================================

## Correlation plot ========================================================
corr_trns = cor(datT1_trns)
corr_res = cor(datT1_res)
corr_all = cor(datT1_all)

testRes_trns = corrplot::cor.mtest(datT1_trns, conf.level = 0.95)
testRes_res = corrplot::cor.mtest(datT1_res, conf.level = 0.95)
testRes_all = corrplot::cor.mtest(datT1_all, conf.level = 0.95)

corrplot::corrplot(corr_trns, p.mat = testRes_trns$p, method = 'circle', type = 'lower', 
                   insig='blank',
                   #addCoef.col ='black', 
                   number.cex = 0.2, order = 'AOE', diag=FALSE)

corrplot::corrplot(corr_res, p.mat = testRes_res$p, method = 'circle', type = 'lower', 
                   insig='blank',
                   #addCoef.col ='black', 
                   number.cex = 0.2, order = 'AOE', diag=FALSE)

corr_trns <- as_tibble(corr_trns)
corr_res <- as_tibble(corr_res)
corr_trns <- corr_trns %>%
  mutate(V1 = names(corr_trns)) %>%
  pivot_longer(cols = inside_sites_total:residence_mean_OUTSIDE,
               names_to = "V2",values_to = "corr") %>%
  filter(V1 != V2) %>%
  arrange(desc(corr))
View(corr_trns)

corr_res <- corr_res %>%
  mutate(V1 = names(corr_res)) %>%
  pivot_longer(cols = inside_sites_total:residence_total_INSIDE,
               names_to = "V2",values_to = "corr") %>%
  filter(V1 != V2) %>%
  arrange(desc(corr))
View(corr_res)

# Remove Variables with 99% correlation (same variable only)
# NONE!
datT2_trns <- datT1_trns
datT2_res <- datT1_res
datT2_all <- datT1_all

# Perform Clustering #######################################################
# K-means Clustering using gap statistic to identify optimal clusters
kmeans_trns <- eclust(datT2_trns, "kmeans", k.max = 10, nboot = 500, verbose = TRUE)
kmeans_res <- eclust(datT2_res, "kmeans", k.max = 10, nboot = 500, verbose = TRUE)
kmeans_all <- eclust(datT2_all, "kmeans", k.max = 10, nboot = 500, verbose = TRUE)

# Jaccard Bootstrap Cluster Validation
library(fpc)
cboot_trns <- fpc::clusterboot(datT2_trns, B = 1000, clustermethod = kmeansCBI,
                               k = 3, count = FALSE, seed = 123)
cboot_res <- fpc::clusterboot(datT2_res, B = 1000, clustermethod = kmeansCBI,
                              k = 4, count = FALSE, seed = 123)
cboot_all <- fpc::clusterboot(datT2_all, B = 1000, clustermethod = kmeansCBI,
                              k = 4, count = FALSE, seed = 123)

# View the means of the Jaccard Bootstrap Validation
# >0.85 = Highly Stable, >0.75 = Valid and stable, >0.6 = patterns exist but
# classification is questionable, <0.6 invalid
(trns_bootmean <- cboot_trns$bootmean)
(res_bootmean <- cboot_res$bootmean)
(all_bootmean <- cboot_all$bootmean)

kmeans_trns <- eclust(datT2_trns, "kmeans", k = 3, nboot = 1000, verbose = TRUE,
                      seed = 123)
kmeans_res <- eclust(datT2_res, "kmeans", k = 4, nboot = 1000, verbose = TRUE,
                     seed = 123)
kmeans_all <- eclust(datT2_all, "kmeans", k = 4, nboot = 1000, verbose = TRUE,
                     seed = 123)

png(file.path("1. Data","Figures","Transient_ClusterOutput.png"), 
    width = 8, height = 6,
    units = "in", res = 600)
fviz_cluster(kmeans_trns, stand = FALSE, geom = "point", axes = c(1,2))+
  theme_classic()
dev.off()

png(file.path("1. Data","Figures","Resident_ClusterOutput.png"), 
    width = 8, height = 6,
    units = "in", res = 600)
fviz_cluster(kmeans_res, stand = FALSE, geom = "point", axes = c(1,2))+
  theme_classic()
dev.off()

# Check the Silhouette Width
(sil_trns <- fviz_nbclust(datT2_trns, kmeans, method = "silhouette", nboot = 1000))
(sil_res <- fviz_nbclust(datT2_res, kmeans, method = "silhouette", nboot = 1000))
png(file.path("1. Data","Figures","Transient_Sil_Plot.png"),
    width = 8, height = 6,
    units = "in", res = 600)
sil_trns
dev.off()

png(file.path("1. Data","Figures","Resident_Sil_Plot.png"),
    width = 8, height = 6,
    units = "in", res = 600)
sil_res
dev.off()
# Silhouette Width supports 2 clusters for residents, good support for 2-8 clusters

# Assign the resulting cluster ID to each 
cluster_assign_trns <- data.frame("cluster" = factor(cboot_trns$result$result$cluster,
                                                levels = sort(unique(cboot_trns$result$result$cluster)))) %>%
  rownames_to_column("TagAge") %>%
  mutate(bootmean = trns_bootmean[cluster]) %>%
  mutate(cluster = factor(paste("T",cluster, sep = ""), 
                          levels = paste("T",sort(unique(cboot_trns$result$result$cluster)),sep ="")))

cluster_assign_res <- data.frame("cluster" = factor(cboot_res$result$result$cluster,
                                                    levels = sort(unique(cboot_res$result$result$cluster)))) %>%
  rownames_to_column("TagAge") %>%
  mutate(bootmean = res_bootmean[cluster]) %>%
  mutate(cluster = factor(paste("R",cluster, sep = ""), 
                          levels = paste("R",sort(unique(cboot_res$result$result$cluster)),sep ="")))

cluster_assign <- bind_rows(cluster_assign_trns,cluster_assign_res)

assigned_dat_res <- dat_res %>%
  select(-outside_sites_total, -all_sites_total, 
         -exit_total,-entry_total, -all_total,
         -residence_mean_OUTSIDE, -residence_total_OUTSIDE,
         -mean_exit_woy, -mean_entry_woy, -time_btwn_mvmts_mean, 
        -residence_total_INSIDE) %>%
  rownames_to_column("TagAge") %>%
  left_join(cluster_assign_res)

assigned_dat_trns <- dat_trns %>%
  rownames_to_column("TagAge") %>%
  left_join(cluster_assign_trns)

assigned_dat <- bind_rows(assigned_dat_trns, assigned_dat_res) %>%
  filter(!is.na(cluster)) %>%
  filter(bootmean >= 0.6) %>%
  mutate(cluster = case_when(
    bootmean < 0.6 ~ "NA",
    TRUE ~ cluster
  ))

# library(plotly)
# p <- plot_ly(assigned_dat, x=~Comp.3, y=~Comp.2, 
#              z=~Comp.1, color=~cluster) %>%
#   add_markers(size=1.5)
# print(p)

# Interpretation ############################################################
## PERMANOVA ================================================================
datT2_assigned <- assigned_dat %>%
  select(-cluster, -TagAge) %>%
  decostand(method = "standardize") %>%
  bind_cols(select(assigned_dat, cluster))
  
## Pairwise PERMANOVA =======================================================

pairwise <- pairwise.adonis2(resp = vegdist(datT2_assigned[,1:(length(datT2_assigned)-1)],
                                            na.rm = TRUE), 
                             fact = factor(datT2_assigned$cluster))

pairwise$p.value
write.csv(pairwise$p.value, file.path("1. Data","Outputs","PairwisePERMANOVA.csv"))

## Kruskal-Wallis ==================================================================
# Run test on all the variables, all vars should be significant

kw_list <- list()
for(i in c("inside_sites_total", "outside_sites_total", 
           "all_sites_total", 
           "weekly_receivers_mean", "time_btwn_detects_total", "time_btwn_mvmts_mean", 
           "mean_distance", "q05_distance", "q25_distance", 
           "q50_distance", "q75_distance", "q99_distance",
           "mean_exit_woy", "mean_entry_woy",
           "exit_total", "all_total", "residence_total_INSIDE", 
           "residence_total_OUTSIDE", "residence_mean_INSIDE", 
           "residence_mean_OUTSIDE")){
  kw_list[[i]] <- kruskal.test(as.formula(paste(i,"~ cluster")), data = datT2_assigned)
}

resp_df <- data.frame()
for(i in names(kw_list)){
  response = i
  tmp <- kw_list[[i]]
  chi_sq = tmp$statistic[1]
  df = tmp$parameter[1]
  p_value = tmp$p.value[1]
  df <- data.frame(response, chi_sq,df,p_value)
  resp_df <- bind_rows(resp_df,df)
}
(resp_df <- resp_df %>% arrange(p_value) %>%
    mutate(p_value = round(p_value,6)) %>%
    remove_rownames())
write.csv(resp_df, file.path("1. Data","Outputs","KW_Response.csv"),
          row.names = FALSE)

dunn_list <- list()
for(i in c("inside_sites_total", "outside_sites_total", 
           "all_sites_total",  
           "weekly_receivers_mean", "time_btwn_detects_total", "time_btwn_mvmts_mean", 
           "mean_distance", "q05_distance", "q25_distance", 
           "q50_distance", "q75_distance", "q99_distance",
           "mean_exit_woy", "mean_entry_woy", "entry_total", 
           "exit_total", "all_total", "residence_total_INSIDE", 
           "residence_total_OUTSIDE", "residence_mean_INSIDE", 
           "residence_mean_OUTSIDE")){
  dunn_list[[i]] <-rstatix::dunn_test(as.formula(paste(i,"~ cluster")), data = datT2_assigned,
                             p.adjust.method = "bonferroni")
}
dunn_test <- bind_rows(dunn_list) %>%
  select(`.y.`,group1, group2, p.adj) %>%
  pivot_wider(names_from = group2, values_from = p.adj)
write.csv(dunn_test,file.path("1. Data","Outputs","Dunn_Test.csv"), 
          row.names = FALSE)

# View the means of each variable for each cluster
assigned_dat %>%
  group_by(cluster) %>%
  summarise(count = n(),
            across(where(~is.numeric(.x)), ~round(mean(.x, na.rm = TRUE),2))) %>%
  t()

### Feature Importance ====================================================
# Identify those features which have the most weight in the clustering process

#devtools::install_github("o1iv3r/FeatureImpCluster")

# FeatureImpCluster requires input from flexclust so convert
res_trns <- flexclust::as.kcca(kmeans_trns, data = datT2_trns, k = 3)
res_res <- flexclust::as.kcca(kmeans_res, data = datT2_res, k = 6)

# Get the Feature Importance
nr_seeds <- 1000
seeds_vec <- sample(1:10000, nr_seeds)

savedImp_trns <- data.frame(matrix(0,nr_seeds,dim(datT2_trns)[2]))
count <- 1
for (s in seeds_vec) {
  set.seed(s)
  res <- flexclust::kcca(datT2_trns, k=2)
  set.seed(s)
  FeatureImp_res <- FeatureImpCluster::FeatureImpCluster(res,as.data.table(datT2_trns),
                                                         sub = 1,biter = 1)
  savedImp_trns[count,] <- FeatureImp_res$featureImp[sort(names(FeatureImp_res$featureImp))]
  count <- count + 1
}
names(savedImp_trns) <- sort(names(FeatureImp_res$featureImp))

savedImp_res <- data.frame(matrix(0,nr_seeds,dim(datT2_res)[2]))
count <- 1
for (s in seeds_vec) {
  set.seed(s)
  res <- flexclust::kcca(datT2_res, k=4)
  set.seed(s)
  FeatureImp_res <- FeatureImpCluster::FeatureImpCluster(res,as.data.table(datT2_res),
                                                         sub = 1,biter = 1)
  savedImp_res[count,] <- FeatureImp_res$featureImp[sort(names(FeatureImp_res$featureImp))]
  count <- count + 1
}
names(savedImp_res) <- sort(names(FeatureImp_res$featureImp))

savedImp_trns <- savedImp_trns  %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "misclassification") %>%
  mutate(group = "transient")

savedImp_res <- savedImp_res  %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "misclassification") %>%
  mutate(group = "resident")

savedImp <- bind_rows(savedImp_trns, savedImp_res) %>%
  group_by(group, variable) %>%
  summarise(ymin = quantile(misclassification, 0.05),
         lower = quantile(misclassification, 0.25),
         median = quantile(misclassification, 0.5),
         upper = quantile(misclassification, 0.75),
         ymax = quantile(misclassification, 0.95)) %>%
  arrange(group, median)

trns_imp_plot <- ggplot(savedImp %>%
                          filter(group == "transient"))+
  geom_boxplot(aes(x = reorder(variable, median, descending = TRUE),
                   ymin = ymin, lower = lower, middle = median, upper = upper,
                   ymax = ymax, group = variable), stat = "identity")+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0),
    plot.margin = unit(c(0,0,0,0),"lines")
  )+
  lims(y = c(0,0.018))+
  coord_flip()+
  labs(x = "", y = "")

res_imp_plot <- ggplot(savedImp %>%
                         filter(group == "resident"))+
  geom_boxplot(aes(x = reorder(variable, median, descending = TRUE),
                   ymin = ymin, lower = lower, middle = median, upper = upper,
                   ymax = ymax, group = variable), stat = "identity")+
  scale_y_continuous(limits = c(0, 0.03))+
  theme_classic()+
  coord_flip()+
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0),
    plot.margin = unit(c(0,1,0,0),"lines")
  )+
  labs(x = "", y = "")

png(file.path("1. Data","Figures","FeatureImportance.png"),
    width = 6, height = 6,  units = "in", res = 600)
print(gridExtra::grid.arrange(trns_imp_plot, res_imp_plot, nrow = 1,
                              padding = unit(0,"line"), 
                              bottom = "Rate of Misclassification"))
dev.off()

## Examine the variables by cluster --------------------------------------

# Residency
residency <- ggplot(assigned_dat %>%
                      select(cluster, residence_total_INSIDE:residence_mean_OUTSIDE) %>%
                      pivot_longer(cols = residence_total_INSIDE:residence_mean_OUTSIDE,
                                   names_to = c("ignore","stat","location"),
                                   values_to = "value",
                                   names_sep = "_") %>%
                      mutate(Statistic = case_when(
                        grepl("mean", stat) ~ "Average",
                        grepl("total", stat) ~ "Total")))+
  geom_boxplot(aes(x = Statistic, y = value, color = location,
                   group = interaction(cluster, Statistic, location)),
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
                        select(cluster,TagAge, inside_sites_total:all_sites_total) %>%
                        pivot_longer(cols = inside_sites_total:all_sites_total,
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
         "group" = c("T1","T2","R1","R2","R3","R6"))

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
                        select(cluster,weekly_receivers_mean))+
  geom_boxplot(aes(y = weekly_receivers_mean, x = ""), linewidth = 0.75, 
               fill = "grey75")+
  facet_grid(cluster ~ ., scales = 'free_y')+
  coord_flip()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5))+
  labs(y = "Average Number of Site Visits per Week",
       x = "")
sites_plot2
ggsave(file.path("1. Data","Figures","AvgWeeklySites.png"), device = "png",
       width = 8, height = 6,        units = "in",dpi = 600)

# Number and Frequency of Movements
movements <- ggplot(assigned_dat %>%
                      select(cluster,entry_total:all_total) %>%
                      pivot_longer(cols = entry_total:all_total,
                                   names_to = c("Direction","stat"),
                                   names_sep = "_",
                                   values_to = "value"))+
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
  end = c(ymd("2018_03-14"),ymd("2018-05-14"),ymd("2018-07-31"),ymd("2018-10-31"),ymd("2018-12-31"))) %>%
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
                filter(exit_total >= 1 & mean_exit_woy > 0),
              aes(y = lubridate::ymd("2017-12-31")+(mean_exit_woy*7), x = 1), 
              fill = "grey50",
              trim = FALSE, bw = 7)+
  geom_text(data = assigned_dat%>%
              separate(TagAge,into = c("TagID","AgeBin"), sep = "-",
                       extra = "merge") %>%
              filter(mean_exit_woy > 1) %>%
              select(TagID, mean_exit_woy, cluster) %>%
              distinct() %>%
              group_by(cluster) %>%
              summarise(mean = mean(mean_exit_woy),
                        count = n()),
            aes(y = lubridate::ymd("2017-12-31")+(mean*7), x = 1.5, 
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
  coord_flip(xlim =c(0.5,2), expand = FALSE)+
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
                filter(mean_entry_woy > 1 & cluster %in% c("T1","T2","T3")),
              aes(y = lubridate::ymd("2017-12-31")+(mean_entry_woy*7), x = 1), 
              fill = "grey50",
              trim = FALSE, bw = 7)+
  geom_text(data = assigned_dat%>%
              separate(TagAge,into = c("TagID","AgeBin"), sep = "-",
                       extra = "merge") %>%
              filter(mean_entry_woy > 1 & cluster %in% c("T1","T2")) %>%
              select(TagID, mean_entry_woy, cluster) %>%
              distinct() %>%
              group_by(cluster) %>%
              summarise(mean = mean(mean_entry_woy),
                        count = n()),
            aes(y = lubridate::ymd("2017-12-31")+(mean*7), x = 1.5, 
                label = paste("N =",count, sep = " ")),
            inherit.aes = FALSE)+
  facet_grid(cluster ~ .,drop = FALSE)+
  coord_flip(xlim =c(0.5,2), expand = FALSE)+
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

mean_crossing_woy <- ggplot()+
  geom_ribbon(data = seasons_timing, aes(x = x, 
                                         ymin = start, 
                                         ymax = end,
                                         fill = season),
              alpha = 0.5)+
  geom_violin(data= assigned_dat%>%
                separate(TagAge,into = c("TagID","AgeBin"), sep = "-",
                         extra = "merge") %>%
                filter(ccf_cross > 0),
              aes(y = lubridate::ymd("2017-12-31")+(ccf_cross*7), x = 1), 
              fill = "grey50",
              trim = FALSE, bw = 7)+
  geom_text(data = assigned_dat%>%
              separate(TagAge,into = c("TagID","AgeBin"), sep = "-",
                       extra = "merge") %>%
              filter(ccf_cross > 1) %>%
              select(TagID, ccf_cross, cluster) %>%
              distinct() %>%
              group_by(cluster) %>%
              summarise(mean = mean(ccf_cross),
                        count = n()),
            aes(y = lubridate::ymd("2017-12-31")+(mean*7), x = 1.5, 
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
  coord_flip(xlim =c(0.5,2), expand = FALSE)+
  theme_classic()+
  scale_y_date(name = "Date", date_labels = "%b",
               breaks = "1 month")+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5))+
  labs(title = "Mean Week of crossing CCF", x = "")
mean_crossing_woy
ggsave(file.path("1. Data","Figures","CrossingWeek.png"), device = "png",
       width = 8, height = 6,        units = "in",dpi = 600)

### Time between Detects Plot ---------------------------------------------
time_detect_plot <- ggplot(assigned_dat %>%
                             select(cluster,TagAge,time_btwn_detects_total) %>%
                             pivot_longer(cols = time_btwn_detects_total,
                                          names_to = "Statistic",
                                          values_to = "value"))+
  geom_boxplot(aes(y = value, x = 1, color = Statistic), linewidth = 0.75,
               fill = "grey75")+
  geom_text(
    aes(label = paste(round(after_stat(count/2),0),"Detections"),
        x = 2,
        y = after_stat(17)),
    stat = "count",
    position = position_nudge(y = -2)
  )+
  facet_grid(cluster ~ ., scales = 'free_y')+
  coord_flip()+
  theme_classic()+
  theme(panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5))+
  scale_x_discrete(breaks = c("time_btwn_detects_total","time_btwn_detects_max"),
                   labels = c("Total","Max"))+
  scale_color_discrete(breaks = c("time_btwn_detects_total","time_btwn_detects_max"),
                       labels = c("Total","Max"))+
  lims(y = c(0,16))+
  labs(x = "", y = "Weeks between Detections")
time_detect_plot
ggsave(file.path("1. Data","Figures","TimeBtwnDetects.png"), device = "png",
       width = 8, height = 6, units = "in",dpi = 600)

### Time between Movements Plot --------------------------------------------
time_move_plot <- ggplot(assigned_dat %>%
                           select(cluster,TagAge,time_btwn_mvmts_mean) %>%
                           pivot_longer(cols = c("time_btwn_mvmts_mean"),
                                        names_to = "Statistic",
                                        values_to = "value") %>%
                           filter(value != 100) %>% filter(!grepl("R",cluster)))+
  geom_boxplot(aes(y = value, x = factor(Statistic), color = Statistic), linewidth = 0.75,
               fill = "grey75")+ 
  geom_text(
    aes(label = paste(round(after_stat(count/2),0),"Movements"),
        x = 1.5,
        y = after_stat(15)),
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
  labs(x = "", y = "Weeks between Transits")+
  lims(y = c(0,16))
time_move_plot
ggsave(file.path("1. Data","Figures","TimeBtwnMovements.png"), device = "png",
       width = 8, height = 6,        units = "in",dpi = 600)

### Age Plot ---------------------------------------------------------------
age_totals <- assigned_dat %>%
  separate(TagAge, into = c("tag","age", "extra"), sep ="-") %>%
  select(age) %>%
  mutate(age = case_when(
    age == "1" ~ "1-2",
    age == "3" ~ "3-5",
    age == "6+" ~ "6+",
    TRUE ~ NA
  )) %>%
  group_by(age) %>%
  tally(name = "age_total")

assigned_dat %>%
  separate(TagAge, into = c("tag","age", "extra"), sep ="-") %>%
  select(cluster, age) %>%
  mutate(age = case_when(
    age == "1" ~ "1-2",
    age == "3" ~ "3-5",
    age == "6+" ~ "6+",
    TRUE ~ NA
  )) %>%
  group_by(cluster, age) %>%
  tally() %>%
  left_join(age_totals) %>%
  filter(!is.na(age)) %>%
  mutate(perc_age = n/age_total) -> ages

age_plot <- ggplot(ages)+
  geom_bar(aes(x = "", y = perc_age, fill = cluster, group = cluster), 
           stat = "identity", position = "stack")+
  facet_grid(age ~ .)+
  coord_flip()+
  scale_fill_manual(name = "Cluster",
                      breaks = c("R3","R4",
                                 "T1","T2","T3"),
                      values= c("#bdd7e7","#6baed6",
                                "#feedde","#fdbe85","#fd8d3c"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5))+
  guides(guide_legend(nrow = 1, ncol = 7))+
  theme(legend.position = "top",
        legend.direction = "horizontal")
age_plot
ggsave(file.path("1. Data","Figures","Cluster_08_Age.png"), device = "png",
       width = 8, height = 6,        units = "in",dpi = 600)

### Water Year Type Plot ---------------------------------------------------
wyt_weights <- data.frame(year = 2013:2022) %>%
  mutate(wyt = get_water_year_type(year),
         duration = 52) %>%
  group_by(wyt) %>%
  summarise(duration = sum(duration)) %>%
  ungroup() %>%
  mutate(total_duration = sum(duration)) %>%
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
                               adj_count = `count`*(duration/total_duration)) %>%
                     distinct())+
  geom_bar(aes(x = "", y = adj_count, fill = cluster, group = cluster), 
           stat = "identity", position = "fill")+
  facet_grid(wyt ~ .)+
  coord_flip()+
  scale_fill_manual(name = "Cluster",
                    breaks = c("R3","R4",
                               "T1","T2","T3"),
                    values= c("#bdd7e7","#6baed6",
                              "#feedde","#fdbe85","#fd8d3c"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5))+
  theme(legend.position = "top")
wyt_plot

### Season Plot-----------------------------------------------------------------
seasons = data.frame(month = c(1:12),
                    season = c("Immigration","Immigration","Spawning","Spawning",
                               "Emmigration","Emmigration","Emmigration",
                               "Residence","Residence","Residence","Immigration",
                               "Immigration"))

season_plot <- ggplot(assigned_dat %>%
                        separate(TagAge,into = c("TagID","AgeBin"), sep = "-",
                                 extra = "merge") %>%
                        select(TagID, AgeBin, cluster) %>%
                        separate("AgeBin", into = c("AgeBin","extra"), sep = "-",
                                 extra = "merge") %>% 
                        separate("extra", into = c("Age2","year","month"), sep = "-") %>%
                        mutate(month = as.numeric(ifelse(is.na(month), year, month))) %>%
                        mutate(year = as.numeric(ifelse(Age2 %in% c(2013:2023), Age2, year))) %>%
                        select(-Age2) %>%
                        left_join(seasons) %>%
                        group_by(cluster,season) %>%
                        tally(name = "count") %>%
                        mutate(season = factor(season,
                                               levels = c("Immigration","Spawning",
                                                          "Emmigration",
                                                          "Residence"))))+
  geom_bar(aes(x = "", y = count, fill = cluster), 
           stat = "identity", position = "fill")+
  coord_flip()+
  facet_grid(season ~ .)+
  scale_fill_manual(name = "Cluster",
                    breaks = c("R3","R4",
                               "T1","T2","T3"),
                    values= c("#bdd7e7","#6baed6",
                              "#feedde","#fdbe85","#fd8d3c"))+
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

assigned_dat %>%
  separate(TagAge,into = c("TagID","AgeBin"), sep = "-",
           extra = "merge") %>%
  select(TagID, AgeBin, cluster) %>%
  separate("AgeBin", into = c("AgeBin","extra"), sep = "-",
           extra = "merge") %>% 
  separate("extra", into = c("Age2","year","month"), sep = "-") %>%
  mutate(month = as.numeric(ifelse(is.na(month), year, month))) %>%
  mutate(year = as.numeric(ifelse(Age2 %in% c(2013:2023), Age2, year))) %>%
  mutate(year = ifelse(month >= 11, year - 1, year)) %>%
  select(-Age2) %>%
  mutate(AgeBin = case_when(
    AgeBin == "1" ~ "1-2",
    AgeBin == "3" ~ "3-5",
    TRUE ~ "6+"
  )) %>%
  select(TagID, AgeBin, cluster, year, month) %>%
  distinct()-> tag_clusters

for(i in (tag_clusters %>% pull(var = "TagID") %>% unique())){
  weeks <- TagBKM_Bin %>%
    mutate(TagID = as.character(TagID)) %>%
    filter(TagID == i) %>%
    pull(Week)
  
  tmp <- TagBKM_Bin %>%
    mutate(TagID = as.character(TagID)) %>%
    select(TagID,Week,SiteCode,last_detect) %>%
    filter(TagID == i) %>%
    arrange(Week, SiteCode) %>%
    full_join(data.frame(TagID = as.character(i),
                         Week = c(min(weeks):max(weeks)))) %>%
    arrange(Week) %>%
    mutate(month = month(studyweek_startdate(Week)),
           year = year(studyweek_startdate(Week))) %>%
    left_join(tag_clusters) %>%
    mutate(cluster = ifelse(is.na(cluster),"NA",cluster))
  
  cl_month <- tmp %>%
    select(year, month, cluster) %>%
    distinct() %>%
    crossing(y = c(0,16))
  
  if(nrow(tmp)>1){
    max_site <- length(unique(tmp$SiteCode))+1
    ggplot()+
      # Season Background
      geom_ribbon(data = cl_month %>%
                    ungroup() %>%
                    mutate(y = ifelse(y == 16, max_site, y)) %>%
                    group_by(year,month,cluster,y) %>%
                    mutate(xmin = ymd(paste(year, month, 1, sep = "-")),
                           xmax = ymd(paste(year, month, 
                                            days_in_month(month), 
                                            sep = "-"))),
                  aes(y = y,
                  xmin = xmin,
                  xmax = xmax,
                  fill = cluster,
                  group = interaction(year, month, cluster)),
                  alpha = 0.5)+
      # Detection History
      geom_point(data = tmp %>% filter(!is.na(SiteCode)),
                 aes(x = studyweek_startdate(Week), y = SiteCode, group = TagID))+
      geom_step(data = tmp  %>% filter(!is.na(SiteCode)),
                aes(x = studyweek_startdate(Week), y = SiteCode, group = TagID,
                    #linetype = factor(SiteCode)
                    ))+
      # Tag Failure Line
      geom_vline(data = tmp %>% filter(SiteCode == "Tag Failure"), 
                 aes(xintercept = studyweek_startdate(Week)), 
                 color = "red", linewidth = 1.5)+
      geom_text(data = tmp %>% filter(SiteCode == "Tag Failure"), 
                aes(x = studyweek_startdate(Week-3),
                    y = "RGD1", label = "End Detection History"),
                angle = 90, color = "red")+
      # CCF Inside/Outside Separating Line
      geom_hline(yintercept = as.numeric(factor("RGD1",
                                                levels = levels(tmp$SiteCode))) + 0.5,
                                                color = "blue", 
                 linewidth = 1.5, linetype = "dashed")+
      # Age
      geom_text(data = tmp %>%
                  ungroup() %>%
                  select(SiteCode, AgeBin, year, month, cluster) %>%
                  mutate(SiteCode = max(SiteCode)) %>%
                  filter(!is.na(AgeBin)) %>%
                  group_by(SiteCode, AgeBin) %>%
                  summarise(year = min(year,na.rm = TRUE),
                            month = min(month, na.rm = TRUE)),
                aes(x = ymd(paste(year, month, 15, sep = "-")),
                    y = "Tag Failure", label = paste("\nAge:", AgeBin),
                    group = AgeBin
                ), check_overlap = TRUE, hjust = 0)+
      # Plot Setup
      scale_x_date(name = "Date",
                   breaks = scales::breaks_pretty(n = 12),
                   date_labels = "%b-%Y")+
      scale_y_discrete(breaks = c("Release","IC3","IC2","IC1","RGD1","RGU1",
                                  "ORS1","WC1","ORN1","WC2","ORN2","GL1","CVP1",
                                  "WC3","ORS3","CLRS"))+
      scale_fill_manual(name = "Cluster",
                        breaks = c("R1","R2","R3","R4",
                                   "T1","T2","T3","NA"),
                        values= c("gray50","gray50","#bdd7e7","#6baed6",
                                  "#feedde","#fdbe85","#fd8d3c","gray50"))+
      scale_linetype_manual(breaks = c("Tag Failure","Release","IC3","IC2","IC1","RGD1","RGU1",
                                       "ORS1","WC1","ORN1","WC2","ORN2","GL1","CVP1",
                                       "WC3","ORS3","CLRS"),
                            values = c("dashed","solid","solid","solid","solid",
                                       "solid","solid","solid","solid","solid",
                                       "solid","solid","solid","solid","solid",
                                       "solid","solid","solid")
                            )+
      theme_classic()+
      theme(panel.spacing = unit(0, "lines"),
            axis.text.x = element_text(angle = -45, hjust = 0),
            axis.title.y = element_blank())+
      labs(subtitle = paste("Tag ID:",tmp$TagID[1]))+
      coord_cartesian(xlim = c(studyweek_startdate(min(tmp$Week)-1),
                               studyweek_startdate(max(tmp$Week)+1)),
                      #ylim = c(1,16)
                      )
    ggsave(file.path("1. Data","Figures","DetectionPlots",paste(tmp$TagID[1],"DetectionHistory_Plain.png",sep = "-")),
           device = "png", width = 10, height = 5, dpi = 600)
  } else {
    print("Tag not in cluster")
  }
}

#Explanatory Variables
assigned_dat %>%
  separate(TagAge,into = c("TagID","AgeBin"), sep = "-",
           extra = "merge") %>%
  select(TagID, AgeBin, cluster, wyt, gate_perc_open) %>%
  separate("AgeBin", into = c("AgeBin","extra"), sep = "-",
           extra = "merge") %>% 
  separate("extra", into = c("Age2","year","month"), sep = "-") %>%
  mutate(month = as.numeric(ifelse(is.na(month), year, month))) %>%
  mutate(year = as.numeric(ifelse(Age2 %in% c(2013:2023), Age2, year))) %>%
  mutate(AgeBin = factor(case_when(
    AgeBin == "1" ~ "1-2",
    AgeBin == "3" ~ "3-5",
    TRUE ~ "6+"), levels = c("1-2","3-5","6+")),
    wyt = factor(wyt, levels = c(1:5),
                 labels = c("Critical","Dry","Below Normal","Above Normal",
                            "Wet"))) %>%
  select(-Age2) %>%
  left_join(seasons) %>%
  mutate(season = factor(season, levels = c("Immigration","Spawning","Emmigration",
                                            "Residence")),
         cluster2 = factor(cluster, levels = c("R3","R4","T3","T2","T1"))
         ) -> dat

age_dunn <- FSA::dunnTest(as.numeric(AgeBin) ~ cluster, data = dat, method = "bonferroni")
wyt_dunn <- FSA::dunnTest(as.numeric(wyt) ~ cluster, data = dat)
season_dunn <- FSA::dunnTest(as.numeric(season) ~ cluster, data = dat)
gate_dunn <- FSA::dunnTest(gate_perc_open ~ cluster, data = dat)

age_cld <- rcompanion::cldList(comparison = age_dunn$res$Comparison,
               p.value    = age_dunn$res$P.adj,
               threshold  = 0.05)[1:2]

wyt_cld <- rcompanion::cldList(comparison = wyt_dunn$res$Comparison,
                               p.value    = wyt_dunn$res$P.adj,
                               threshold  = 0.05)[1:2]

season_cld <- rcompanion::cldList(comparison = season_dunn$res$Comparison,
                               p.value    = season_dunn$res$P.adj,
                               threshold  = 0.05)[1:2]

gate_cld <- rcompanion::cldList(comparison = gate_dunn$res$Comparison,
                                  p.value    = gate_dunn$res$P.adj,
                                  threshold  = 0.05)[1:2]

names(age_cld) <- c("cluster","Letter")
names(wyt_cld) <- c("cluster","Letter")
names(season_cld) <- c("cluster","Letter")
names(gate_cld) <- c("cluster","Letter")

ggplot(dat) +
  geom_boxplot(aes(x = cluster, y = gate_perc_open))+
  scale_color_brewer(palette = "Dark2")+
  theme_bw() +
  geom_text(data = gate_cld, aes(label = Letter, 
                                 y = max(dat$gate_perc_open), 
                                 x = cluster,
                                 color = Letter), 
            vjust = -0.5,
            hjust= 0.5,
            size=3.5,
            check_overlap = F)+
  theme(legend.position = "none")+
  guides(color = "none")+
  labs(y = "Percent of Time Gate was Open",
       x = "Cluster")

ggsave("GatePercentage_Plot.png", device = "png", width = 8, height = 6)

ggplot(dat) +
  geom_bar(aes(x = cluster, fill = AgeBin, group = AgeBin),
           position = "fill")+
  scale_fill_discrete(name = "Age")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw() +
  geom_text(data = age_cld, aes(label = Letter, 
                                 y = 1, 
                                 x = cluster,
                                 color = Letter), 
            vjust = -0.5,
            hjust= 0.5,
            size=3.5,
            check_overlap = F)+
  guides(color = "none")+
  labs(y = "Percent of Behavior States at Age", x = "Cluster")
ggsave("Age_Plot.png", device = "png", width = 8, height = 6)

ggplot(dat) +
  geom_bar(aes(x = cluster, fill = season, group = season),
           position = "fill")+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_discrete(name = "season")+
  theme_bw() +
  geom_text(data = season_cld, aes(label = Letter, 
                                y = 1, 
                                x = cluster,
                                color = Letter), 
            vjust = -0.5,
            hjust= 0.5,
            size=3.5,
            check_overlap = F)+
  guides(color = "none")+
  labs(x = "Cluster",
       y = "Percent of Behavior States in each Season")
ggsave("Season_Plot.png", device = "png", width = 8, height = 6)

ggplot(dat) +
  geom_bar(aes(x = cluster, fill = wyt, group = wyt),
           position = "fill")+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_discrete(name = "Water Year Type")+
  theme_bw() +
  geom_text(data = wyt_cld, aes(label = Letter, 
                                   y = 1, 
                                   x = cluster,
                                   color = Letter), 
            vjust = -0.5,
            hjust= 0.5,
            size=3.5,
            check_overlap = F)+
  guides(color = "none")+
  labs(x = "Cluster",
       y = "Percent of Behavior States in each Water Year Type")
ggsave("WYT_Plot.png", device = "png", width = 8, height = 6)
