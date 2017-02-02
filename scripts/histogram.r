####################################################################################
####################################################################################
# expression ratios histograms

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

radiation.data <- read.table("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/FINAL_master_list_of_gene_counts_MIN.sense.radiation_renamed.sorted_ratios.txt", header = TRUE)
# radiation.data <- read.table("~/Dropbox/workspace/radiation/results/dose_response_curves/FINAL_master_list_of_gene_counts_MIN.antisense.radiation_renamed.sorted_ratios_updated.txt", header = TRUE)

radiation.data %>%
  gather(Sample_id,Ratio,Proton5.S2.Proton0.S1:Gamma200.S20.Gamma0.S11) %>%
  mutate(Dose=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\2",Sample_id),
         Dose = as.numeric(Dose),
         Radiation_type=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\1",Sample_id)) %>%
  
  # ggplot(aes(x=Ratio)) + geom_histogram(binwidth=0.5, aes(fill=Radiation_type), position = "dodge") + #for 20 bins
  # ggplot(aes(x=Ratio)) + geom_histogram(binwidth=0.03, aes(fill=Radiation_type), position = "dodge") + #for 300 bins
  ggplot(aes(x=Ratio)) + geom_histogram(binwidth=0.001, aes(fill=Radiation_type), position = "dodge") + #for 1000 bins
    # geom_vline(aes(xintercept=mean(Ratio)), colour="blue") +
    # geom_vline(aes(xintercept=median(Ratio)), colour="green") +
    geom_freqpoly(aes(color=Radiation_type, size=10)) + 
    scale_x_log10() +
    # facet_grid(.~Dose)
    facet_wrap(~Dose, ncol = 3) +
    theme(text = element_text(size=40)) + 
    labs(y="Gene expression", x="Expression ratio")

ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/all_sense_expression_ratios_comparison_1000bins_wlog_3x3.tiff", height = 80, width = 80, units ="cm")
# ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/antisense_radiation_ratios_comparison_1000bins_wlog.tiff", height = 15, width = 100, units ="cm")

# CDF (cumulative distribution)
radiation.data %>%
  gather(Sample_id,Ratio,Proton5.S2.Proton0.S1:Gamma200.S20.Gamma0.S11) %>%
  mutate(Dose=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\2",Sample_id),
         Dose = as.numeric(Dose),
         Radiation_type=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\1",Sample_id)) %>%
  
  ggplot(aes(x=Ratio, color=Radiation_type)) + geom_step(stat="ecdf", aes(size=8)) +
  scale_x_log10() +
  facet_wrap(~Dose, ncol = 3) +
  theme(text = element_text(size=40)) + 
  labs(y=" ", x="Expression ratio")
  
ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/all_sense_cdf_expression_ratios_comparison_wlog_3x3.tiff", height = 80, width = 80, units ="cm")


####################################################################################
####################################################################################
# expression ratios histograms for genes with FDR<0.1

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)


radiation.data <- read.table("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/all_sense_0.1_ratios.txt", sep="\t", header = TRUE)
radiation.data <- read.table("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/all_PCR_0.1_ratios.txt", sep="\t", header = TRUE)

radiation.data %>%
gather(Sample_id,Ratio,Proton5.S2.baseline:Gamma200.S20.baseline) %>%
mutate(Dose=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\2",Sample_id),
      Dose = as.numeric(Dose),
      Radiation_type=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\1",Sample_id)) %>%

ggplot(aes(x=Ratio)) + 
  geom_histogram(binwidth=0.001, aes(fill=Radiation_type), position = "dodge") +
  geom_freqpoly(aes(color=Radiation_type, size=10)) + 
  scale_x_log10() +
  facet_wrap(~Dose, ncol = 3) +
  theme(text = element_text(size=40)) + 
  labs(y="Gene expression", x="Expression ratio")

# ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/all_sense_0.1_ratios_comparison_1000bins_wlog.tiff", height = 60, width = 60, units ="cm")
ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/all_PCR_0.1_ratios_comparison_1000bins_wlog_3x3.tiff", height = 60, width = 60, units ="cm")


# CDF (cumulative distribution)
radiation.data %>%
  gather(Sample_id,Ratio,Proton5.S2.baseline:Gamma200.S20.baseline) %>%
  mutate(Dose=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\2",Sample_id),
         Dose = as.numeric(Dose),
         Radiation_type=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\1",Sample_id)) %>%
  
  ggplot(aes(x=Ratio, color=Radiation_type, size=6)) + geom_step(stat="ecdf") +
  scale_x_log10() +
  facet_wrap(~Dose, ncol = 3) +
  theme(text = element_text(size=40)) + 
  labs(y=" ", x="Expression ratio")

ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/all_sense_0.1_cdf_expression_ratios_comparison_wlog_3x3.tiff", height = 80, width = 80, units ="cm")
ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/all_PCR_0.1_cdf_expression_ratios_comparison_wlog_3x3.tiff", height = 80, width = 80, units ="cm")


# Dose - Quartile expression ratio (treated/control)
radiation.data <- read.table("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/all_sense_0.1_ratios.txt", sep="\t", header = TRUE)

radiation.data %>%
  gather(Sample_id,Ratio,Proton5.S2.baseline:Gamma200.S20.baseline) %>%
  mutate(Dose=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\2",Sample_id),
         Dose = as.numeric(Dose),
         Radiation_type=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\1",Sample_id)) %>%

  ggplot(aes(x=Dose, y=Ratio, color=Radiation_type, size = 5)) +
  theme(text = element_text(size=40)) + 
  # labs(y="Median of expression radiated/control") + 
  # stat_summary(fun.y = median, geom = "line", position = "identity")
  # labs(y="Q1 of expression radiated/control") +
  # stat_summary(fun.y = function(z) { quantile(z,0.25) }, geom = "line", position = "identity")
  # labs(y="Q3 of expression radiated/control") +
  # stat_summary(fun.y = function(z) { quantile(z,0.75) }, geom = "line", position = "identity")

  #all in one
  labs(y="Quartiles of expression radiated/control") +
  stat_summary(fun.ymin = function(z) { quantile(z,0.25) },
               fun.ymax = function(z) { quantile(z,0.75) },
               fun.y = median,
               position = "identity")

  ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/Q1_expression_ratios.tiff", height = 80, width = 80, units ="cm")
  ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/Q3_expression_ratios.tiff", height = 80, width = 80, units ="cm")
  ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/Median_expression_ratios.tiff", height = 80, width = 80, units ="cm")
  
  ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/summary_stats_expression_ratios.tiff", height = 80, width = 80, units ="cm")



####################################################################################
####################################################################################
# expression for genes with FDR<0.1 per dose (1..10), per radiation (P, G), per experiment (RNA-Seq, RT-PCR)

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)


radiation.data <- read.table("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/common_genes_across_experiments_0.1.txt", sep="\t", header = TRUE)
radiation.data <- read.table("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/union_genes_across_experiments_0.1.txt", sep="\t", header = TRUE)

radiation.data %>%
  gather(Sample_id,Expression,RNA.Gamma0.S11:PCR.Proton200.S10) %>%
  mutate(Dose=gsub("^(RNA|PCR).(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\3",Sample_id),
         Dose = factor(Dose, levels=seq(from=0,to=200,by=5), labels=seq(from=0,to=200,by=5)),
         # Dose = as.numeric(Dose),
         Radiation_type=gsub("^(RNA|PCR).(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\2",Sample_id),
         Experiment=gsub("^(RNA|PCR).(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\1",Sample_id),
         Experiment_Radiation=gsub("^(RNA.Gamma|RNA.Proton|PCR.Gamma|PCR.Proton)([[:digit:]]+)\\.S[[:digit:]]+.*","\\1",Sample_id)) %>%

  ggplot(aes(x=Dose, y=Expression, fill=Radiation_type) ) + 
  geom_boxplot(position = "dodge", na.rm = TRUE) +
  facet_grid(.~Experiment) +
  theme(text = element_text(size=40)) 
  # scale_x_discrete(drop=FALSE, breaks=c(0,50,100,150,200))
  
  ## or plot 2x2 facets
  # ggplot(aes(x=Dose, y=Expression, group=Dose, fill=Radiation_type)) + 
  # geom_boxplot(position = "dodge") +
  # facet_grid(Experiment~Radiation_type)

ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/common_genes_expression_across_experiments_0.1.tiff", height = 60, width = 60, units ="cm")
ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/union_genes_expression_across_experiments_0.1.tiff", height = 60, width = 60, units ="cm")




####################################################################################
####################################################################################
# signed area comparison per experiment (RNA-Seq, RT-PCR)

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(plotly)


radiation.data <- read.table("~/Dropbox/workspace/radiation/results/signed_areas/signed_areas_summary.txt", sep="\t", header = TRUE)

radiation.data %>%
  gather(Domination,GenesCount,Gamma:No.change) %>%
  mutate(Expression=gsub("^(Proton|Gamma|No-change)","\\1",Domination)) %>%
  
  ggplot(aes(x=Expression, y=GenesCount, fill=Expression) ) +
  scale_y_log10() +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_grid(.~Experiment) +
  theme(text = element_text(size=40)) 

ggsave("~/Dropbox/workspace/radiation/results/signed_areas/signed_areas_comparison.tiff", height = 60, width = 60, units ="cm")



# heatmap for signed area
library(scales)
signed.data <-  read.table("~/Dropbox/workspace/radiation/results/signed_areas/RNA-Seq_signed_areas_reformatted.txt", sep="\t", header = TRUE)
# signed.data <-  read.table("~/Dropbox/workspace/radiation/results/signed_areas/RT-PCR_signed_areas_reformatted.txt", sep="\t", header = TRUE)

signed.data %>% 
  arrange(desc(signed_area)) %>%
  # mutate(tissue=gsub("^+_(+)","\\1",id)) %>%
  
  ggplot(aes(x = "", y = id, fill = signed_area)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("darkblue",  "blue", "purple", "green", "yellow", "orange", "red", "darkred"), values = rescale(c(41849,20000,1,-20000,-40000,-60000,-80000,-96010))) + # RNA
  # scale_fill_gradientn(colours = c("darkblue",  "blue", "purple", "green", "yellow", "orange", "red", "darkred"), values = rescale(c(211435,-300000,-800000, -1300000, -1800000, -2300000, -2800000, -3865063))) + # PCR
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        panel.grid = element_blank(),
        strip.background = element_blank())

# ggsave("~/Dropbox/workspace/radiation/results/signed_areas/RT-PCR_signed_areas_heatmap.tiff", height = 25, width = 25, units ="cm")
ggsave("~/Dropbox/workspace/radiation/results/signed_areas/RNA-Seq_signed_areas_heatmap.tiff", height = 25, width = 25, units ="cm")



####################################################################################
####################################################################################
# dose response curve


library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(drc)
library(reshape2)

# radiation.data <- read.table("~/Dropbox/workspace/radiation/expression_data/FINAL_master_list_of_gene_counts_MIN.sense.radiation_renamed.sorted.txt", header = TRUE)
radiation.data <- read.table("~/Dropbox/workspace/radiation/trend_data/RT-PCR_gene_expressions_edited.txt", header = TRUE)

### Dose response box plots  

radiation.data %>% 
gather(Sample_id,Expression,Proton0.S1:Gamma200.S20) %>% 
mutate(Dose=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+","\\2",Sample_id),
       Dose = as.numeric(Dose),
       Radiation_type=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+","\\1",Sample_id)) %>% 
   
ggplot(aes(x=Dose, y=Expression, fill=Radiation_type, group=Dose)) + 
  scale_y_log10() +
  geom_boxplot(position = "dodge")
  # geom_smooth(method = "drm", fct = LL.4(), se = TRUE)

# ggsave("~/Dropbox/workspace/radiation/results/dose_response_curves/sense_dose_expression_boxplots.tiff", height = 40, width = 40, units ="cm")
ggsave("~/Dropbox/workspace/radiation/results/dose_response_curves/PCR_dose_expression_boxplots.tiff", height = 40, width = 40, units ="cm")


### Dose response curves

radiation.data %>% 
gather(Sample_id,Expression,Proton0.S1:Gamma200.S20) %>% 
  mutate(Dose=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+","\\2",Sample_id),
         Dose = as.numeric(Dose),
         Radiation_type=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+","\\1",Sample_id)) %>% 
  
  ggplot(aes(x=Dose, y=Expression)) + 
  geom_line(aes(color=Radiation_type))

# ggsave("~/Dropbox/workspace/radiation/results/dose_response_curves/sense_dose_response_curves.tiff", height = 40, width = 60, units ="cm")
ggsave("~/Dropbox/workspace/radiation/results/dose_response_curves/PCR_dose_response_curves.tiff", height = 40, width = 40, units ="cm")


#### Plot summary statistics

radiation.data %>% 
  gather(Sample_id,Expression,Proton0.S1:Gamma200.S20) %>% 
  mutate(Dose=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+","\\2",Sample_id),
         Dose = as.numeric(Dose),
         Radiation_type=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+","\\1",Sample_id)) %>% 
  
  ggplot(aes(x=Dose, y=Expression, color=Radiation_type)) + 
  stat_summary(geom = "line", position = "identity", fun.y = "mean") # geom = "point"

# ggsave("~/Dropbox/workspace/radiation/results/dose_response_curves/sense_mean_expressions.tiff", height = 40, width = 40, units ="cm")
ggsave("~/Dropbox/workspace/radiation/results/dose_response_curves/PCR_mean_expressions.tiff", height = 40, width = 40, units ="cm")




##### Pathway heatmaps

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(scales)

ipadata <- read.csv("~/Dropbox/workspace/radiation/IPA/non-curated_pathways_fdr0.1.csv")
ipadata <- read.csv("~/Dropbox/workspace/radiation/IPA/curated_pathways_fdr0.1.csv")

ipadata.for_hclust = ipadata %>% column_to_rownames("Pathway") %>% select(Proton:Gamma) %>% as.matrix
ipadata.hclust_results = hclust(as.dist((1 - cor(t(ipadata.for_hclust)))))

ipadata %>% 
  arrange(Proton) %>%
  mutate(Pathway=factor(Pathway,levels=as.character(Pathway))) %>% 
  gather(Radiation,logPvalue,Proton:Gamma) %>%
  mutate(Radiation=factor(Radiation,levels=c("Proton","Gamma")),
         Pathway = factor(Pathway, levels=ipadata.hclust_results$labels[ipadata.hclust_results$order])) %>%
  ggplot(aes(x = Radiation, y = Pathway, fill=logPvalue)) +
  geom_tile() +
  # scale_fill_gradient(low = "lightblue", high="darkblue") +
  scale_fill_gradientn(colours = c("darkred", "red", "orange", "yellow", "green", "blue", "darkblue"), values = rescale(c(6,5,4,3,2,1,0))) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 13),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 15),
      panel.grid = element_blank())

ggsave("~/Dropbox/workspace/radiation/results/pathways/non-curated_pathways_fdr0.1_corclust.tiff", height = 40, width = 40, units ="cm")
ggsave("~/Dropbox/workspace/radiation/results/pathways/curated_pathways_fdr0.1_corclust.tiff", height = 40, width = 40, units ="cm")


