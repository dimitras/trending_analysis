####################################################################################
####################################################################################
# expression ratios histograms

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# radiation.data <- read.table("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/FINAL_master_list_of_gene_counts_MIN.sense.radiation_renamed.sorted_ratios.txt", header = TRUE)
radiation.data <- read.table("~/Dropbox/workspace/radiation/results/dose_response_curves/FINAL_master_list_of_gene_counts_MIN.antisense.radiation_renamed.sorted_ratios_updated.txt", header = TRUE)

radiation.data %>%
  gather(Sample_id,Ratio,Proton5.S2.Proton0.S1:Gamma200.S20.Gamma0.S11) %>%
  mutate(Dose=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\2",Sample_id),
         Dose = as.numeric(Dose),
         Radiation_type=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\1",Sample_id)) %>%
  
  # mutate(mean5 = mean(subset(radiation.data, select = c(Proton5.S2.Proton0.S1, Gamma5.S12.Gamma0.S11)), na.rm=TRUE),
  #        mean5 = as.numeric(mean5),
  #        median5 = median(subset(radiation.data, select = c(Proton5.S2.Proton0.S1, Gamma5.S12.Gamma0.S11)), na.rm=TRUE),
  #        median5 = as.numeric(median5)
  # sd?
  # quantiles?
  # ) %>%
  
  # ggplot(aes(x=Ratio)) + geom_histogram(binwidth=0.5, aes(fill=Radiation_type), position = "dodge") + #for 20 bins
  # ggplot(aes(x=Ratio)) + geom_histogram(binwidth=0.03, aes(fill=Radiation_type), position = "dodge") + #for 300 bins
  ggplot(aes(x=Ratio)) + geom_histogram(binwidth=0.001, aes(fill=Radiation_type), position = "dodge") + #for 1000 bins
    # geom_vline(aes(xintercept=mean(Ratio)), colour="blue") +
    # geom_vline(aes(xintercept=median(Ratio)), colour="green") +
    geom_freqpoly(aes(color=Radiation_type)) + 
    scale_x_log10() +
    facet_grid(.~Dose)
    
# ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/radiation_ratios_comparison_1000bins_wlog.tiff", height = 15, width = 100, units ="cm")
ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/antisense_radiation_ratios_comparison_1000bins_wlog.tiff", height = 15, width = 100, units ="cm")



####################################################################################
####################################################################################
# expression ratios histograms for genes with FDR<0.1

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)


# radiation.data <- read.table("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/all_sense_0.1_ratios.txt", sep="\t", header = TRUE)
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
ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/all_PCR_0.1_ratios_comparison_1000bins_wlog.tiff", height = 60, width = 60, units ="cm")



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
radiation.data <- read.table("~/Dropbox/workspace/radiation/trend_data/RT-PCR_gene_expressions.txt", header = TRUE)

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


