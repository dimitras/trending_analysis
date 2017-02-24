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




