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


#
signed.data <-  read.table("~/Dropbox/workspace/radiation/results/signed_areas/RNA-Seq_signed_areas_reformatted.txt", sep="\t", header = TRUE)
sid = signed.data %>% 
  arrange(desc(signed_area))
pval = t.test(sid$signed_area)$p.value
sid.mean = mean(sid$signed_area)

ggplot(aes(x=signed_area) ) +
  geom_histogram(binwidth=0.01, position = "identity") +
  theme(text = element_text(size=40)) +
  annotate("text",x=sid.mean,y=1.1,label= pvals[1], size=12) +


radiation.data <- read.table("~/Dropbox/workspace/radiation/results/signed_areas/signed_areas_up_and_downs.txt", sep="\t", header = TRUE)
radiation.data %>%
  gather(Domination,GenesCount,Gamma:Proton) %>%
  mutate(Expression=gsub("^(Proton|Gamma)","\\1",Domination)) %>%
  
  ggplot(aes(x=Expression, y=GenesCount, fill=Expression) ) +
  scale_y_log10() +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(text = element_text(size=40)) 

ggsave("~/Dropbox/workspace/radiation/results/signed_areas/signed_areas_number_of_dominants.tiff", height = 60, width = 60, units ="cm")




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



