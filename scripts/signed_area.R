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
signed.data <-  read.table("~/Documents/dimitra/Workspace/RNA-Seq/radiation/results/signed_areas/RNA-Seq_signed_areas_reformatted.txt", sep="\t", header = TRUE)
sid = signed.data %>% 
  arrange(desc(signed_area))
  
pval = t.test(sid$signed_area)$p.value
sid.mean = mean(sid$signed_area)

sid %>%
  ggplot(aes(x=signed_area) ) +
    geom_histogram(binwidth=0.01, position = "identity") +
    theme(text = element_text(size=40))
    # annotate("text",x=sid.mean,y=1.1,label= pval, size=12) 


rad.data = read.table("~/Documents/dimitra/Workspace/RNA-Seq/radiation/results/signed_areas/signed_areas_up_and_downs.txt", sep="\t", header = TRUE)
rad.data %>%
  gather(Domination,GenesCount,Gamma:Proton) %>%
  mutate(Expression=gsub("^(Proton|Gamma)","\\1",Domination)) %>%
  
  ggplot(aes(x=Expression, y=GenesCount, fill=Expression) ) +
  scale_y_log10() +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(text = element_text(size=40)) 

ggsave("~/Dropbox/workspace/radiation/results/signed_areas/signed_areas_number_of_dominants.tiff", height = 60, width = 60, units ="cm")




# heatmap for signed area
library(ggplot2)
library(data.table)
library(scales)

# signed.data = fread("~/Documents/dimitra/Workspace/RNA-Seq/radiation/results/signed_areas/signed_areas_for_19_common_genes.txt")
signed.data = fread("~/Documents/dimitra/Workspace/RNA-Seq/radiation/results/signed_areas/RNA-Seq_signed_areas_reformatted.txt")

signed.data = signed.data[,.(id, signed_area, signed_area_cp = signed_area)]

signed.data %>%
  arrange(signed_area) %>% 
  mutate(id=factor(id,levels=id)) %>%
  
  ggplot(aes(x = "", y = id, fill = signed_area)) +
  geom_tile() +
  # scale_fill_gradientn(colours = c("turquoise1", "white", "indianred1"), values = rescale(c(1300,0,-1300))) +
  scale_fill_gradientn(colours = c("turquoise1", "white", "indianred1"), values = rescale(c(42000,0,-96000))) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # axis.text.y = element_text(size = 20),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        panel.grid = element_blank())

# ggsave("~/Documents/dimitra/Workspace/RNA-Seq/radiation/results/signed_areas/signed_areas_for_19_common_genes_heatmap.tiff", height = 35, width = 15, units ="cm")
ggsave("~/Documents/dimitra/Workspace/RNA-Seq/radiation/results/signed_areas/RNA-Seq_signed_areas_heatmap.tiff", height = 45, width = 15, units ="cm")
