###########################
###########################
##### Pathway heatmaps ####

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


# Simpler heatmap
library(scales)
# ipadata <- read.csv("~/Dropbox/workspace/radiation/IPA/curated_pathways_fdr0.1.csv")
# ipadata <- read.csv("~/Dropbox/workspace/radiation/IPA/curated300genes_at.least.1.larger.than.neglog.05.csv")
ipadata <- read.csv("~/Documents/dimitra/Workspace/RNA-Seq/radiation/IPA/curated300genes_at.least.1.larger.than.neglog.01.csv")

# ipadata.for_hclust = ipadata %>% column_to_rownames("Pathway") %>% select(Gamma:Proton) %>% as.matrix
# ipadata.hclust_results = hclust(as.dist((1 - cor(t(ipadata.for_hclust)))))

ipadata %>%
  # arrange(Gamma) %>%
  arrange(desc(Sorted)) %>%
  mutate(Pathway=factor(Pathway,levels=as.character(Pathway))) %>% 
  gather(Radiation,logPvalue,Proton:Gamma) %>%
  mutate(Radiation=factor(Radiation,levels=c("Gamma", "Proton"))) %>%
         # Pathway = factor(Pathway, levels=ipadata.hclust_results$labels[ipadata.hclust_results$order])) %>%
  mutate(logPvalue=replace(logPvalue, logPvalue<2, NA)) %>%
  
  ggplot(aes(x = Radiation, y = Pathway, fill=logPvalue)) +
  geom_tile() +
  # scale_fill_gradient(low = "white", high="darkgreen", name="-logPvalue") +
  # scale_fill_gradient(low = "white", high="darkgreen", name = "-logPvalue", labels = c("<2", "3", "4", "5"), breaks = c(2,3,4,5)) +
  # scale_fill_gradientn(name = "-logPvalue", colours = c("darkgreen" ,"green", "lightgreen", "white"), values = rescale(c(4,3,2)), labels = c("<2", "3", "4","5"), breaks = c(2,3,4,5)) +
  scale_fill_gradientn(name = "-logPvalue", colours = c("darkgreen" ,"green", "lightgreen", "white"), values = rescale(c(4,3,2)), labels = c("2", "3", "4","5"), breaks = c(2.05,3,4,5), na.value = "grey80") +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        panel.grid = element_blank())

# ggsave("~/Dropbox/workspace/radiation/results/pathways/curated300genes_at.least.1.larger.than.neglog.05.tiff", height = 35, width = 30, units ="cm")
# ggsave("~/Dropbox/workspace/radiation/results/pathways/curated300genes_at.least.1.larger.than.neglog.05.eps", height = 35, width = 30, units ="cm")
# ggsave("~/Dropbox/workspace/radiation/results/pathways/curated300genes_at.least.1.larger.than.neglog.01.tiff", height = 35, width = 30, units ="cm")
# ggsave("~/Dropbox/workspace/radiation/results/pathways/curated300genes_at.least.1.larger.than.neglog.01.eps", height = 35, width = 30, units ="cm")
# ggsave("~/Dropbox/workspace/radiation/results/pathways/curated300genes_at.least.1.larger.than.neglog.01_distclust.tiff", height = 35, width = 30, units ="cm")
ggsave("~/Documents/dimitra/Workspace/RNA-Seq/radiation/results/pathways/curated300genes_at.least.1.larger.than.neglog.01_custom.ordered.pdf", height = 35, width = 30, units ="cm")
# ggsave("~/Dropbox/workspace/radiation/results/pathways/curated300genes_at.least.1.larger.than.neglog.01_custom.ordered.eps", height = 35, width = 30, units ="cm")

# my_palette <- colorRampPalette(c("green", "yellow","red" ))(n = 299)    # n  number of colors to be generated
# col_breaks = c(seq(0, 0.2, length=100),  # for red
#                seq(0.2,0.8,length=100),  # for yellow
#                seq(0.8,1,length=100))    # for green



#Updated heatmap: 300 genes, pvalue<0.01
library(scales)
ipadata <- read.csv("~/Documents/dimitra/Workspace/RNA-Seq/radiation/IPA/top300common_unique_lt.01.csv")

group_order = c("apoptosis p53 dependent", "p53 independent", "DNA damage.cellular stress", "inflammation", "other")
ipadata.ordered = ipadata[with(ipadata, order(factor(Group, levels=unique(group_order)))), ]
myPalette = c("indianred1" ,"orchid3", "red3", "salmon1", "steelblue3")
names(myPalette) = levels(ipadata.ordered$Group)

cols = c("apoptosis p53 dependent" = "indianred1", "p53 independent" = "orchid3", "DNA damage.cellular stress" = "red3", "inflammation" = "salmon1", "other" = "steelblue3")
shapes = c("apoptosis p53 dependent" = 15, "p53 independent" = 15, "DNA damage.cellular stress" = 15, "inflammation" = 15, "other" = 15)

ipd = ipadata.ordered %>%
  arrange(desc(Sorted)) %>%
  mutate(Pathway=factor(Pathway,levels=as.character(Pathway))) %>% 
  gather(Radiation,logPvalue,Proton:Gamma) %>%
  mutate(Radiation=factor(Radiation,levels=c("Gamma", "Proton"))) %>%
  mutate(logPvalue=replace(logPvalue, logPvalue<2, NA)) %>%
  mutate(Group=factor(Group,levels=unique(group_order))) 

ipd %>%   
  ggplot(aes(x = Radiation, y = Pathway, fill=logPvalue, colour=myPalette[ipd$Group])) +
  geom_tile() +
  scale_fill_gradientn(name = "-logPvalue", colours = c("darkolivegreen" ,"darkolivegreen4", "darkolivegreen3", "darkolivegreen1"), values = rescale(c(8,6,4,2)), labels = c("2", "4", "6","8"), breaks = c(2.05,4,6,8), na.value = "grey80") +
  scale_colour_manual(name='Groups', values=c("indianred1" ,"orchid3", "red3", "salmon1", "steelblue3"), labels=group_order) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 16, colour=myPalette[ipd$Group]),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.key.height = unit(1.5,"line"),
        legend.key.width = unit(1.5,"line"),
        panel.grid = element_blank()) +
  guides(colour = guide_legend(override.aes = list(fill=c("indianred1" ,"orchid3", "red3", "salmon1", "steelblue3"))))
  
ggsave("~/Documents/dimitra/Workspace/RNA-Seq/radiation/results/pathways/top300common_unique_lt.01_grouped.pdf", height = 35, width = 30, units ="cm")




# Doesn't work!
# Signed areas heatmap
sigdata <- read.table("~/Dropbox/workspace/radiation/results/signed_areas/all_signed_areas.txt", sep="\t", header = TRUE)

sigdata = sigdata[sigdata$experiment=="RNA-Seq", 1-2]

sigdata %>%
  arrange(desc(signed_area)) %>%
  mutate(id=factor(id,levels=as.character(id))) %>%
  
  ggplot(aes(x = experiment, y = id, fill=signed_area)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high="green", name="signed_area") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_blank(),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        panel.grid = element_blank())

ggsave("~/Dropbox/workspace/radiation/results/signed_areas/RNASeq.signed.areas.heatmap.tiff", height = 50, width = 30, units ="cm")
