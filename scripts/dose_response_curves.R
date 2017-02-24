###########################
###########################
### Dose response curves###

# radiation.data <- read.csv("~/Dropbox/workspace/radiation/results/venn_diagrams/common_genes_from_sense_trend_analysis_0.1.csv")
radiation.data <- read.csv("~/Dropbox/workspace/radiation/results/venn_diagrams/common.genes.sense.0.1.PORTvalues.txt", sep="\t", header = TRUE)

radiation.data %>%
  # gather(Sample_id,Expression,Gamma0.S11:Proton200.S10) %>%
  gather(Sample_id,Expression,Proton0.S1:Gamma200.S20) %>%
  mutate(Dose=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+","\\2",Sample_id),
         Dose = as.numeric(Dose),
         Radiation_type=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+","\\1",Sample_id)) %>% 
  
  ggplot(aes(x=Dose, y=Expression, color=Radiation_type)) + 
  geom_line(size=2) +
  # facet_wrap(~ID, ncol=4) +
  facet_wrap(Location~id, ncol=4, labeller = label_wrap_gen(multi_line=FALSE)) +
  scale_y_log10() +
  theme_bw(base_size = 32) + 
  theme(legend.position = c(0.89,0.04), legend.title = element_text(size = 28),
        legend.text = element_text(size = 26), legend.key.height = unit(1,"cm")) +
  labs(y="log gene expression", x="Dose (cGy)")

# ggsave("~/Dropbox/workspace/radiation/results/dose_response_curves/19common.genes.sense.0.1.drc.tiff", height = 50, width = 35, units ="cm")
# ggsave("~/Dropbox/workspace/radiation/results/dose_response_curves/19common.genes.sense.0.1.drc.eps", height = 50, width = 35, units ="cm")
ggsave("~/Dropbox/workspace/radiation/results/dose_response_curves/19common.genes.sense.0.1.fromPORT.drc.tiff", height = 60, width = 45, units ="cm")
ggsave("~/Dropbox/workspace/radiation/results/dose_response_curves/19common.genes.sense.0.1.fromPORT.drc.eps", height = 50, width = 35, units ="cm")




# http://stackoverflow.com/questions/24169675/multiple-colors-in-a-facet-strip-background
###### Change the facet bg color by group #####
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(scales)
library(gtable)
library(grid)
library(gridExtra)
library(cowplot)

# radiation.data <- read.csv("~/Dropbox/workspace/radiation/results/venn_diagrams/common_genes_from_sense_trend_analysis_0.1.csv")
radiation.data <- read.csv("~/Dropbox/workspace/radiation/results/venn_diagrams/common.genes.sense.0.1.PORTvalues.txt", sep="\t", header = TRUE)

# radiation.data$Location = factor(radiation.data$Location, levels=unique(c('Cytoplasm','Nucleus','Plasma membrane','Other')))

# radiation.data.ordered = radiation.data[with(radiation.data, order(Location, factor(ID, sort(order(Location)[ID])))), ]
# id_ordered = radiation.data.ordered$ID
radiation.data.ordered = radiation.data[with(radiation.data, order(Location, factor(id, sort(order(Location)[id])))), ]
id_ordered = radiation.data.ordered$id

rd = radiation.data.ordered %>%
  # gather(Sample_id,Expression,Gamma0.S11:Proton200.S10) %>%
  gather(Sample_id,Expression,Proton0.S1:Gamma200.S20) %>%
  mutate(Dose=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+","\\2",Sample_id),
         Dose = as.numeric(Dose),
         Radiation_type=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+","\\1",Sample_id), 
         # ID = factor(ID, levels=unique(id_ordered)))
         id = factor(id, levels=unique(id_ordered)))

p1 = ggplot(data = rd, aes(x=Dose, y=Expression, color=Radiation_type)) + 
  geom_line(size=2) +
  # facet_wrap(~ID, ncol=4) +
  facet_wrap(~id, ncol=4) +
  scale_y_log10() +
  # scale_x_log10() + # also change the xlabel
  theme_bw(base_size = 32) + 
  theme(legend.position = c(0.89,0.1), legend.title = element_text(size = 28),
        legend.text = element_text(size = 26), legend.key.height = unit(1,"cm")) +
  # labs(y="log Normalized gene expression", x="Dose (cGy)")
  labs(y="Log Gene expression", x="Dose (cGy)")
  # labs(y="Log Gene expression", x="Log Dose (cGy)")

dummy = ggplot(data = rd, aes(x=Dose, y=Expression)) + 
  # facet_wrap(~ID, ncol=4) +
  facet_wrap(~id, ncol=4) +
  geom_rect(aes(fill=Location), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
  scale_fill_manual(values=c("#4DAF4A", "#984EA3", "#FF7F00", "#F781BF")) +
  theme_minimal() + 
  theme(strip.text.x = element_text(size = 25), legend.position = c(0.89,0.04), legend.title = element_text(size = 28),
        legend.text = element_text(size = 26), legend.key.height = unit(1,"cm"))

##grab the legend only
# g_legend <- function(dummy){ 
#   tmp <- ggplot_gtable(ggplot_build(dummy)) 
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
#   legend <- tmp$grobs[[leg]] 
#   return(legend)} 
# facet.legend <- g_legend(dummy)

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(dummy)

gtable_select <- function (x, ...) 
{
  matches <- c(...)
  x$layout <- x$layout[matches, , drop = FALSE]
  x$grobs <- x$grobs[matches]
  x
}

panels <- grepl(pattern="panel", g2$layout$name)
strips <- grepl(pattern="strip-t", g2$layout$name)
g2$layout$t[panels] <- g2$layout$t[panels] - 1
g2$layout$b[panels] <- g2$layout$b[panels] - 1

new_strips <- gtable_select(g2, panels | strips)
grid.newpage()
grid.draw(new_strips)

gtable_stack <- function(g1, g2){
  g1$grobs <- c(g1$grobs, g2$grobs)
  g1$layout <- transform(g1$layout, z= z-max(z), name="g2")
  g1$layout <- rbind(g1$layout, g2$layout)
  g1
}

new_plot <- gtable_stack(g1, new_strips)
grid.newpage()
grid.draw(new_plot)

facet.legend <- g2$grobs[[which(sapply(g2$grobs, function(x) x$name) %in% "guide-box")]]

combined_plot = ggdraw() +
  draw_plot(new_plot, 0, 0, 1, 1) +
  draw_plot(facet.legend, 0.63, 0.05, .3, .3)

ggsave("~/Dropbox/workspace/radiation/results/dose_response_curves/19common.genes.sense.0.1.fromPORT.drc_colorgrouped.tiff", combined_plot, height = 60, width = 45, units ="cm")
ggsave("~/Dropbox/workspace/radiation/results/dose_response_curves/19common.genes.sense.0.1.fromPORT.drc_colorgrouped_wLogDose.tiff", combined_plot, height = 60, width = 45, units ="cm")


# tiff("~/Dropbox/workspace/radiation/results/dose_response_curves/19common.genes.sense.0.1.fromPORT.drc_colorgrouped_w_legend.tiff", width = 45, height = 60, units = "cm", res = 300)
# print(grid.arrange(new_plot, facet.legend, nrow = 2, ncol=2, widths = c(45, 2), heights = c(60, 2)))
# dev.off()

# ggsave("~/Dropbox/workspace/radiation/results/dose_response_curves/19common.genes.sense.0.1.drc_colorgrouped_.tiff", new_plot, height = 60, width = 45, units ="cm")
# ggsave("~/Dropbox/workspace/radiation/results/dose_response_curves/19common.genes.sense.0.1.drc_colorgrouped_.eps", new_plot, height = 60, width = 45, units ="cm")
# ggsave("~/Dropbox/workspace/radiation/results/dose_response_curves/19common.genes.sense.0.1.fromPORT.drc_colorgrouped_.tiff", new_plot, height = 60, width = 45, units ="cm")

