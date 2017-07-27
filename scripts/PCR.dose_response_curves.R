library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(scales)
library(gtable)
library(grid)
library(gridExtra)
library(cowplot)

radiation.data <- read.csv("~/Documents/dimitra/Workspace/RNA-Seq/radiation/results/dose_response_curves/RT-PCR_gene_expressions_edited_w.genesymbols.txt", sep="\t", header = TRUE)

radiation.data$Tissue = factor(radiation.data$Tissue, levels=unique(c('Aorta','Kidney','Liver','Lung','Heart')))

radiation.data.ordered = radiation.data[with(radiation.data, order(Tissue, factor(ID, sort(order(Tissue)[ID])))), ]
id_ordered = radiation.data.ordered$ID

rd = radiation.data.ordered %>% 
  gather(Sample_id,Expression,Proton0.S1:Gamma200.S20) %>%
  mutate(Dose=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+","\\2",Sample_id),
         Dose = as.numeric(Dose),
         Radiation_type=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+","\\1",Sample_id),
         Gene=gsub("^(.+)_(.+)","\\1",ID),
         ID = factor(ID, levels=unique(id_ordered)))

p1 = ggplot(data = rd, aes(x=Dose, y=Expression, color=Radiation_type)) + 
  geom_line(size=2) +
  facet_wrap(~ID, ncol=4) +
  scale_y_log10() +
  theme_bw(base_size = 32) + 
  theme(legend.position = c(0.9,0.14), legend.title = element_text(size = 28),
        legend.text = element_text(size = 26), legend.key.height = unit(1,"cm")) +
  labs(y="Log Gene expression", x="Dose (cGy)")

dummy = ggplot(data = rd, aes(x=Dose, y=Expression)) + 
  facet_wrap(~ID, ncol=4) + 
  geom_rect(aes(fill=Tissue), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
  scale_fill_manual(values=c("#4DAF4A", "#984EA3", "#FF7F00", "#F781BF", "#E69F00")) +
  theme_minimal() + 
  theme(strip.text.x = element_text(size = 25), legend.position = c(0.88,0.18), legend.title = element_text(size = 28),
        legend.text = element_text(size = 26), legend.key.height = unit(1,"cm"))

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
  draw_plot(facet.legend, 0.6, 0.04, .3, .3)

ggsave("~/Documents/dimitra/Workspace/RNA-Seq/radiation/results/dose_response_curves/approved/RT-PCR_gene_expressions.drc_colorgrouped.pdf", combined_plot, height = 60, width = 45, units ="cm")
