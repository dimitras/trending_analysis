####################################################################################
####################################################################################
# expression ratios histograms FOR ALL DATA

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(broom)

radiation.data <- read.table("~/Documents/dimitra/Workspace/RNA-Seq/radiation/results/expression_ratios_histograms/FINAL_master_list_of_gene_counts_MIN.sense.radiation_renamed.sorted_ratios.txt", header = TRUE)
# radiation.data <- read.table("~/Dropbox/workspace/radiation/results/dose_response_curves/FINAL_master_list_of_gene_counts_MIN.antisense.radiation_renamed.sorted_ratios_updated.txt", header = TRUE)

rd = radiation.data %>%
  gather(Sample_id,Ratio,Proton5.S2.Proton0.S1:Gamma200.S20.Gamma0.S11) %>%
  mutate(Dose=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\2",Sample_id),
         Dose = factor(paste(Dose,"cGy"), levels = c("5 cGy","10 cGy","25 cGy","50 cGy","75 cGy","100 cGy","125 cGy","150 cGy","200 cGy")),
         Radiation_type=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\1",Sample_id))

means = aggregate(Ratio~Dose+Radiation_type, rd, mean)
qs = aggregate(Ratio~Dose+Radiation_type, rd, fivenum)

rd %>%  
  ggplot(aes(x=Ratio)) + geom_histogram(binwidth=0.001, aes(fill=Radiation_type), position = "identity") + #for 1000 bins
    geom_vline(data = means, aes(xintercept=Ratio, colour=Radiation_type), size =1) +
    geom_vline(data = qs, aes(xintercept=Ratio[,2], colour=Radiation_type), size =1) +
    geom_vline(data = qs, aes(xintercept=Ratio[,4], colour=Radiation_type), size =1) +
    geom_freqpoly(aes(color=Radiation_type), size=5, alpha = 0.7) + 
    scale_x_log10(breaks = c(0.2, 1, 7)) +
    coord_cartesian(xlim=c(0.2, 7)) +
    facet_wrap(~Dose, ncol = 3) +
    theme_bw(base_size = 60) + 
    guides(colour = guide_legend(override.aes = list(size = 10)))+
    theme(legend.position = c(0.92,0.95), legend.title = element_text(size = 43), legend.spacing.y = unit(14, "cm"),
        legend.text = element_text(size = 48), legend.background=element_blank(), legend.key=element_rect(fill=NA), legend.key.size = unit(1, "cm"), legend.key.height = unit(2,"cm")) +
    labs(y="Number of genes", x="Log Expression ratio (irradiated/control)")

ggsave("~/Documents/dimitra/Workspace/RNA-Seq/radiation/results/expression_ratios_histograms/approved?/all_sense_expression_ratios_comparison_1000bins_wlog_3x3.pdf", height = 80, width = 80, units ="cm")
# ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/all_sense_expression_ratios_comparison_1000bins_wlog_3x3.eps", height = 80, width = 80, units ="cm")
# ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/antisense_radiation_ratios_comparison_1000bins_wlog.tiff", height = 15, width = 100, units ="cm")

# OLD
# ggplot(aes(x=Ratio)) + geom_histogram(binwidth=0.5, aes(fill=Radiation_type), position = "dodge") + #for 20 bins
# ggplot(aes(x=Ratio)) + geom_histogram(binwidth=0.03, aes(fill=Radiation_type), position = "dodge") + #for 300 bins
#
# geom_vline(aes(xintercept=mean(Ratio)), colour="green", size =1) +
# geom_vline(aes(xintercept=quantile(Ratio,0.25)), colour="purple", size =1) +
# geom_vline(aes(xintercept=quantile(Ratio,0.75)), colour="purple", size =1) +
# geom_vline(aes(xintercept=quantile(Ratio,0.90)), colour="purple", size =1) +

############################################################
############################################################
# CDF (cumulative distribution) FOR ALL DATA

radiation.data %>%
  gather(Sample_id,Ratio,Proton5.S2.Proton0.S1:Gamma200.S20.Gamma0.S11) %>%
  mutate(Dose=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\2",Sample_id),
         Dose = factor(paste(Dose,"cGy"), levels = c("5 cGy","10 cGy","25 cGy","50 cGy","75 cGy","100 cGy","125 cGy","150 cGy","200 cGy")),
         Radiation_type=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\1",Sample_id)) %>%
  
  ggplot(aes(x=Ratio, color=Radiation_type)) + 
  geom_step(stat="ecdf", size=3) +
  # scale_x_log10(breaks = c(0, 1, 5)) +
  scale_x_log10() +
  # coord_cartesian(xlim=c(0, 20)) +
  facet_wrap(~Dose, ncol = 3) +
  theme_bw(base_size = 40) + 
  theme(legend.position = c(0.93,0.05), legend.title = element_text(size = 32),
        legend.text = element_text(size = 30), legend.key.height = unit(1,"cm")) +
  labs(y="Empirical cumulative density function F(Expression ratio)", x="Expression ratio (irradiated/control)")

ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/all_sense_cdf_expression_ratios_comparison_wlog_3x3_zoomed.tiff", height = 60, width = 60, units ="cm")
# ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/all_sense_cdf_expression_ratios_comparison_wlog_3x3.eps", height = 80, width = 80, units ="cm")




####################################################################################
####################################################################################
# expression ratios histograms for genes with FDR<0.1

radiation.data <- read.table("~/Documents/dimitra/Workspace/RNA-Seq/radiation/results/expression_ratios_histograms/all_sense_0.1_ratios.txt", sep="\t", header = TRUE)
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




####################################################################################
####################################################################################
# CDF (cumulative distribution) for genes with FDR<0.1

radiation.data %>%
  gather(Sample_id,Ratio,Proton5.S2.baseline:Gamma200.S20.baseline) %>%
  mutate(Dose=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\2",Sample_id),
         Dose = factor(paste(Dose,"cGy"), levels = c("5 cGy","10 cGy","25 cGy","50 cGy","75 cGy","100 cGy","125 cGy","150 cGy","200 cGy")),
         Radiation_type=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\1",Sample_id)) %>%
  
  ggplot(aes(x=Ratio, color=Radiation_type)) + 
  geom_step(stat="ecdf", size=3) +
  # scale_x_log10() +
  facet_wrap(~Dose, ncol=3) + 
  theme_bw(base_size = 60) + 
  guides(colour = guide_legend(override.aes = list(size = 10)))+
  theme(legend.position = c(0.9,0.05), legend.title = element_text(size = 50), legend.spacing.y = unit(14, "cm"),
        legend.text = element_text(size = 50), legend.background=element_blank(), legend.key=element_rect(fill=NA), legend.key.size = unit(1, "cm"), legend.key.height = unit(2,"cm")) +
  labs(y="Empirical cumulative density function F(Expression ratio)", x="Expression ratio (irradiated/control)")

ggsave("~/Documents/dimitra/Workspace/RNA-Seq/radiation/results/expression_ratios_histograms/approved?/all_sense_0.1_cdf_expression_ratios_comparison_wlog_3x3.pdf", height = 80, width = 80, units ="cm")
# ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/all_sense_0.1_cdf_expression_ratios_comparison_wlog_3x3.eps", height = 80, width = 80, units ="cm")

# ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/all_PCR_0.1_cdf_expression_ratios_comparison_wlog_3x3.tiff", height = 80, width = 80, units ="cm")
# ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/all_PCR_0.1_cdf_expression_ratios_comparison_wlog_3x3.eps", height = 80, width = 80, units ="cm")




############################################################################################
############################################################################################
# Dose - Quartile expression ratio (treated/control) for ALL data AND genes with FDR<0.1

# radiation.data <- read.table("~/Documents/dimitra/Workspace/RNA-Seq/radiation/results/expression_ratios_histograms/FINAL_master_list_of_gene_counts_MIN.sense.radiation_renamed.sorted_ratios.txt", header = TRUE)
radiation.data <- read.table("~/Documents/dimitra/Workspace/RNA-Seq/radiation/results/expression_ratios_histograms/all_sense_0.1_ratios.txt", header = TRUE)
# FOR TESTING
# radiation.data <- read.table("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/test_ratios.txt", sep="\t", header = TRUE)

# t.test(radiation.data$Proton5.S2.baseline,radiation.data$Gamma5.S12.baseline, paired = TRUE)
# colMeans(subset(radiation.data,select = c(2:4)), na.rm = TRUE)

# wt = c(
#   wilcox.test(radiation.data$Proton5.S2.Proton0.S1,radiation.data$Gamma5.S12.Gamma0.S11)$p.value,
#   wilcox.test(radiation.data$Proton10.S3.Proton0.S1,radiation.data$Gamma10.S13.Gamma0.S11)$p.value,
#   wilcox.test(radiation.data$Proton25.S4.Proton0.S1,radiation.data$Gamma25.S14.Gamma0.S11)$p.value,
#   wilcox.test(radiation.data$Proton50.S5.Proton0.S1,radiation.data$Gamma50.S15.Gamma0.S11)$p.value,
#   wilcox.test(radiation.data$Proton75.S6.Proton0.S1,radiation.data$Gamma75.S16.Gamma0.S11)$p.value,
#   wilcox.test(radiation.data$Proton100.S7.Proton0.S1,radiation.data$Gamma100.S17.Gamma0.S11)$p.value,
#   wilcox.test(radiation.data$Proton125.S8.Proton0.S1,radiation.data$Gamma125.S18.Gamma0.S11)$p.value,
#   wilcox.test(radiation.data$Proton150.S9.Proton0.S1,radiation.data$Gamma150.S19.Gamma0.S11)$p.value,
#   wilcox.test(radiation.data$Proton200.S10.Proton0.S1,radiation.data$Gamma200.S20.Gamma0.S11)$p.value
# )
#7.862887e-02  5.793485e-06  7.495997e-48  1.538752e-29  5.490702e-28 5.975821e-136  2.123698e-07  2.522544e-38  1.443209e-27

wt = c(
wilcox.test(radiation.data$Proton5.S2.baseline,radiation.data$Gamma5.S12.baseline)$p.value,
wilcox.test(radiation.data$Proton10.S3.baseline,radiation.data$Gamma10.S13.baseline)$p.value,
wilcox.test(radiation.data$Proton25.S4.baseline,radiation.data$Gamma25.S14.baseline)$p.value,
wilcox.test(radiation.data$Proton50.S5.baseline,radiation.data$Gamma50.S15.baseline)$p.value,
wilcox.test(radiation.data$Proton75.S6.baseline,radiation.data$Gamma75.S16.baseline)$p.value,
wilcox.test(radiation.data$Proton100.S7.baseline,radiation.data$Gamma100.S17.baseline)$p.value,
wilcox.test(radiation.data$Proton125.S8.baseline,radiation.data$Gamma125.S18.baseline)$p.value,
wilcox.test(radiation.data$Proton150.S9.baseline,radiation.data$Gamma150.S19.baseline)$p.value,
wilcox.test(radiation.data$Proton200.S10.baseline,radiation.data$Gamma200.S20.baseline)$p.value
)
# 2.280435e-03 1.263370e-01 6.391689e-06 2.916258e-05 3.231177e-03 8.465336e-04 6.494024e-05 6.130252e-04 1.297570e-04

pvals=c()
for (pvalue in wt){
  if (pvalue < 0.05) {
    pvals <- c(pvals, "*")
  }
  else if (pvalue >= 0.05) {
    pvals <- c(pvals, "")
  }
}

radiation.data %>%
  # gather(Sample_id,Ratio,Proton5.S2.Proton0.S1:Gamma200.S20.Gamma0.S11) %>%
  gather(Sample_id,Ratio,Proton5.S2.baseline:Gamma200.S20.baseline) %>%
  mutate(Dose=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\2",Sample_id),
         Dose = as.numeric(Dose),
         Radiation_type=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\1",Sample_id)) %>%
  # order(Radiation_type, Dose, ID) %>%
  # group_by(ID) %>%
  # do(tidy(wilcox.test(Ratio ~ Radiation_type+Dose+ID, paired=TRUE))) %>%
  
  ggplot(aes(x=Dose, y=Ratio, color=Radiation_type)) +
  theme_bw(base_size = 60) + 
  theme(legend.position = c(0.9,0.04), legend.title = element_text(size = 50), legend.key.size = unit(1, "cm"), legend.spacing.y = unit(15, "cm"), legend.spacing.x = unit(10, "cm"),
        legend.text = element_text(size = 50), legend.key.height = unit(2,"cm"), legend.background=element_blank()) +
  # theme(legend.position = c(0.93,0.04), legend.title = element_text(size = 32), legend.key.size = unit(0.4, "cm"),
  #       legend.text = element_text(size = 30), legend.key.height = unit(1.5,"cm"), legend.background=element_blank()) +
  #separate graphs
  # labs(y="Median of expression radiated/control") + 
  # stat_summary(fun.y = median, geom = "line", position = "identity")
  # labs(y="Q1 of expression radiated/control") +
  # stat_summary(fun.y = function(z) { quantile(z,0.25) }, geom = "line", position = "identity")
  # labs(y="Q3 of expression radiated/control") +
  # stat_summary(fun.y = function(z) { quantile(z,0.75) }, geom = "line", position = "identity")

  #all in one
  labs(y="Quartiles of Expression ratio (irradiated/control)", x="Dose (cGy)") +
  stat_summary(fun.ymin = function(z) { quantile(z,0.25) },
               fun.ymax = function(z) { quantile(z,0.75) },
               fun.y = median,
               size=4, 
               position = position_dodge(width=3)) +
  # df1 <- data.frame(a = c(1, 1:3,3), b = c(39, 40, 40, 40, 39)) +
  # geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 2, y = 42, label = "*", size = 8)
  
  # annotate("text",x=6,y=1.12,label= pvals[1], size=12, angle=90) +
  # annotate("text",x=11,y=1.09,label= pvals[2], size=12, angle=90) +
  # annotate("text",x=26,y=1.12,label= pvals[3], size=12, angle=90) +
  # annotate("text",x=51,y=1.1,label= pvals[4], size=12, angle=90) +
  # annotate("text",x=76,y=1.18,label= pvals[5], size=12, angle=90) +
  # annotate("text",x=101,y=1.14,label= pvals[6], size=12, angle=90) +
  # annotate("text",x=126,y=1.12,label= pvals[7], size=12, angle=90) +
  # annotate("text",x=151,y=1.14,label= pvals[8], size=12, angle=90) +
  # annotate("text",x=201,y=1.12,label= pvals[9], size=12, angle=90)
  
  annotate("text",x=6,y=1.025,label= pvals[1], size=12, angle=90) +
  annotate("text",x=11,y=1.025,label= pvals[2], size=12, angle=90) +
  annotate("text",x=26,y=1.045,label= pvals[3], size=12, angle=90) +
  annotate("text",x=51,y=1.06,label= pvals[4], size=12, angle=90) +
  annotate("text",x=76,y=1.085,label= pvals[5], size=12, angle=90) +
  annotate("text",x=101,y=1.075,label= pvals[6], size=12, angle=90) +
  annotate("text",x=126,y=1.11,label= pvals[7], size=12, angle=90) +
  annotate("text",x=151,y=1.13,label= pvals[8], size=12, angle=90) +
  annotate("text",x=201,y=1.16,label= pvals[9], size=12, angle=90)
  
  # guides(color=guide_legend(overrride.aes=list(size=1, linetype=1)))

  # ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/Q1_expression_ratios.tiff", height = 80, width = 80, units ="cm")
  # ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/Q3_expression_ratios.tiff", height = 80, width = 80, units ="cm")
  # ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/Median_expression_ratios.tiff", height = 80, width = 80, units ="cm")

    ggsave("~/Documents/dimitra/Workspace/RNA-Seq/radiation/results/expression_ratios_histograms/approved?/summary_stats_all_expression_ratios.wilcox.pdf", height = 80, width = 80, units ="cm")
  # ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/summary_stats_all_expression_ratios.tiff", height = 80, width = 80, units ="cm")
  # ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/summary_stats_expression_ratios_0.1.tiff", height = 80, width = 80, units ="cm")
  # ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/summary_stats_expression_ratios_0.1.jpg", height = 80, width = 80, units ="cm")
  ggsave("~/Documents/dimitra/Workspace/RNA-Seq/radiation/results/expression_ratios_histograms/approved?/summary_stats_expression_ratios_0.1.wilcox.pdf", height = 80, width = 80, units ="cm")
  



########################################################################################################################
########################################################################################################################
# expression for genes with FDR<0.1 per dose (1..10), per radiation (P, G), per experiment (RNA-Seq, RT-PCR)
# compare the expression stats across RNA and PCR with boxplots to see the trend of the expression as the dose increases.

radiation.data <- read.table("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/common_genes_across_experiments_0.1.txt", sep="\t", header = TRUE)
# radiation.data <- read.table("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/union_genes_across_experiments_0.1.txt", sep="\t", header = TRUE)

#compare between experiments
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
# ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/union_genes_expression_across_experiments_0.1.tiff", height = 60, width = 60, units ="cm")


#RNA-Seq experiment only
radiation.data <- read.table("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/common_genes_across_experiments_0.1.txt", sep="\t", header = TRUE)
radiation.data %>%
  gather(Sample_id,Expression,RNA.Gamma0.S11:RNA.Proton200.S10) %>%
  mutate(Dose=gsub("^(RNA).(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\3",Sample_id),
         Dose = factor(Dose, levels=seq(from=0,to=200,by=5), labels=seq(from=0,to=200,by=5)),
         # Dose = as.numeric(Dose),
         Radiation_type=gsub("^(RNA).(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\2",Sample_id)) %>%
  
  ggplot(aes(x=Dose, y=Expression, fill=Radiation_type)) + 
  geom_boxplot(position = "dodge", na.rm = TRUE, width=10) +
  scale_x_discrete(limits=0:200, labels=c("0","5","10","25","50","75","100","125","150","175","200"), breaks=c(0,5,10,25,50,75,100,125,150,175,200)) +
  theme_bw(base_size = 40) + 
  theme(legend.position = c(0.93,0.1), legend.title = element_text(size = 32), legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 30), legend.key.height = unit(1.5,"cm"), legend.background=element_blank()) +
  labs(y="Normalized gene expression", x="Dose (cGy)")
  
  # ggplot(aes(x=interval, y=OR, colour=Drug)) + 
  # geom_point() + 
  # geom_line()+
  # geom_ribbon(aes(ymin=data$lower, ymax=data$upper), linetype=2, alpha=0.1)

ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/common_genes_expression_trend_0.1.tiff", height = 60, width = 60, units ="cm")


#RNA-Seq experiment only - shaded summary_stat
radiation.data <- read.table("~/Documents/dimitra/Workspace/RNA-Seq/radiation/results/expression_ratios_histograms/common_genes_across_experiments_0.1.txt", sep="\t", header = TRUE)

wt = c(
wilcox.test(radiation.data$RNA.Gamma0.S11,radiation.data$RNA.Proton0.S1)$p.value,
wilcox.test(radiation.data$RNA.Gamma5.S12,radiation.data$RNA.Proton5.S2)$p.value,
wilcox.test(radiation.data$RNA.Gamma10.S13,radiation.data$RNA.Proton10.S3)$p.value,
wilcox.test(radiation.data$RNA.Gamma25.S14,radiation.data$RNA.Proton25.S4)$p.value,
wilcox.test(radiation.data$RNA.Gamma50.S15,radiation.data$RNA.Proton50.S5)$p.value,
wilcox.test(radiation.data$RNA.Gamma75.S16,radiation.data$RNA.Proton75.S6)$p.value,
wilcox.test(radiation.data$RNA.Gamma100.S17,radiation.data$RNA.Proton100.S7)$p.value,
wilcox.test(radiation.data$RNA.Gamma125.S18,radiation.data$RNA.Proton125.S8)$p.value,
wilcox.test(radiation.data$RNA.Gamma150.S19,radiation.data$RNA.Proton150.S9)$p.value,
wilcox.test(radiation.data$RNA.Gamma200.S20,radiation.data$RNA.Proton200.S10)$p.value
)

pvals=c()
for (pvalue in wt){
  if (pvalue < 0.05) {
    pvals <- c(pvals, "*")
  }
  else if (pvalue >= 0.05) {
    pvals <- c(pvals, "")
  }
}

radiation.data %>%
  gather(Sample_id,Expression,RNA.Gamma0.S11:RNA.Proton200.S10) %>%
  mutate(Dose=gsub("^(RNA).(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\3",Sample_id),
         Dose = factor(Dose, levels=seq(from=0,to=200,by=5), labels=seq(from=0,to=200,by=5)),
         Radiation_type=gsub("^(RNA).(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\2",Sample_id)) %>%
  
  ggplot(aes(x=Dose, y=Expression, group=Radiation_type, color=Radiation_type)) + # fill=Radiation_type)) + 
  stat_summary(fun.ymin = function(z) { quantile(z,0.25) },
               fun.ymax = function(z) { quantile(z,0.75) },
               fun.y = median, size=2, geom = "smooth", aes(fill=Radiation_type), alpha=0.3) +
  scale_x_discrete(limits=0:200, labels=c("0","5","10","25","50","75","100","125","150","175","200"), breaks=c(0,5,10,25,50,75,100,125,150,175,200)) +
  theme_bw(base_size = 60) + 
  guides(colour = guide_legend(override.aes = list(size = 10)))+
  theme(legend.position = c(0.85,0.1), legend.title = element_text(size = 50), legend.spacing.y = unit(5, "cm"),
        legend.text = element_text(size = 50), legend.background=element_blank(), legend.key=element_rect(fill=NA), legend.key.size = unit(2, "cm"), legend.key.height = unit(2,"cm")) +
  labs(y="Normalized gene expression", x="Dose (cGy)") +
  annotate("text",x=100,y=1.45,label= pvals[7], size=40)

ggsave("~/Documents/dimitra/Workspace/RNA-Seq/radiation/results/expression_ratios_histograms/approved?/common_genes_expression_trend_0.1_shaded_w.wilcox.pdf", height = 60, width = 60, units ="cm")












# violin plot

radiation.data <- read.table("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/common_genes_across_experiments_0.1.txt", sep="\t", header = TRUE)

# median.quartile <- function(radiation.data){
#   out <- quantile(radiation.data, probs = c(0.25,0.5,0.75))
#   names(out) <- c("ymin","y","ymax")
#   return(out) 
# }

source("~/Dropbox/workspace/radiation/scripts/geom_flat_violin.R")

radiation.data %>%
  gather(Sample_id,Expression,RNA.Gamma0.S11:RNA.Proton200.S10) %>%
  mutate(Dose=gsub("^(RNA).(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\3",Sample_id),
         Dose = factor(Dose, levels=seq(from=0,to=200,by=5), labels=seq(from=0,to=200,by=5)),
         Radiation_type=gsub("^(RNA).(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\2",Sample_id)) %>%
  
  ggplot(aes(x=Dose, y=Expression, fill=Radiation_type)) + 
  # geom_violin(aes(y=Expression, fill=Radiation_type), draw_quantiles=c(0.25, 0.5, 0.75)) +
  geom_flat_violin(aes(y=Expression, fill=Radiation_type)) +
  # geom_boxplot(width=0.1) +
  labs(y="Normalized gene expression", x="Dose (cGy)") +
  theme_bw(base_size = 40) + 
  theme(legend.position = c(0.93,0.1), legend.title = element_text(size = 32), legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 30), legend.key.height = unit(1.5,"cm"), legend.background=element_blank()) 
  
ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/common_genes_expression_trend_0.1_violin.tiff", height = 60, width = 60, units ="cm")  
  

require(vioplot)
require(devtools)
require(digest)
source("~/Dropbox/workspace/radiation/scripts/vioplot2.R")

rd = radiation.data %>%
  gather(Sample_id,Expression,RNA.Gamma0.S11:RNA.Proton200.S10) %>%
  mutate(Dose=gsub("^(RNA).(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\3",Sample_id),
         Dose = factor(Dose, levels=seq(from=0,to=200,by=5), labels=seq(from=0,to=200,by=5)),
         Radiation_type=gsub("^(RNA).(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\2",Sample_id))


plot(x=NULL, y=NULL, xlim = c(0, 200), ylim=c(min(rd), max(rd)), type="n", ann=FALSE, axes=F) +
axis(1, at=c(0,200), labels=c(0,200)) +
axis(2)


for (i in unique(rd$Dose)) {
  for (j in unique(rd$Radiation_type)){
    vioplot2(values[which(rd$Dose == i & rd$Radiation_type == j)],
             at = ifelse(i == "0", 1, 2),
             side = ifelse(j == 1, "left", "right"),
             col = ifelse(j == 1, "purple", "lightblue"),
             add = T)
  }
}

rd = rd %>%
title("Violin plot", xlab="Dose") +
legend("bottomright", fill = c("purple", "lightblue"),
       fill = Radiation_type, box.lty=0)
