library(PMCMR)
library(FSA)
library(reshape2)

############################################################################################
############################################################################################
# Kruskal Wallis Post Hoc test for genes with FDR<0.1
# Group by Radiation and Dose

# radiation.data <- read.table("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/FINAL_master_list_of_gene_counts_MIN.sense.radiation_renamed.sorted_ratios.txt", header = TRUE)
radiation.data <- read.table("~/Documents/dimitra/Workspace/RNA-Seq/radiation/results/expression_ratios_histograms/all_sense_0.1_ratios.txt", sep="\t", header = TRUE)

rd = radiation.data %>%
  gather(Sample_id,Ratio,Proton5.S2.baseline:Gamma200.S20.baseline) %>%
  mutate(Dose=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\2",Sample_id),
         Radiation_type=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\1",Sample_id),
         Dose = factor(Dose, levels=unique(Dose)),
         Group = gsub("^([[:alpha:]]+[[:digit:]]+)\\.S[[:digit:]]+.*","\\1",Sample_id),
         Group = factor(Group, levels=unique(Group)))

#Kruskal Wallis test
kt_all = kruskal.test(Ratio~Group,data=rd)$p.value

#Dunn test
dt = posthoc.kruskal.dunn.test(x=rd$Ratio, g=rd$Group, p.adjust.method="bonferroni")
dt_long <- melt(dt$p.value, measure.vars=c("5","10","25","50","75","100","125","150","200"))
colnames(dt_long) <- c("DoseY", "DoseX", "Pvalue")
dt_long$Pair = paste(dt_long$DoseX,"-",dt_long$DoseY)
dt_long = na.omit(dt_long) 

selected = c("Proton5 - Gamma5", "Proton10 - Gamma10", "Proton25 - Gamma25", "Proton50 - Gamma50", "Proton75 - Gamma75", "Proton100 - Gamma100", "Proton125 - Gamma125", "Proton150 - Gamma150", "Proton200 - Gamma200")
# dt_long$Selected = gsub("^(Proton[[:digit:]]+.*Gamma[[:digit:]]+)","\\1",dt_long$Pair)
dt_long_selected = dt_long[dt_long$Pair %in% selected,]

# dt_long %>%
dt_long_selected %>%  
  ggplot(aes(x=Pvalue, y=Pair)) +
  geom_point(size=6) +
  geom_vline(xintercept = 0.05) +
  # grid() +
  theme_bw(base_size = 20) +
  labs(y="Dose-pairs", x="Dunn's test adjusted p-value (Bonferroni)") +
  scale_x_continuous(breaks=c(0,.05,.25,.5,.75,1))

ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/kruskal.dunn.pvalues_w.selected.groups.tiff", height = 20, width = 40, units ="cm")


#OR with FSA
dt2 = dunnTest(Ratio ~ Group, data=rd,method="bonferroni")
dt2.res = dt2$res 
dt2.res.selected = dt2.res[dt2.res$Comparison %in% selected,] # for selected group


### Compact letter display
# When there are many p-values to evaluate, it is useful to condense a table of p-values to a compact letter display format. In the output, groups are separated by letters. Groups sharing the same letter are not significantly different. Compact letter displays are a clear and succinct way to present results of multiple comparisons
library(rcompanion)

cldList(comparison = dt2.res.selected$Comparison,
        p.value    = dt2.res.selected$P.adj,
        threshold  = 0.05)

cldList(comparison = dt2.res$Comparison,
        p.value    = dt2.res$P.adj,
        threshold  = 0.05)






############################################################################################
############################################################################################
# Kruskal Wallis Post Hoc test for genes with FDR<0.1
# Group by Dose 

# radiation.data <- read.table("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/FINAL_master_list_of_gene_counts_MIN.sense.radiation_renamed.sorted_ratios.txt", header = TRUE)
radiation.data <- read.table("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/all_sense_0.1_ratios.txt", sep="\t", header = TRUE)

rd = radiation.data %>%
  # gather(Sample_id,Ratio,Proton5.S2.Proton0.S1:Gamma200.S20.Gamma0.S11) %>%
  gather(Sample_id,Ratio,Proton5.S2.baseline:Gamma200.S20.baseline) %>%
  mutate(Dose=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\2",Sample_id),
         # Dose = as.numeric(Dose),
         Radiation_type=gsub("^(Proton|Gamma)([[:digit:]]+)\\.S[[:digit:]]+.*","\\1",Sample_id),
         Dose = factor(Dose, levels=unique(Dose))) #c("5","10","25","50","75","100","125","150","200")))

# kt = c(
#   kruskal.test(radiation.data$Proton5.S2.Proton0.S1,radiation.data$Gamma5.S12.Gamma0.S11)$p.value,
#   kruskal.test(radiation.data$Proton10.S3.Proton0.S1,radiation.data$Gamma10.S13.Gamma0.S11)$p.value,
#   kruskal.test(radiation.data$Proton25.S4.Proton0.S1,radiation.data$Gamma25.S14.Gamma0.S11)$p.value,
#   kruskal.test(radiation.data$Proton50.S5.Proton0.S1,radiation.data$Gamma50.S15.Gamma0.S11)$p.value,
#   kruskal.test(radiation.data$Proton75.S6.Proton0.S1,radiation.data$Gamma75.S16.Gamma0.S11)$p.value,
#   kruskal.test(radiation.data$Proton100.S7.Proton0.S1,radiation.data$Gamma100.S17.Gamma0.S11)$p.value,
#   kruskal.test(radiation.data$Proton125.S8.Proton0.S1,radiation.data$Gamma125.S18.Gamma0.S11)$p.value,
#   kruskal.test(radiation.data$Proton150.S9.Proton0.S1,radiation.data$Gamma150.S19.Gamma0.S11)$p.value,
#   kruskal.test(radiation.data$Proton200.S10.Proton0.S1,radiation.data$Gamma200.S20.Gamma0.S11)$p.value
# )
# [1] 1.107851e-189 7.519225e-136  8.741057e-25  1.132323e-32  2.621398e-45  7.894216e-83  2.182544e-25  0.000000e+00 5.328510e-317

# kt = c(
#   kruskal.test(radiation.data$Proton5.S2.baseline,radiation.data$Gamma5.S12.baseline)$p.value,
#   kruskal.test(radiation.data$Proton10.S3.baseline,radiation.data$Gamma10.S13.baseline)$p.value,
#   kruskal.test(radiation.data$Proton25.S4.baseline,radiation.data$Gamma25.S14.baseline)$p.value,
#   kruskal.test(radiation.data$Proton50.S5.baseline,radiation.data$Gamma50.S15.baseline)$p.value,
#   kruskal.test(radiation.data$Proton75.S6.baseline,radiation.data$Gamma75.S16.baseline)$p.value,
#   kruskal.test(radiation.data$Proton100.S7.baseline,radiation.data$Gamma100.S17.baseline)$p.value,
#   kruskal.test(radiation.data$Proton125.S8.baseline,radiation.data$Gamma125.S18.baseline)$p.value,
#   kruskal.test(radiation.data$Proton150.S9.baseline,radiation.data$Gamma150.S19.baseline)$p.value,
#   kruskal.test(radiation.data$Proton200.S10.baseline,radiation.data$Gamma200.S20.baseline)$p.value
# )
# [1] 0.5305351 0.4793370 0.3430635 0.4255500 0.4556526 0.3507277 0.3592291 0.3951031 0.4556526

kt_all = kruskal.test(Ratio~Dose,data=rd)$p.value
#Kruskal-Wallis chi-squared = 1643.7, df = 8, p-value < 2.2e-16
#5.13025e-217

#Dunn test
dt = posthoc.kruskal.dunn.test(x=rd$Ratio, g=rd$Dose, p.adjust.method="bonferroni")
# dun = posthoc.kruskal.dunn.test(x=rd$Ratio, g=rd$Dose, p.adjust.method="none")
library(reshape2)
dt_long <- melt(dt$p.value, measure.vars=c("5","10","25","50","75","100","125","150","200"))
colnames(dt_long) <- c("DoseY", "DoseX", "Pvalue")
dt_long$Pair = paste(dt_long$DoseX,"-",dt_long$DoseY)
dt_long = na.omit(dt_long) 

dt_long %>%
  ggplot(aes(x=Pvalue, y=Pair)) +
  geom_point(size=6) +
  geom_vline(xintercept = 0.05) +
  grid() +
  theme_bw(base_size = 30) +
  labs(y="Dose-pairs", x="Dunn's test adjusted p-value (Bonferroni)") +
  scale_x_continuous(breaks=c(0,.05,.25,.5,.75,1))

ggsave("~/Dropbox/workspace/radiation/results/expression_ratios_histograms/kruskal.dunn.pvalues.tiff", height = 50, width = 50, units ="cm")


#OR with FSA
dt = dunnTest(Ratio ~ Dose, data=rd,method="bonferroni")
dt.res = dt$res 

library(rcompanion)
cldList(comparison = dt.res$Comparison,
        p.value    = dt.res$P.adj,
        threshold  = 0.05)

