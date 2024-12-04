# all the R code used in the somatic mutations manuscript for plotting and analysis
library(table1)
library(ggplot2)
library(dplyr)
library(ggstatsplot)
library(tidyverse)
library(reshape2)
library(Polychrome)
library(ggpubr)
library(pals)
library(lme4) # for multilevel models
library(compositions)
library(tidyverse) # for data manipulation and plots
library(sjmisc)

# three different input dataframes, All mutations, all LP/P mutations and the filtered mutations which are corrected for gene length.
df_somatic_mutations <- read_csv('/Users/iwan/Research/Somatic_mutations/output_data/Somatic_mutations_with_controls.csv')
df_pathogenic <- read_csv('/Users/iwan/Research/Somatic_mutations/output_data/Pathogenic_somatic_mutations_with_controls.csv')
df_filtered_pathogenic <- read_csv('/Users/iwan/Research/Somatic_mutations/output_data/Pathogenic_somatic_mutations_with_controls_filtered.csv')





# plotting healthy vs controls pathogenic only, also done for the full dataframe
som_cd <- row_sums(df_pathogenic[df_pathogenic['Diagnosis'] != 'Control',2:28916],n=5)
som_control <-  row_sums(df_pathogenic[df_pathogenic['Diagnosis'] == 'Control',2:28916],n=5)


cd_sums <- as.data.frame(som_cd$rowsums)
colnames(cd_sums) <- c('Summed_mutations')
cd_sums$group <- 'IBD'

control_sums <- as.data.frame(som_control$rowsums)
colnames(control_sums) <- c('Summed_mutations')
control_sums$group <- 'Control'

data <- rbind(cd_sums, control_sums)
x <- ggplot(data, aes(x = factor(group, level=c('IBD', 'Control')), y = Summed_mutations,fill=group) ) +
  geom_violin(width=0.5,adjust=1.3, alpha=1 , color ='white') +
  geom_boxplot(width=0.05, alpha=0.7, color='white', alpha=0.7, outlier.colour ='black') + #, color = c("#339900", "#663399"

  xlab("Diagnosis") +
  ylab("Number of somatic mutations per biopsy") +
  ggtitle("Detected LP somatic mutations per biopsy split over diagnosis") +
  theme_classic() + scale_fill_identity() +
  scale_fill_manual(values=c("#339900", "#FF6633"))   # Color for violin plot

x <- x + scale_x_discrete(labels = c(paste0("IBD, n=", length(cd_sums$Summed_mutations)),
                                      paste0("Control, n=", length(control_sums$Summed_mutations))))
pval <- wilcox.test(data[data$group == 'IBD',]$Summed_mutations,  data[data$group == 'Control',]$Summed_mutations)
x <- x +
  geom_signif(comparisons = list(c("IBD", "Control")), annotations = round(pval$p.value, 10),
              map_signif_level = TRUE, step_increase=0.1)
x
x <- x +
  geom_signif(comparisons = list(c("IBD", "Control")),
              map_signif_level = TRUE, step_increase=0.1)
x




# load in gene panels:


gwas_ibd_genes = read_csv('/Users/iwan/Research/Somatic_mutations/metadata/paper_gwas_genes.txt', col_names=FALSE)
veo_ibd_genes = read_csv('/Users/iwan/Research/Somatic_mutations/metadata/veo_ibd_genes.txt', col_names=FALSE)
PID_ibd_genes = read_csv('/Users/iwan/Research/Somatic_mutations/metadata/immunodeficiency_genes.txt', col_names=FALSE)
loci_gwas_ibd_genes = read_csv('/Users/iwan/Research/Somatic_mutations/metadata/IBD_loci_gene_names_uniq.txt', col_names=FALSE)



# check per panel which genes are present in the data to prevent errors
present_gwas_genes <- c()
for (gene in gwas_ibd_genes$X1){
  if (is.element(gene,colnames(df_pathogenic))  ){
    present_gwas_genes <- c(present_gwas_genes, gene)
  }
}


present_veo_genes <- c()
for (gene in veo_ibd_genes$X1){
  if (is.element(gene,colnames(df_pathogenic))  ){
    present_veo_genes <- c(present_veo_genes, gene)
  }
}

present_pid_genes <- c()
for (gene in PID_ibd_genes$X1){
  if (is.element(gene,colnames(df_pathogenic))  ){
    present_pid_genes <- c(present_pid_genes, gene)
  }
}

present_non_genes <- c()
for (gene in non_ibd_genes$X1){
  if (is.element(gene,colnames(df_pathogenic))  ){
    present_non_genes <- c(present_non_genes, gene)
  }
}

# Take only the present veo-ibd genes
# this is repeated for all gene panels and all somatic mutations as well as only pathogenic
som_IBD <- row_sums(df_pathogenic[df_pathogenic['Diagnosis'] != 'Control', present_veo_genes],n=1)
som_control <-  row_sums(df_pathogenic[df_pathogenic['Diagnosis'] == 'Control',present_veo_genes],n=1)


IBD_sums <- as.data.frame(som_IBD$rowsums)
colnames(IBD_sums) <- c('Summed_mutations')
IBD_sums$group <- 'IBD'
control_sums <- as.data.frame(som_control$rowsums)
colnames(control_sums) <- c('Summed_mutations')
control_sums$group <- 'Control'


data <- rbind(IBD_sums, control_sums)


x <- ggplot(data, aes(x = factor(group, level=c('IBD', 'Control')), y = Summed_mutations, fill=group)) +

  geom_violin(width=0.5,adjust=1.3, alpha=1 , color ='white') +
  geom_boxplot(width=0.05, alpha=0.7, color='white', alpha=0.7, outlier.colour ='black') +#) +

  xlab("Diagnosis") +
  ylab("Number of somatic mutations per biopsy") +
  ggtitle("Detected LP somatic mutations per biopsy in VEO-IBD genes over diagnosis") +
  theme_classic()+ scale_fill_identity() +
  scale_fill_manual(values=c("#339900", "#FF6633"))

x <- x + scale_x_discrete(labels = c(paste0("IBD, n=", length(IBD_sums$Summed_mutations)),
                                      paste0("Control, n=", length(control_sums$Summed_mutations))))


pval <- wilcox.test(data[data$group == 'IBD',]$Summed_mutations,  data[data$group == 'Control',]$Summed_mutations)

x <- x +
  geom_signif(comparisons = list(c("IBD", "Control")), annotations = round(pval$p.value, 7),
              map_signif_level = TRUE, step_increase=0.1)
x


x <- x +
  geom_signif(comparisons = list(c("IBD", "Control")),
              map_signif_level = TRUE, step_increase=0.1)
x




# here we plot the mutations within IBD only
# take gene panel divide number of mutations bu total number of genes in the panel
# repeated for all mutations and pathogenic only.
veo_ibd_subset <- rowSums(all_som_mut_df_no_control[,present_veo_genes], na.rm = TRUE) / length(veo_ibd_genes$X1)
pid_subset <-  rowSums(all_som_mut_df_no_control[,present_pid_genes], na.rm = TRUE) / length(PID_ibd_genes$X1)
gwas_ibd_subset <-  rowSums(all_som_mut_df_no_control[,present_gwas_genes], na.rm = TRUE) / length(gwas_ibd_genes$X1)
non_ibd_subset <- rowSums(all_som_mut_df_no_control[,present_non_genes], na.rm = TRUE) / 16688
loci_gwas_ibd_subset <-  rowSums(all_som_mut_df_no_control[,present_loci_gwas_genes], na.rm = TRUE) / length(loci_gwas_ibd_genes$X1)
data <- data.frame(
    Group = c(rep(c("VEO-IBD Genes"), each=length(veo_ibd_subset)), rep(c("PID Genes"), each = length(pid_subset)),
              rep(c("GWAS Candidate Genes"), each = length(gwas_ibd_subset)), rep(c("Non-IBD Genes"), each = length(non_ibd_subset)), rep(c("GWAS Loci Genes"), each = length(loci_gwas_ibd_subset))),
    Mutations = c(veo_ibd_subset, pid_subset, gwas_ibd_subset, non_ibd_subset, loci_gwas_ibd_subset)
)


x <- ggplot(data, aes(x = factor(Group, level=c('VEO-IBD Genes', 'PID Genes', 'GWAS Candidate Genes', 'GWAS Loci Genes', 'Non-IBD Genes')), y = Mutations, fill=Group)) +

  geom_violin(width=0.5,adjust=1.3, alpha=1 , color ='white') +
  geom_boxplot(width=0.05, alpha=0.7, color='white', outlier.colour ='black') +#) +

  xlab("Gene set") +
  ylab("Log number of LP mutations per base in gene set") +
  ggtitle("Somatic LP mutations per gene set") +
  theme_classic() +
  scale_y_continuous(trans='log10') + theme(axis.text.x = element_text(angle = 45, vjust = 0.6, hjust=0.5))


x <- x + scale_x_discrete(labels = c(paste0("VEO-IBD Genes\n n=", length(veo_ibd_genes$X1)),
                                      paste0("PID Genes\n n=", length(PID_ibd_genes$X1)),
                                      paste0("GWAS Candidate Genes\n n=", length(gwas_ibd_genes$X1)),
                                     paste0("GWAS Loci Genes\n n=", length(loci_gwas_ibd_genes$X1)),
                                      paste0("Non-IBD Genes\n n=", '16688')))

x

# p values are calculated and then manually inserted for plotting purposes
pval <- wilcox.test(data[data$Group == 'VEO-IBD Genes',]$Mutations,  data[data$Group == 'Non-IBD Genes',]$Mutations)

x <- x +
  geom_signif(comparisons = list(c('VEO-IBD Genes', 'Non-IBD Genes')), annotations = '< 2.2e-16',
              map_signif_level = TRUE, step_increase=0.1, vjust=0, margin_top=0.1)

pval <- wilcox.test(data[data$Group == 'Loci GWAS-IBD Genes',]$Mutations,  data[data$Group == 'Non-IBD Genes',]$Mutations)

x <- x +
  geom_signif(comparisons = list(c('Loci GWAS-IBD Genes', 'Non-IBD Genes')), annotations = '< 2.2e-16',
              map_signif_level = TRUE, step_increase=0.1, vjust=0, margin_top=0.3)


pval <- wilcox.test(data[data$Group == 'PID Genes',]$Mutations,  data[data$Group == 'Non-IBD Genes',]$Mutations)

x <- x +
  geom_signif(comparisons = list(c('PID Genes', 'Non-IBD Genes')), annotations = '< 2.2e-16',
              map_signif_level = TRUE, step_increase=0.1, vjust=0, margin_top=0.2)

pval <- wilcox.test(data[data$Group == 'GWAS-IBD Genes',]$Mutations,  data[data$Group == 'Non-IBD Genes',]$Mutations)

x <- x +
  geom_signif(comparisons = list(c('GWAS-IBD Genes', 'Non-IBD Genes')), annotations = '< 2.2e-16',
              map_signif_level = TRUE, step_increase=0.1, vjust=0, margin_top=0.4)



## here we check the association with our selected cancer genes:
cancer_gene_list <- c('APC', 'TP53', 'TTN', 'KRAS', 'MUC16','SYNE1', 'FAT4', 'PIK3CA', 'OBSCN', 'ZFHX4', 'RYR2', 'CSMD3')

data_age <- all_som_mut_df_no_control$Age_at_biopsy
data_cancer_gene <- as.data.frame(data_age)
data_cancer_gene$mutations_counts <-  rowSums(all_som_mut_df_no_control[,cancer_gene_list], na.rm = TRUE)
data_cancer_gene <- drop_na(data_cancer_gene)


lm_eqn <- function(df){
    m <- lm(mutations_counts ~ data_age, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
         list(a = format(unname(coef(m)[1]), digits = 2),
              b = format(unname(coef(m)[2]), digits = 2),
             r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
}
x <- ggplot(data_cancer_gene, aes(x=data_age, y=mutations_counts)) +
    geom_point(size=2, color='blue' ,fill = "lightblue", alpha=0.7) +
    geom_smooth(method=lm, color='darkred') +
    theme_bw()+
    labs(title = "Somatic P/LP mutations in cancer genes vs age at biopsy",
        x = "Age at biopsy",
        y = "Number of somatic P/LP mutations per biopsy")+
  geom_text(x = 50, y = 0.6, label = lm_eqn(data_cancer_gene), parse = TRUE)
x



## Emmeans plotting code:
library(emmeans)
loci_inflam.lm <- lm('loci_ibd ~Inflammation +  Sequencing_Year + Location_rough +  (1|umcg_id)', data=all_som_mut_df_no_control, na.action = na.exclude)
loci_inflam.emm.s <- emmeans(loci_inflam.lm, data=all_som_mut_df_no_control, 'Inflammation')

gwas_inflam.lm <- lm('gwas_ibd ~Inflammation +  Sequencing_Year + Location_rough +  (1|umcg_id)', data=all_som_mut_df_no_control, na.action = na.exclude)
gwas_inflam.emm.s <- emmeans(gwas_inflam.lm, data=all_som_mut_df_no_control, 'Inflammation')

veo_inflam.lm <- lm('veo_ibd ~Inflammation +  Sequencing_Year + Location_rough +  (1|umcg_id)', data=all_som_mut_df_no_control, na.action = na.exclude)
veo_inflam.emm.s <- emmeans(veo_inflam.lm, data=all_som_mut_df_no_control, 'Inflammation')

pid_inflam.lm <- lm('pid_ibd ~Inflammation +  Sequencing_Year + Location_rough +  (1|umcg_id)', data=all_som_mut_df_no_control, na.action = na.exclude)
pid_inflam.emm.s <- emmeans(pid_inflam.lm, data=all_som_mut_df_no_control, 'Inflammation')


plot(loci_inflam.emm.s, comparisons = TRUE) +
  theme_bw() +
  ylab('Inflammation status') +
  ggtitle('GWAS Loci genes')

plot(gwas_inflam.emm.s, comparisons = TRUE) +
  theme_bw() +
  ylab('Inflammation status') +
  ggtitle('GWAS Candidate genes')

plot(veo_inflam.emm.s, comparisons = TRUE) +
  theme_bw() +
  ylab('Inflammation status') +
  ggtitle('VEO-IBD genes')

plot(pid_inflam.emm.s, comparisons = TRUE) +
  theme_bw() +
  ylab('Inflammation status') +
  ggtitle('PID genes')

# these individual plots are then combined in illustrator








## this code is used for the gene mutation plots:
## Jak 1 as an example

# first we define the location of each mutation and the counts
mut <- c(56, 143, 150, 223, 325, 339, 346, 591, 658, 692, 755, 777, 893, 995)
counts <- c(1, 1, 1, 1, 1,  18, 5, 1, 1, 1, 1, 1, 1, 1)

data <- data.frame(x=mut,y=counts)

# Plot as a lollypop plot
p <- ggplot(data, aes(x=x, y=y)) +
  geom_segment( aes(x=x, xend=x, y=0, yend=y),
                size=0.9, alpha=0.8) +
  geom_point(color='black',
             size=3) +
  theme_classic() +
  ylab('Number of observed mutations') +
  xlab('Protein position')+
  xlim(0,1154)+
  labs(title="JAK1")
p

# then we code the important regions in the genome from InterPro https://www.ebi.ac.uk/interpro/
genelength <- c(0:1154, 0:1154)
data=data.frame(x=genelength)
region_of_interest1 <- c(34:420,34:420,34:420)
data_region_of_interest1=data.frame(x=region_of_interest1)
region_of_interest2 <- c(437:544, 437:544,437:544)
data_region_of_interest2=data.frame(x=region_of_interest2)
region_of_interest3 <- c(583:855,583:855, 583:855, 875:1153, 875:1153, 875:1153)
data_region_of_interest3=data.frame(x=region_of_interest3)

# then we code this as a horizontal bar plot highlighting the important areas
ggplot(data, aes(x = x)) +
    geom_histogram(fill = "grey", bins=1000, binwidth=1, alpha = 1) +  #, aes(color = "JAK1")
    geom_histogram(data = data_region_of_interest1, fill = "red", bins=1000, binwidth=1, alpha = 1) + #, aes(color = "FERM domain")
    geom_histogram(data = data_region_of_interest2, fill = "lightblue",bins=1000, binwidth=1,  alpha = 1) +#, aes(color = "SH2 domain")
    geom_histogram(data = data_region_of_interest3, fill = "lightgreen",bins=1000, binwidth=1,  alpha = 1) +#, aes(color = "Protein Kinase domain")
  theme_classic()+
  ylim(0,10)+
  xlim(0,1154)

 # finally we manually overlap these two plots and annotate the mutations in illustrator




 ## here we check our model statistics for phenotypes affecting mutation counts:

model_df <- df_single_all
formula <- 'mutation_sum ~ Diagnosis+ Inflammation+ Location_rough + Sex + Age_at_biopsy +(1|umcg_id)'
model <- lmer(formula, data =model_df)

summary(model)


# and here for individual gene panels :


formula <- 'veo_ibd ~ Diagnosis+ Inflammation+ Location_rough + Sex + Age_at_biopsy +(1|umcg_id)'
model <- lmer(formula, data =model_df)
summary(model)


formula <- 'gwas_ibd ~ Diagnosis+ Inflammation+ Location_rough + Sex + Age_at_biopsy +(1|umcg_id)'
model <- lmer(formula, data =model_df)
summary(model)


formula <- 'loci_ibd ~ Diagnosis+ Inflammation+ Location_rough + Sex + Age_at_biopsy +(1|umcg_id)'
model <- lmer(formula, data =model_df)
summary(model)


formula <- 'pid_ibd ~ Diagnosis+ Inflammation+ Location_rough + Sex + Age_at_biopsy +(1|umcg_id)'
model <- lmer(formula, data =model_df)
summary(model)

formula <- 'non_ibd ~ Diagnosis+ Inflammation+ Location_rough + Sex + Age_at_biopsy +(1|umcg_id)'
model <- lmer(formula, data =model_df)
summary(model)
