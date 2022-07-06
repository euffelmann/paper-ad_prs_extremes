

# set-up
rm(list = ls())
library(tidyverse)
library(data.table)
library(ggpubr)

# load data
data_list <- list()
alz_fpr_fpn_list <- list() # fpn = False positive number
files <- list.files(path = "simulation_results", pattern = glob2rx("run_summary_*.rds"), full.names = T)

##### main figure 1: prs_threshold = 0.05; fpr_threshold = 5e-8 ####

prs_threshold = 0.05; fpr_threshold = 5e-8

for (i in 1:length(files)) {
  
  run_summary <- readRDS(files[i])[[as.character(prs_threshold)]]
  alz_sumstats <- run_summary$gwas_alz_sumstats
  n_all_null_snps <- alz_sumstats %>% filter(str_detect(SNP, "null")) %>% nrow()
  n_null_snps_in_prs <- alz_sumstats %>% filter(str_detect(SNP, "null") & P <= prs_threshold) %>% nrow()
  n_null_snps_not_in_prs <- alz_sumstats %>% filter(str_detect(SNP, "null") & P > prs_threshold) %>% nrow()
  run_summary[["prs_extremes_fpr"]]$snp_set_number <- rep(x = c(n_all_null_snps, n_all_null_snps, 
                                                                n_null_snps_in_prs, n_null_snps_in_prs, 
                                                                n_null_snps_not_in_prs, n_null_snps_not_in_prs))
  run_summary[["prs_extremes_fpr"]]$fpn <- as.numeric(run_summary[["prs_extremes_fpr"]]$fpr) * run_summary[["prs_extremes_fpr"]]$snp_set_number
  data_list[[i]] <- run_summary[["prs_extremes_fpr"]]
  
  alz_fpr <- alz_sumstats %>% filter(str_detect(SNP, "null")) %>% {sum(.$P < fpr_threshold)} / n_all_null_snps
  alz_fpn <- alz_sumstats %>% filter(str_detect(SNP, "null")) %>% {sum(.$P < fpr_threshold)}
  alz_fpr_fpn_list[[i]] <- c(alz_fpr, alz_fpn)
  
}

fpr_dfr <- do.call("rbind", data_list)
fpr_dfr <- fpr_dfr %>%
  mutate(overlap = as.numeric(overlap),
         overlap2 = fct_reorder(paste0(overlap, "%"), overlap),
         fpr = as.numeric(fpr), 
         sig_threshold = as.numeric(sig_threshold),
         snp_set_name = case_when(snp_set_name == "all_null_snps" ~ "All null-SNPs",
                                  snp_set_name == "null_snps_in_prs" ~ "Null-SNPs used in PRS",
                                  snp_set_name == "null_snps_not_in_prs" ~ "Null-SNPs not used in PRS"))
fpr_dfr$fpr[is.na(fpr_dfr$fpr)] <- 0

n_all_null_snps <- fpr_dfr$snp_set_number[fpr_dfr$snp_set_name == "All null-SNPs"][1]
n_null_snps_in_prs <- fpr_dfr$snp_set_number[fpr_dfr$snp_set_name == "Null-SNPs used in PRS"][1]
n_null_snps_not_in_prs <- fpr_dfr$snp_set_number[fpr_dfr$snp_set_name == "Null-SNPs not used in PRS"][1]

p1 <- fpr_dfr %>%
  filter(sig_threshold == fpr_threshold,
         snp_set_name == "All null-SNPs") %>%
  ggplot(aes(x = overlap2, y = fpr)) +
  #geom_violin() +
  geom_boxplot(width=0.3, color="black", alpha=0.2) +
  geom_point(color="grey", alpha=0.4, size = 0.5) + 
  theme_classic() +
  scale_y_continuous(sec.axis = sec_axis(~ . * n_all_null_snps, name="False positive number")) +
  theme(
    legend.position="none",
    strip.text.x = element_text(size = 8),
    text=element_text(family="Charter"),
    strip.background = element_rect(
      color="#EBEBEB", fill="#EBEBEB", size=1.5, linetype="solid"),
    plot.title = element_text(hjust = 0.5, size = 10)
  ) +
  xlab("Sample overlap") +
  ylab("False positive rate") +
  ggtitle("All null-SNPs") +
  geom_hline(yintercept=5e-8, linetype="dashed", 
             color = "red", size=0.5) +
  geom_text(aes(1, 5e-8, label = 5e-8, colour = "red"),
            size = 3, vjust = -0.5, hjust = 0.5)

p2 <- fpr_dfr %>%
  filter(sig_threshold == fpr_threshold,
         snp_set_name == "Null-SNPs used in PRS") %>%
  ggplot(aes(x = overlap2, y = fpr)) +
  #geom_violin() +
  geom_boxplot(width=0.3, color="black", alpha=0.2) +
  geom_point(color="grey", alpha=0.4, size = 0.5) +
  scale_y_continuous(sec.axis = sec_axis(~ . * n_null_snps_in_prs, name="")) +
  theme_classic() +
  theme(
    legend.position="none",
    strip.text.x = element_text(size = 8),
    text=element_text(family="Charter"),
    strip.background = element_rect(
      color="#EBEBEB", fill="#EBEBEB", size=1.5, linetype="solid"),
    plot.title = element_text(hjust = 0.5, size = 10)
  ) +
  xlab("") +
  ylab("") +
  ggtitle("Null-SNPs used in PRS") +
  geom_hline(yintercept=fpr_threshold, linetype="dashed", 
             color = "red", size=0.5)

p3 <- fpr_dfr %>%
  filter(sig_threshold == fpr_threshold,
         snp_set_name == "Null-SNPs not used in PRS") %>%
  ggplot(aes(x = overlap2, y = fpr)) +
  #  geom_violin() +
  geom_boxplot(width=0.3, color="black", alpha=0.2) +
  geom_point(color="grey", alpha=0.4, size = 0.5) +
  scale_y_continuous(sec.axis = sec_axis(~ . * n_null_snps_not_in_prs, name="False positive number", breaks=c(0, 1))) +
  theme_classic() +
  theme(
    legend.position="none",
    strip.text.x = element_text(size = 8),
    text=element_text(family="Charter"),
    strip.background = element_rect(
      color="#EBEBEB", fill="#EBEBEB", size=1.5, linetype="solid"),
    plot.title = element_text(hjust = 0.5, size = 10)
  ) +
  xlab("Sample overlap") +
  ylab("False positive rate") +
  ggtitle("Null-SNPs not used in PRS") +
  geom_hline(yintercept=fpr_threshold, linetype="dashed", 
             color = "red", size=0.5)

alz_fpr_fpn_dfr <- as.data.frame(do.call("rbind", alz_fpr_fpn_list))
names(alz_fpr_fpn_dfr) <- c("fpr", "fpn")
p4 <- alz_fpr_fpn_dfr %>%
  ggplot(aes(y = fpr)) +
  #geom_violin(aes(x = 0)) +
  geom_boxplot(width=0.3, color="black", alpha=0.2) +
  scale_x_discrete() +
  geom_point(aes(x = 0), color="grey", alpha = 1, size = 1) +
  scale_y_continuous(sec.axis = sec_axis(~ . * n_all_null_snps, name="", breaks=c(0, 1))) +
  theme_classic() +
  theme(
    legend.position="none",
    strip.text.x = element_text(size = 8),
    text=element_text(family="Charter"),
    strip.background = element_rect(
      color="#EBEBEB", fill="#EBEBEB", size=1.5, linetype="solid"),
    plot.title = element_text(hjust = 0.5, size = 10)
  ) +
  xlab("") +
  ylab("") +
  ggtitle("GWAS of AD") +
  geom_hline(yintercept=fpr_threshold, linetype="dashed", 
             color = "red", size=0.5)

figure <- ggarrange(p1, p2, p3, p4,
                    labels = c("A", "B", "C", "D"),
                    ncol = 2, nrow = 2, font.label = list(size = 10, family="Charter"), hjust = c(-0.5, -3.5))

ggsave(plot = figure, filename = paste0("plots/prs_threshold", prs_threshold, "_fpr_threshold", fpr_threshold, ".png"),
       width = 15, height = 15, units = "cm", dpi = 500, bg = "white")


##### suppl figure 1: prs_threshold = 0.05; fpr_threshold = 0.05 ####

prs_threshold = 0.05; fpr_threshold = 0.05

for (i in 1:length(files)) {
  
  run_summary <- readRDS(files[i])[[as.character(prs_threshold)]]
  alz_sumstats <- run_summary$gwas_alz_sumstats
  n_all_null_snps <- alz_sumstats %>% filter(str_detect(SNP, "null")) %>% nrow()
  n_null_snps_in_prs <- alz_sumstats %>% filter(str_detect(SNP, "null") & P <= prs_threshold) %>% nrow()
  n_null_snps_not_in_prs <- alz_sumstats %>% filter(str_detect(SNP, "null") & P > prs_threshold) %>% nrow()
  run_summary[["prs_extremes_fpr"]]$snp_set_number <- rep(x = c(n_all_null_snps, n_all_null_snps, 
                                                                n_null_snps_in_prs, n_null_snps_in_prs, 
                                                                n_null_snps_not_in_prs, n_null_snps_not_in_prs))
  run_summary[["prs_extremes_fpr"]]$fpn <- as.numeric(run_summary[["prs_extremes_fpr"]]$fpr) * run_summary[["prs_extremes_fpr"]]$snp_set_number
  data_list[[i]] <- run_summary[["prs_extremes_fpr"]]
  
  alz_fpr <- alz_sumstats %>% filter(str_detect(SNP, "null")) %>% {sum(.$P < fpr_threshold)} / n_all_null_snps
  alz_fpn <- alz_sumstats %>% filter(str_detect(SNP, "null")) %>% {sum(.$P < fpr_threshold)}
  alz_fpr_fpn_list[[i]] <- c(alz_fpr, alz_fpn)  
}

fpr_dfr <- do.call("rbind", data_list)
fpr_dfr <- fpr_dfr %>%
  mutate(overlap = as.numeric(overlap),
         overlap2 = fct_reorder(paste0(overlap, "%"), overlap),
         fpr = as.numeric(fpr), 
         sig_threshold = as.numeric(sig_threshold),
         snp_set_name = case_when(snp_set_name == "all_null_snps" ~ "All null-SNPs",
                                  snp_set_name == "null_snps_in_prs" ~ "Null-SNPs used in PRS",
                                  snp_set_name == "null_snps_not_in_prs" ~ "Null-SNPs not used in PRS"))
fpr_dfr$fpr[is.na(fpr_dfr$fpr)] <- 0

n_all_null_snps <- fpr_dfr$snp_set_number[fpr_dfr$snp_set_name == "All null-SNPs"][1]
n_null_snps_in_prs <- fpr_dfr$snp_set_number[fpr_dfr$snp_set_name == "Null-SNPs used in PRS"][1]
n_null_snps_not_in_prs <- fpr_dfr$snp_set_number[fpr_dfr$snp_set_name == "Null-SNPs not used in PRS"][1]

p1 <- fpr_dfr %>%
  filter(sig_threshold == fpr_threshold,
         snp_set_name == "All null-SNPs") %>%
  ggplot(aes(x = overlap2, y = fpr)) +
  #geom_violin() +
  geom_boxplot(width=0.3, color="black", alpha=0.2) +
  geom_point(color="grey", alpha=0.4, size = 0.5) + 
  theme_classic() +
  scale_y_continuous(sec.axis = sec_axis(~ . * n_all_null_snps, name="False positive number")) +
  theme(
    legend.position="none",
    strip.text.x = element_text(size = 8),
    text=element_text(family="Charter"),
    strip.background = element_rect(
      color="#EBEBEB", fill="#EBEBEB", size=1.5, linetype="solid"),
    plot.title = element_text(hjust = 0.5, size = 10)
  ) +
  xlab("Sample overlap") +
  ylab("False positive rate") +
  ggtitle("All null-SNPs") +
  geom_hline(yintercept = fpr_threshold, linetype="dashed", 
             color = "red", size=0.5)

p2 <- fpr_dfr %>%
  filter(sig_threshold == fpr_threshold,
         snp_set_name == "Null-SNPs used in PRS") %>%
  ggplot(aes(x = overlap2, y = fpr)) +
  #geom_violin() +
  geom_boxplot(width=0.3, color="black", alpha=0.2) +
  geom_point(color="grey", alpha=0.4, size = 0.5) +
  scale_y_continuous(sec.axis = sec_axis(~ . * n_null_snps_in_prs, name="")) +
  theme_classic() +
  theme(
    legend.position="none",
    strip.text.x = element_text(size = 8),
    text=element_text(family="Charter"),
    strip.background = element_rect(
      color="#EBEBEB", fill="#EBEBEB", size=1.5, linetype="solid"),
    plot.title = element_text(hjust = 0.5, size = 10)
  ) +
  xlab("") +
  ylab("") +
  ggtitle("Null-SNPs used in PRS") +
  geom_hline(yintercept = fpr_threshold, linetype="dashed", 
             color = "red", size = 0.5)

p3 <- fpr_dfr %>%
  filter(sig_threshold == fpr_threshold,
         snp_set_name == "Null-SNPs not used in PRS") %>%
  ggplot(aes(x = overlap2, y = fpr)) +
  #  geom_violin() +
  geom_boxplot(width=0.3, color="black", alpha=0.2) +
  geom_point(color="grey", alpha=0.4, size = 0.5) +
  scale_y_continuous(sec.axis = sec_axis(~ . * n_null_snps_not_in_prs, name="False positive number")) +
  theme_classic() +
  theme(
    legend.position="none",
    strip.text.x = element_text(size = 8),
    text=element_text(family="Charter"),
    strip.background = element_rect(
      color="#EBEBEB", fill="#EBEBEB", size=1.5, linetype="solid"),
    plot.title = element_text(hjust = 0.5, size = 10)
  ) +
  xlab("Sample overlap") +
  ylab("False positive rate") +
  ggtitle("Null-SNPs not used in PRS") +
  geom_hline(yintercept = fpr_threshold, linetype = "dashed", 
             color = "red", size = 0.5)

alz_fpr_fpn_dfr <- as.data.frame(do.call("rbind", alz_fpr_fpn_list))
names(alz_fpr_fpn_dfr) <- c("fpr", "fpn")
p4 <- alz_fpr_fpn_dfr %>%
  ggplot(aes(y = fpr)) +
  #geom_violin(aes(x = 0)) +
  geom_boxplot(width = 0.3, color="black", alpha=0.2) +
  scale_x_discrete() +
  geom_point(aes(x = 0), color="grey", alpha = 1, size = 1) +
  scale_y_continuous(sec.axis = sec_axis(~ . * n_all_null_snps, name="")) +
  theme_classic() +
  theme(
    legend.position="none",
    strip.text.x = element_text(size = 8),
    text=element_text(family="Charter"),
    strip.background = element_rect(
      color="#EBEBEB", fill="#EBEBEB", size=1.5, linetype="solid"),
    plot.title = element_text(hjust = 0.5, size = 10)
  ) +
  xlab("") +
  ylab("") +
  ggtitle("GWAS of ALZ") +
  geom_hline(yintercept=fpr_threshold, linetype="dashed", 
             color = "red", size=0.5)

figure <- ggarrange(p1, p2, p3, p4,
                    labels = c("A", "B", "C", "D"),
                    ncol = 2, nrow = 2, font.label = list(size = 10, family="Charter"), hjust = c(-0.5, -3.5))

ggsave(plot = figure, filename = paste0("plots/prs_threshold", prs_threshold, "_fpr_threshold", fpr_threshold, ".png"),
       width = 15, height = 15, units = "cm", dpi = 500, bg = "white")





##### suppl figure 2: prs_threshold = 1; fpr_threshold = 5e-8 ####

prs_threshold = 1; fpr_threshold = 5e-8

for (i in 1:length(files)) {
  
  run_summary <- readRDS(files[i])[[as.character(prs_threshold)]]
  n_all_null_snps <- alz_sumstats %>% filter(str_detect(SNP, "null")) %>% nrow()
  n_null_snps_in_prs <- alz_sumstats %>% filter(str_detect(SNP, "null") & P <= prs_threshold) %>% nrow()
  n_null_snps_not_in_prs <- alz_sumstats %>% filter(str_detect(SNP, "null") & P > prs_threshold) %>% nrow()
  run_summary[["prs_extremes_fpr"]]$snp_set_number <- rep(x = c(n_all_null_snps, n_all_null_snps, 
                                                                n_null_snps_in_prs, n_null_snps_in_prs, 
                                                                n_null_snps_not_in_prs, n_null_snps_not_in_prs))
  run_summary[["prs_extremes_fpr"]]$fpn <- as.numeric(run_summary[["prs_extremes_fpr"]]$fpr) * run_summary[["prs_extremes_fpr"]]$snp_set_number
  data_list[[i]] <- run_summary[["prs_extremes_fpr"]]
  
}

fpr_dfr <- do.call("rbind", data_list)
fpr_dfr <- fpr_dfr %>%
  mutate(overlap = as.numeric(overlap),
         overlap2 = fct_reorder(paste0(overlap, "%"), overlap),
         fpr = as.numeric(fpr), 
         sig_threshold = as.numeric(sig_threshold),
         snp_set_name = case_when(snp_set_name == "all_null_snps" ~ "All null-SNPs",
                                  snp_set_name == "null_snps_in_prs" ~ "Null-SNPs used in PRS",
                                  snp_set_name == "null_snps_not_in_prs" ~ "Null-SNPs not used in PRS"))
fpr_dfr$fpr[is.na(fpr_dfr$fpr)] <- 0

n_all_null_snps <- fpr_dfr$snp_set_number[fpr_dfr$snp_set_name == "All null-SNPs"][1]
n_null_snps_in_prs <- fpr_dfr$snp_set_number[fpr_dfr$snp_set_name == "Null-SNPs used in PRS"][1]
n_null_snps_not_in_prs <- fpr_dfr$snp_set_number[fpr_dfr$snp_set_name == "Null-SNPs not used in PRS"][1]

p1 <- fpr_dfr %>%
  filter(sig_threshold == fpr_threshold,
         snp_set_name == "All null-SNPs") %>%
  ggplot(aes(x = overlap2, y = fpr)) +
  #geom_violin() +
  geom_boxplot(width=0.3, color="black", alpha=0.2) +
  geom_point(color="grey", alpha=0.4, size = 0.5) + 
  theme_classic() +
  scale_y_continuous(sec.axis = sec_axis(~ . * n_all_null_snps, name="False positive number")) +
  theme(
    legend.position="none",
    strip.text.x = element_text(size = 8),
    text=element_text(family="Charter"),
    strip.background = element_rect(
      color="#EBEBEB", fill="#EBEBEB", size=1.5, linetype="solid"),
    plot.title = element_text(hjust = 0.5, size = 10)
  ) +
  xlab("Sample overlap") +
  ylab("False positive rate") +
  ggtitle(expression("P-value threshold "<=" 1")) +
  geom_hline(yintercept = fpr_threshold, linetype="dashed", 
             color = "red", size=0.5) +
  geom_text(aes(0, 5e-8, label = 5e-8, vjust = -0.5, hjust = -0.2, colour = "red"),
            size = 3)

prs_threshold = 5e-8; fpr_threshold = 5e-8

for (i in 1:length(files)) {
  
  run_summary <- readRDS(files[i])[[as.character(prs_threshold)]]
  n_all_null_snps <- alz_sumstats %>% filter(str_detect(SNP, "null")) %>% nrow()
  n_null_snps_in_prs <- alz_sumstats %>% filter(str_detect(SNP, "null") & P <= prs_threshold) %>% nrow()
  n_null_snps_not_in_prs <- alz_sumstats %>% filter(str_detect(SNP, "null") & P > prs_threshold) %>% nrow()
  run_summary[["prs_extremes_fpr"]]$snp_set_number <- rep(x = c(n_all_null_snps, n_all_null_snps, 
                                                                n_null_snps_in_prs, n_null_snps_in_prs, 
                                                                n_null_snps_not_in_prs, n_null_snps_not_in_prs))
  run_summary[["prs_extremes_fpr"]]$fpn <- as.numeric(run_summary[["prs_extremes_fpr"]]$fpr) * run_summary[["prs_extremes_fpr"]]$snp_set_number
  data_list[[i]] <- run_summary[["prs_extremes_fpr"]]
  
}

fpr_dfr <- do.call("rbind", data_list)
fpr_dfr <- fpr_dfr %>%
  mutate(overlap = as.numeric(overlap),
         overlap2 = fct_reorder(paste0(overlap, "%"), overlap),
         fpr = as.numeric(fpr), 
         sig_threshold = as.numeric(sig_threshold),
         snp_set_name = case_when(snp_set_name == "all_null_snps" ~ "All null-SNPs",
                                  snp_set_name == "null_snps_in_prs" ~ "Null-SNPs used in PRS",
                                  snp_set_name == "null_snps_not_in_prs" ~ "Null-SNPs not used in PRS"))
fpr_dfr$fpr[is.na(fpr_dfr$fpr)] <- 0

n_all_null_snps <- max(fpr_dfr$snp_set_number[fpr_dfr$snp_set_name == "All null-SNPs"])
n_null_snps_in_prs <- max(fpr_dfr$snp_set_number[fpr_dfr$snp_set_name == "Null-SNPs used in PRS"])
n_null_snps_not_in_prs <- max(fpr_dfr$snp_set_number[fpr_dfr$snp_set_name == "Null-SNPs not used in PRS"])

p2 <- fpr_dfr %>%
  filter(sig_threshold == fpr_threshold,
         snp_set_name == "All null-SNPs") %>%
  ggplot(aes(x = overlap2, y = fpr)) +
  #geom_violin() +
  geom_boxplot(width=0.3, color="black", alpha=0.2) +
  geom_point(color="grey", alpha=0.4, size = 0.5) + 
  theme_classic() +
  scale_y_continuous(sec.axis = sec_axis(~ . * n_all_null_snps, name="False positive number", breaks = c(0, 1))) +
  theme(
    legend.position="none",
    strip.text.x = element_text(size = 8),
    text=element_text(family="Charter"),
    strip.background = element_rect(
      color="#EBEBEB", fill="#EBEBEB", size=1.5, linetype="solid"),
    plot.title = element_text(hjust = 0.5, size = 10)
  ) +
  xlab("Sample overlap") +
  ylab("False positive rate") +
  ggtitle(expression("P-value threshold "<=" 5e-8")) +
  geom_hline(yintercept=5e-8, linetype="dashed", 
             color = "red", size=0.5) +
  geom_text(aes(0, 5e-8, label = 5e-8, vjust = -0.5, hjust = -0.2, colour = "red"),
            size = 3)

figure <- ggarrange(p1, p2,
                    labels = c("A", "B"),
                    ncol = 2, font.label = list(size = 10, family="Charter"), hjust = c(-0.5, -3.5))

ggsave(plot = figure, filename = paste0("plots/prs_threshold", prs_threshold, "_and1_fpr_threshold", fpr_threshold, ".png"),
       width = 17, height = 10, units = "cm", dpi = 500, bg = "white")



##### suppl figure 3: variance of Z-scores #####

prs_threshold = 0.05

for (i in 1:length(files)) {
  
  run_summary <- readRDS(files[i])[[as.character(prs_threshold)]]
  data_list[[i]] <- run_summary[["prs_extremes_z_stat"]]
  
}

z_var_dfr <- do.call("rbind", data_list)
z_var_dfr <- z_var_dfr %>%
  mutate()

p1 <- z_var_dfr %>%
  mutate(z_stat_var = as.numeric(z_stat_var),
         overlap = as.numeric(overlap),
         overlap2 = fct_reorder(paste0(overlap, "%"), overlap),
         snp_set_name = case_when(snp_set_name == "all_null_snps" ~ "All null-SNPs",
                                  snp_set_name == "null_snps_in_prs" ~ "Null-SNPs used in PRS",
                                  snp_set_name == "null_snps_not_in_prs" ~ "Null-SNPs not used in PRS")) %>%
  ggplot(aes(x = overlap2, y = z_stat_var)) +
  #geom_violin() +
  geom_boxplot(width=0.5, color="black", alpha=0.2) +
  geom_point(color="grey", alpha=0.4, size = 0.5) + 
  facet_wrap( ~ snp_set_name, scales = "free") +
  theme_classic() +
  theme(
    strip.text.x = element_text(size = 8),
    text=element_text(family="Charter")) +
  xlab("Sample overlap") +
  ylab("Var(Z)") +
  geom_hline(yintercept=1, linetype="dashed", 
             color = "red", size=0.5)

ggsave(plot = p1, filename = paste0("plots/z_stat_var.png"),
       width = 15, height = 10, units = "cm", dpi = 500, bg = "white")
