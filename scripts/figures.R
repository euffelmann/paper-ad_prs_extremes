

##### set-up data #####

rm(list = ls())
library(tidyverse)
library(data.table)
library(ggpubr)
library(gt)

data_list <- list()
alz_fpr_fpn_list <- list() # fpn = False positive number
files <- list.files(path = "simulation_results", pattern = glob2rx("run_summary_*.rds"), full.names = T)

##### main figure 1: prs_threshold = 0.05; fpr_threshold = 5e-8 ####

## false positive rate

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
         snp_set_name = case_when(snp_set_name == "all_null_snps" ~ "all null-SNPs",
                                  snp_set_name == "null_snps_in_prs" ~ "null-SNPs used in PRS",
                                  snp_set_name == "null_snps_not_in_prs" ~ "null-SNPs not used in PRS"))
fpr_dfr$fpr[is.na(fpr_dfr$fpr)] <- 0

n_all_null_snps <- fpr_dfr$snp_set_number[fpr_dfr$snp_set_name == "all null-SNPs"][1]
n_null_snps_in_prs <- fpr_dfr$snp_set_number[fpr_dfr$snp_set_name == "null-SNPs used in PRS"][1]
n_null_snps_not_in_prs <- fpr_dfr$snp_set_number[fpr_dfr$snp_set_name == "null-SNPs not used in PRS"][1]

p1 <- fpr_dfr %>%
  filter(sig_threshold == fpr_threshold,
         snp_set_name == "all null-SNPs") %>%
  group_by(overlap2, snp_set_name) %>%
  summarise(mean_fpr = mean(fpr),
            se_fpr = sd(fpr) / sqrt(100),
            max_fpr = max(fpr),
            mean_fpn = mean(fpn),
            se_fpn = sd(fpn) / sqrt(100),
            max_fpn = max(fpn),
            .groups = 'drop') %>%
  ggplot(aes(x = overlap2, y = mean_fpr)) +
  geom_point(size = 1, colour = "black") + 
  geom_errorbar(aes(ymin = mean_fpr - 3.291 * se_fpr, ymax = mean_fpr + 3.291 * se_fpr), width=.3) +
  theme_classic() +
  scale_y_continuous(sec.axis = sec_axis(~ . * n_all_null_snps, name="")) +
  theme(
    legend.position="none",
    strip.text.x = element_text(size = 8),
    text=element_text(family="Helvetica"),
    strip.background = element_rect(
      color="#EBEBEB", fill="#EBEBEB", size=1.5, linetype="solid"),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.title=element_text(size=10)
  ) +
  xlab("sample overlap") +
  ylab("false positive rate per SNP") +
  ggtitle("all null-SNPs") +
  geom_hline(yintercept=5e-8, linetype="dashed", 
             color = "red", size=0.5) +
  geom_text(aes(1, 5e-8, label = 5e-8, colour = "red"),
            size = 3, vjust = -0.5, hjust = 0.5)

p2 <- fpr_dfr %>%
  filter(sig_threshold == fpr_threshold,
         snp_set_name == "null-SNPs used in PRS") %>%
  group_by(overlap2, snp_set_name) %>%
  summarise(mean_fpr = mean(fpr),
            se_fpr = sd(fpr) / sqrt(100),
            max_fpr = max(fpr),
            mean_fpn = mean(fpn),
            se_fpn = sd(fpn) / sqrt(100),
            max_fpn = max(fpn),
            .groups = 'drop') %>%
  ggplot(aes(x = overlap2, y = mean_fpr)) +
  geom_point(size = 1, colour = "black") + 
  geom_errorbar(aes(ymin = mean_fpr - 3.291 * se_fpr, ymax = mean_fpr + 3.291 * se_fpr), width=.3) +
  theme_classic() +
  scale_y_continuous(sec.axis = sec_axis(~ . * n_null_snps_in_prs, name="")) +
  theme(
    legend.position="none",
    strip.text.x = element_text(size = 8),
    text=element_text(family="Helvetica"),
    strip.background = element_rect(
      color="#EBEBEB", fill="#EBEBEB", size=1.5, linetype="solid"),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.title=element_text(size=10)
  ) +
  xlab("sample overlap") +
  ylab("") +
  ggtitle("null-SNPs used in PRS") +
  geom_hline(yintercept=5e-8, linetype="dashed", 
             color = "red", size=0.5)

p3 <- fpr_dfr %>%
  filter(sig_threshold == fpr_threshold,
         snp_set_name == "null-SNPs not used in PRS") %>%
  group_by(overlap2, snp_set_name) %>%
  summarise(mean_fpr = mean(fpr),
            se_fpr = sd(fpr) / sqrt(100),
            max_fpr = max(fpr),
            mean_fpn = mean(fpn),
            se_fpn = sd(fpn) / sqrt(100),
            max_fpn = max(fpn),
            .groups = 'drop') %>%
  ggplot(aes(x = overlap2, y = mean_fpr)) +
  geom_point(size = 1, colour = "black") + 
  geom_errorbar(aes(ymin = 0, ymax = mean_fpr + 3.291 * se_fpr), width=.3) +
  scale_y_continuous(sec.axis = sec_axis(~ . * n_null_snps_not_in_prs, name="# false positives per study", breaks = c(0, 1)),
                     limits = c(0, 6.2e-6)) +
  theme_classic() +
  theme(
    legend.position="none",
    strip.text.x = element_text(size = 8),
    text=element_text(family="Helvetica"),
    strip.background = element_rect(
      color="#EBEBEB", fill="#EBEBEB", size=1.5, linetype="solid"),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.title=element_text(size=10)
  ) +
  xlab("sample overlap") +
  ylab("") +
  ggtitle("null-SNPs not used in PRS") +
  geom_hline(yintercept=5e-8, linetype="dashed", 
             color = "red", size=0.5)

figure <- ggarrange(p1, p2, p3,
                    labels = c("a", "b", "c"),
                    ncol = 3, nrow = 1, font.label = list(size = 10, family="Helvetica", face = "bold"),
                    align = "v")

ggsave(plot = figure, filename = paste0("plots/fpr_fpn_prs_threshold", prs_threshold, "_fpr_threshold", fpr_threshold, ".png"),
       width = 20, height = 6, units = "cm", dpi = 500, bg = "white")


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
         snp_set_name = case_when(snp_set_name == "all_null_snps" ~ "all null-SNPs",
                                  snp_set_name == "null_snps_in_prs" ~ "null-SNPs used in PRS",
                                  snp_set_name == "null_snps_not_in_prs" ~ "null-SNPs not used in PRS"))
fpr_dfr$fpr[is.na(fpr_dfr$fpr)] <- 0

n_all_null_snps <- fpr_dfr$snp_set_number[fpr_dfr$snp_set_name == "all null-SNPs"][1]
n_null_snps_in_prs <- fpr_dfr$snp_set_number[fpr_dfr$snp_set_name == "null-SNPs used in PRS"][1]
n_null_snps_not_in_prs <- fpr_dfr$snp_set_number[fpr_dfr$snp_set_name == "null-SNPs not used in PRS"][1]

p1 <- fpr_dfr %>%
  filter(sig_threshold == fpr_threshold,
         snp_set_name == "all null-SNPs") %>%
  group_by(overlap2, snp_set_name) %>%
  summarise(mean_fpr = mean(fpr),
            se_fpr = sd(fpr) / sqrt(100),
            max_fpr = max(fpr),
            mean_fpn = mean(fpn),
            se_fpn = sd(fpn) / sqrt(100),
            max_fpn = max(fpn),
            .groups = 'drop') %>%
  ggplot(aes(x = overlap2, y = mean_fpr)) +
  geom_point(size = 1, colour = "black") + 
  geom_errorbar(aes(ymin = mean_fpr - 3.291 * se_fpr, ymax = mean_fpr + 3.291 * se_fpr), width=.3) +
  theme_classic() +
  scale_y_continuous(sec.axis = sec_axis(~ . * n_all_null_snps, name="")) +
  theme(
    legend.position="none",
    strip.text.x = element_text(size = 8),
    text=element_text(family="Helvetica"),
    strip.background = element_rect(
      color="#EBEBEB", fill="#EBEBEB", size=1.5, linetype="solid"),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.title=element_text(size=10)
  ) +
  xlab("sample overlap") +
  ylab("false positive rate per SNP") +
  ggtitle("all null-SNPs") +
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "red", size=0.5)

p2 <- fpr_dfr %>%
  filter(sig_threshold == fpr_threshold,
         snp_set_name == "null-SNPs used in PRS") %>%
  group_by(overlap2, snp_set_name) %>%
  summarise(mean_fpr = mean(fpr),
            se_fpr = sd(fpr) / sqrt(100),
            max_fpr = max(fpr),
            mean_fpn = mean(fpn),
            se_fpn = sd(fpn) / sqrt(100),
            max_fpn = max(fpn),
            .groups = 'drop') %>%
  ggplot(aes(x = overlap2, y = mean_fpr)) +
  geom_point(size = 1, colour = "black") + 
  geom_errorbar(aes(ymin = mean_fpr - 3.291 * se_fpr, ymax = mean_fpr + 3.291 * se_fpr), width=.3) +
  theme_classic() +
  scale_y_continuous(sec.axis = sec_axis(~ . * n_null_snps_in_prs, name=""),
                     breaks = c(0.05, 0.25, 0.5, 0.75, 1)) +
  theme(
    legend.position="none",
    strip.text.x = element_text(size = 8),
    text=element_text(family="Helvetica"),
    strip.background = element_rect(
      color="#EBEBEB", fill="#EBEBEB", size=1.5, linetype="solid"),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.title=element_text(size=10)
  ) +
  xlab("sample overlap") +
  ylab("") +
  ggtitle("null-SNPs used in PRS") +
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "red", size=0.5)

p3 <- fpr_dfr %>%
  filter(sig_threshold == fpr_threshold,
         snp_set_name == "null-SNPs not used in PRS") %>%
  group_by(overlap2, snp_set_name) %>%
  summarise(mean_fpr = mean(fpr),
            se_fpr = sd(fpr) / sqrt(100),
            max_fpr = max(fpr),
            mean_fpn = mean(fpn),
            se_fpn = sd(fpn) / sqrt(100),
            max_fpn = max(fpn),
            .groups = 'drop') %>%
  ggplot(aes(x = overlap2, y = mean_fpr)) +
  geom_point(size = 1, colour = "black") + 
  geom_errorbar(aes(ymin = mean_fpr - 3.291 * se_fpr, ymax = mean_fpr + 3.291 * se_fpr), width=.3) +
  theme_classic() +
  scale_y_continuous(sec.axis = sec_axis(~ . * n_null_snps_not_in_prs, name="# false positives per study")) +
  theme(
    legend.position="none",
    strip.text.x = element_text(size = 8),
    text=element_text(family="Helvetica"),
    strip.background = element_rect(
      color="#EBEBEB", fill="#EBEBEB", size=1.5, linetype="solid"),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.title=element_text(size=10)
  ) +
  xlab("sample overlap") +
  ylab("") +
  ggtitle("null-SNPs not used in PRS") +
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "red", size=0.5)


figure <- ggarrange(p1, p2, p3,
                    labels = c("a", "b", "c"),
                    ncol = 3, nrow = 1, font.label = list(size = 10, family="Helvetica", face = "bold"), hjust = c(-0.5, -3.5))

ggsave(plot = figure, filename = paste0("plots/prs_threshold", prs_threshold, "_fpr_threshold", fpr_threshold, ".png"),
       width = 20, height = 7, units = "cm", dpi = 500, bg = "white")


##### suppl figure 2: prs_threshold = 1 or 5e-8 #####

## false positive rate
prs_threshold = 1; fpr_threshold = 5e-8

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
  
}

fpr_dfr <- do.call("rbind", data_list)
fpr_dfr_prs1 <- fpr_dfr %>%
  mutate(overlap = as.numeric(overlap),
         overlap2 = fct_reorder(paste0(overlap, "%"), overlap),
         fpr = as.numeric(fpr), 
         sig_threshold = as.numeric(sig_threshold),
         snp_set_name = case_when(snp_set_name == "all_null_snps" ~ "all null-SNPs",
                                  snp_set_name == "null_snps_in_prs" ~ "null-SNPs used in PRS",
                                  snp_set_name == "null_snps_not_in_prs" ~ "null-SNPs not used in PRS"))
fpr_dfr_prs1$fpr[is.na(fpr_dfr_prs1$fpr)] <- 0

n_all_null_snps <- fpr_dfr_prs1$snp_set_number[fpr_dfr_prs1$snp_set_name == "all null-SNPs"][1]
n_null_snps_in_prs <- fpr_dfr_prs1$snp_set_number[fpr_dfr_prs1$snp_set_name == "null-SNPs used in PRS"][1]
n_null_snps_not_in_prs <- fpr_dfr_prs1$snp_set_number[fpr_dfr_prs1$snp_set_name == "null-SNPs not used in PRS"][1]

p1 <- fpr_dfr_prs1 %>%
  filter(sig_threshold == fpr_threshold,
         snp_set_name == "all null-SNPs") %>%
  group_by(overlap2, snp_set_name) %>%
  summarise(mean_fpr = mean(fpr),
            se_fpr = sd(fpr) / sqrt(100),
            max_fpr = max(fpr),
            mean_fpn = mean(fpn),
            se_fpn = sd(fpn) / sqrt(100),
            max_fpn = max(fpn),
            .groups = 'drop') %>%
  ggplot(aes(x = overlap2, y = mean_fpr)) +
  geom_point(size = 1, colour = "black") + 
  geom_errorbar(aes(ymin = mean_fpr - 3.291 * se_fpr, ymax = mean_fpr + 3.291 * se_fpr), width=.3) +
  theme_classic() +
  scale_y_continuous(sec.axis = sec_axis(~ . * n_all_null_snps, name="")) +
  theme(
    legend.position="none",
    strip.text.x = element_text(size = 8),
    text=element_text(family="Helvetica"),
    strip.background = element_rect(
      color="#EBEBEB", fill="#EBEBEB", size=1.5, linetype="solid"),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.title=element_text(size=10)
  ) +
  xlab("sample overlap") +
  ylab("false positive rate per SNP") +
  ggtitle("p-value threshold: 1") +
  geom_hline(yintercept=5e-8, linetype="dashed", 
             color = "red", size=0.5) +
  geom_text(aes(1, 5e-8, label = 5e-8, colour = "red"),
            size = 3, vjust = -0.5, hjust = 0.5)

prs_threshold = 5e-8; fpr_threshold = 5e-8

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
  
}

fpr_dfr <- do.call("rbind", data_list)
fpr_dfr_prsbonf <- fpr_dfr %>%
  mutate(overlap = as.numeric(overlap),
         overlap2 = fct_reorder(paste0(overlap, "%"), overlap),
         fpr = as.numeric(fpr), 
         sig_threshold = as.numeric(sig_threshold),
         snp_set_name = case_when(snp_set_name == "all_null_snps" ~ "all null-SNPs",
                                  snp_set_name == "null_snps_in_prs" ~ "null-SNPs used in PRS",
                                  snp_set_name == "null_snps_not_in_prs" ~ "null-SNPs not used in PRS"))
fpr_dfr_prsbonf$fpr[is.na(fpr_dfr_prsbonf$fpr)] <- 0

n_all_null_snps <- fpr_dfr_prsbonf$snp_set_number[fpr_dfr_prsbonf$snp_set_name == "all null-SNPs"][1]
n_null_snps_in_prs <- fpr_dfr_prsbonf$snp_set_number[fpr_dfr_prsbonf$snp_set_name == "null-SNPs used in PRS"][1]
n_null_snps_not_in_prs <- fpr_dfr_prsbonf$snp_set_number[fpr_dfr_prsbonf$snp_set_name == "null-SNPs not used in PRS"][1]

p2 <- fpr_dfr_prsbonf %>%
  filter(sig_threshold == fpr_threshold,
         snp_set_name == "all null-SNPs") %>%
  group_by(overlap2, snp_set_name) %>%
  summarise(mean_fpr = mean(fpr),
            se_fpr = sd(fpr) / sqrt(100),
            max_fpr = max(fpr),
            mean_fpn = mean(fpn),
            se_fpn = sd(fpn) / sqrt(100),
            max_fpn = max(fpn),
            .groups = 'drop') %>%
  ggplot(aes(x = overlap2, y = mean_fpr)) +
  geom_point(size = 1, colour = "black") + 
  geom_errorbar(aes(ymin = 0, ymax = mean_fpr + 3.291 * se_fpr), width=.3) +
  scale_y_continuous(sec.axis = sec_axis(~ . * n_all_null_snps, name="# false positives per study", breaks = c(0, 1)),
                     limits = c(0, 6.2e-6)) +
  theme_classic() +
  theme(
    legend.position="none",
    strip.text.x = element_text(size = 8),
    text=element_text(family="Helvetica"),
    strip.background = element_rect(
      color="#EBEBEB", fill="#EBEBEB", size=1.5, linetype="solid"),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.title=element_text(size=10)
  ) +
  xlab("sample overlap") +
  ylab("") +
  ggtitle("p-value threshold: 5e-8") +
  geom_hline(yintercept=5e-8, linetype="dashed", 
             color = "red", size=0.5) +
  geom_text(aes(1, 5e-8, label = 5e-8, colour = "red"),
            size = 3, vjust = -0.5, hjust = 3)


figure <- ggarrange(p1, p2,
                    labels = c("a", "b"),
                    ncol = 2, nrow = 1, font.label = list(size = 10, family="Helvetica", face = "bold"), hjust = c(-0.5, -3.5), align = "v")

ggsave(plot = figure, filename = paste0("plots/prs_threshold1_and_5e-8.png"),
       width = 15, height = 6, units = "cm", dpi = 500, bg = "white")



##### supplementary tables: set-up #####

print_w <- function(vector, format = "SE"){
  vector <- as.numeric(vector)
  mean <- formatC(mean(vector),format="e",digits=2)
  I95 <- paste(formatC(quantile(vector,0.025,na.rm=TRUE),format="e",digits=2)," to ",formatC(quantile(vector,0.975,na.rm=TRUE),format="e",digits=2),sep="")
  SE <- formatC(sd(vector)/sqrt(length((vector))),format="e",digits=2)
  if(format=="SE"){out <- paste(mean," (",SE,")",sep="")}
  if(format=="I95"){out <- paste(mean," (",I95,")",sep="")}
  return(out)
}

prs_threshold = 0.05
sim_results <- NULL
fpr_temp <- NULL

#### supplementary tables 1: variance of test statistics, false positive rate, number of false positives ####

for (i in 1:length(files)) {
  
  sim_results[[i]] <- readRDS(files[i])[[as.character(prs_threshold)]]
  alz_sumstats <- sim_results[[i]]$gwas_alz_sumstats
  n_all_null_snps <- alz_sumstats %>% filter(str_detect(SNP, "null")) %>% nrow()
  n_null_snps_in_prs <- alz_sumstats %>% filter(str_detect(SNP, "null") & P <= prs_threshold) %>% nrow()
  n_null_snps_not_in_prs <- alz_sumstats %>% filter(str_detect(SNP, "null") & P > prs_threshold) %>% nrow()
  sim_results[[i]][["prs_extremes_fpr"]]$snp_set_number <- rep(x = c(n_all_null_snps, n_all_null_snps, 
                                                                  n_null_snps_in_prs, n_null_snps_in_prs, 
                                                                  n_null_snps_not_in_prs, n_null_snps_not_in_prs))
  sim_results[[i]][["prs_extremes_fpr"]]$fpn <- as.numeric(sim_results[[i]][["prs_extremes_fpr"]]$fpr) * sim_results[[i]][["prs_extremes_fpr"]]$snp_set_number
  fpr_temp[[i]] <- sim_results[[i]][["prs_extremes_fpr"]]
}

fpr_dfr <- as.data.frame(do.call("rbind", fpr_temp))
fpr_dfr <- fpr_dfr %>%
  filter(sig_threshold == 5e-8) %>%
  select(!sig_threshold)

fpr_summary <- fpr_dfr %>%
  mutate(overlap = as.numeric(overlap),
         overlap2 = fct_reorder(paste0(overlap, "%"), overlap),
         fpr = as.numeric(fpr), 
         snp_set_name = case_when(snp_set_name == "all_null_snps" ~ "all null-SNPs",
                                  snp_set_name == "null_snps_in_prs" ~ "null-SNPs used in PRS",
                                  snp_set_name == "null_snps_not_in_prs" ~ "null-SNPs not used in PRS"),
         fpr = replace_na(fpr, 0)) %>%
  group_by(overlap2, snp_set_name) %>%
  summarise(fpr = print_w(fpr),
            max_fpr = formatC(max(fpr), format="e",digits=2),
            fpn = print_w(fpn),
            max_fpn = formatC(max(fpn), format="e",digits=2),
            .groups = 'drop')

gt_tbl <- fpr_summary %>%
  select(!overlap2) %>%
  gt(rowname_col = "snp_set_name") %>%
  tab_stubhead(label = "overlap") %>%
  tab_spanner(
    label = "false positive rate per SNP",
    columns = c(fpr, max_fpr)) %>%
  tab_spanner(
    label = "# false positives per study",
    columns = c(fpn, max_fpn)) %>%
  tab_row_group(
    label = "0%",
    rows = 1:3
  ) %>%
  tab_row_group(
    label = "50%",
    rows = 4:6
  ) %>%
  tab_row_group(
    label = "100%",
    rows = 7:9
  )  %>%
  cols_label(
    fpr = "mean (s.e.m.)",
    max_fpr = "max",
    fpn = "mean (s.e.m.)",
    max_fpn = "max"
  ) %>%
  tab_options(table.font.names = "Helvetica") %>%
  tab_options(
    table.border.top.color = "white",
    heading.title.font.size = px(16),
    column_labels.border.top.width = 3,
    column_labels.border.top.color = "black",
    column_labels.border.bottom.width = 3,
    column_labels.border.bottom.color = "black",
    table_body.border.bottom.color = "black",
    table.border.bottom.color = "white",
    table.width = pct(100),
    table.background.color = "white"
  ) %>%
  cols_align(align="center") %>%
  tab_style(
    style = list(
      cell_borders(
        sides = c("top", "bottom"),
        color = "white",
        weight = px(1)
      ),
      cell_text(
        align="center"
      ),
      cell_fill(color = "white", alpha = NULL)
    ),
    locations = cells_body(
      columns = everything(),
      rows = everything()
    )
  ) %>%
  opt_align_table_header(align = "left") %>%
gtsave(filename = "tables/table1.html", inline_css = TRUE)


