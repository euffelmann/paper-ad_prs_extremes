GWAS of polygenic risk extremes for Alzheimer’s disease suffers from
high false positive rate
================
Emil Uffelmann
31/05/2022

------------------------------------------------------------------------

## ‘Genome-wide association of polygenic risk extremes for Alzheimer’s disease in the UK Biobank’ suffers from high false positive rate

Here we show that defining a phenotype based on the extremes of a
PRS-distribution, and then running a GWAS on this newly-defined
phenotype leads to strongly inflated false positive rates. This is
because one effectively regresses the dependent variable on itself. In
other words, SNPs used to predict the dependent variable were used to
define it.

-   Original paper: <https://www.nature.com/articles/s41598-022-12391-2>

------------------------------------------------------------------------

``` r
rm(list = ls())
library(tidyverse)
library(data.table)

# load data
data_list <- list()
files <- list.files(path = "simulation_results/", pattern = glob2rx("run_summary_*.rds"), full.names = T)
```

``` r
prs_threshold = 0.05
sim_results <- NULL

for (i in 1:length(files)) {
  
  sim_results[[i]] <- readRDS(files[i])[[as.character(prs_threshold)]]
}
```

``` r
## extract fpr and fpn data
temp1 <- NULL

for (i in 1:length(sim_results)) {
  sim_results_i <- sim_results[[i]]
  alz_sumstats <- sim_results_i$gwas_alz_sumstats
  n_all_null_snps <- alz_sumstats %>% filter(str_detect(SNP, "null")) %>% nrow()
  n_null_snps_in_prs <- alz_sumstats %>% filter(str_detect(SNP, "null") & P <= prs_threshold) %>% nrow()
  n_null_snps_not_in_prs <- alz_sumstats %>% filter(str_detect(SNP, "null") & P > prs_threshold) %>% nrow()
  sim_results_i[["prs_extremes_fpr"]]$snp_set_number <- rep(x = c(n_all_null_snps, n_all_null_snps, 
                                                                  n_null_snps_in_prs, n_null_snps_in_prs, 
                                                                  n_null_snps_not_in_prs, n_null_snps_not_in_prs))
  sim_results_i[["prs_extremes_fpr"]]$fpn <- as.numeric(sim_results_i[["prs_extremes_fpr"]]$fpr) * sim_results_i[["prs_extremes_fpr"]]$snp_set_number
  temp1[[i]] <- sim_results_i[["prs_extremes_fpr"]]
}

fpr_dfr <- do.call("rbind", temp1)
fpr_dfr <- fpr_dfr %>%
  mutate(overlap = as.numeric(overlap),
         overlap2 = fct_reorder(paste0(overlap, "%"), overlap),
         fpr = as.numeric(fpr), 
         sig_threshold = as.numeric(sig_threshold),
         snp_set_name = case_when(snp_set_name == "all_null_snps" ~ "All null-SNPs",
                                  snp_set_name == "null_snps_in_prs" ~ "Null-SNPs used in PRS",
                                  snp_set_name == "null_snps_not_in_prs" ~ "Null-SNPs not used in PRS"),
         fpr = replace_na(fpr, 0))

fpr_prs_0.05_sig_gws <- fpr_dfr[fpr_dfr$sig_threshold == 5e-8,]

fpr_prs_0.05_sig_gws_table <- fpr_prs_0.05_sig_gws %>% 
  group_by(overlap, snp_set_name) %>%
  summarise(mean_fpr = mean(fpr),
            se_fpr = sd(fpr) / sqrt(100),
            max_fpr = max(fpr),
            mean_fpn = mean(fpn),
            se_fpn = sd(fpn) / sqrt(100),
            max_fpn = max(fpn),
            .groups = 'drop')

gt_tbl <- fpr_prs_0.05_sig_gws_table %>%
  select(!overlap) %>%
  gt(rowname_col = "snp_set_name") %>%
  tab_stubhead(label = "Overlap") %>%
  tab_spanner(
    label = "False-positive rate",
    columns = c(mean_fpr, se_fpr, max_fpr)) %>%
  tab_spanner(
    label = "Number of false-positives",
    columns = c(mean_fpn, se_fpn, max_fpn)) %>%
  tab_row_group(
    label = "0%",
    rows = 1:3
  ) %>%
  tab_row_group(
    label = "50%",
    rows = 4:6
  ) %>%
  tab_row_group(
    label = "sqrt(100)%",
    rows = 7:9
  )  %>%
  cols_label(
    mean_fpr = "mean",
    se_fpr = "s.e.m.",
    max_fpr = "max",
    mean_fpn = "mean",
    se_fpn = "s.e.m.",
    max_fpn = "max"
  ) %>%
  tab_options(table.font.names = "charter") %>%
  fmt_number(
    columns = c(mean_fpr, max_fpr),
    decimals = 7,
    use_seps = FALSE
  ) %>%
  fmt_number(
    columns = se_fpr,
    decimals = 8,
    use_seps = FALSE
  ) %>%
  fmt_number(
    columns = c(mean_fpn, se_fpn, max_fpn),
    decimals = 3,
    use_seps = FALSE,
    drop_trailing_zeros = T
  ) #%>%
  #gtsave(filename = "/tables/prs_threshold0.05_fpr_threshold5e-08.html", inline_css = TRUE)

gt_tbl
```

<div id="kjlxvxwnxy" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: charter;
}

#kjlxvxwnxy .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#kjlxvxwnxy .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#kjlxvxwnxy .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#kjlxvxwnxy .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#kjlxvxwnxy .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#kjlxvxwnxy .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#kjlxvxwnxy .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#kjlxvxwnxy .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#kjlxvxwnxy .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#kjlxvxwnxy .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#kjlxvxwnxy .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#kjlxvxwnxy .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
}

#kjlxvxwnxy .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#kjlxvxwnxy .gt_from_md > :first-child {
  margin-top: 0;
}

#kjlxvxwnxy .gt_from_md > :last-child {
  margin-bottom: 0;
}

#kjlxvxwnxy .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#kjlxvxwnxy .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#kjlxvxwnxy .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#kjlxvxwnxy .gt_row_group_first td {
  border-top-width: 2px;
}

#kjlxvxwnxy .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#kjlxvxwnxy .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#kjlxvxwnxy .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#kjlxvxwnxy .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#kjlxvxwnxy .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#kjlxvxwnxy .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#kjlxvxwnxy .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#kjlxvxwnxy .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#kjlxvxwnxy .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#kjlxvxwnxy .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-left: 4px;
  padding-right: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#kjlxvxwnxy .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#kjlxvxwnxy .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#kjlxvxwnxy .gt_left {
  text-align: left;
}

#kjlxvxwnxy .gt_center {
  text-align: center;
}

#kjlxvxwnxy .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#kjlxvxwnxy .gt_font_normal {
  font-weight: normal;
}

#kjlxvxwnxy .gt_font_bold {
  font-weight: bold;
}

#kjlxvxwnxy .gt_font_italic {
  font-style: italic;
}

#kjlxvxwnxy .gt_super {
  font-size: 65%;
}

#kjlxvxwnxy .gt_two_val_uncert {
  display: inline-block;
  line-height: 1em;
  text-align: right;
  font-size: 60%;
  vertical-align: -0.25em;
  margin-left: 0.1em;
}

#kjlxvxwnxy .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 75%;
  vertical-align: 0.4em;
}

#kjlxvxwnxy .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#kjlxvxwnxy .gt_slash_mark {
  font-size: 0.7em;
  line-height: 0.7em;
  vertical-align: 0.15em;
}

#kjlxvxwnxy .gt_fraction_numerator {
  font-size: 0.6em;
  line-height: 0.6em;
  vertical-align: 0.45em;
}

#kjlxvxwnxy .gt_fraction_denominator {
  font-size: 0.6em;
  line-height: 0.6em;
  vertical-align: -0.05em;
}
</style>
<table class="gt_table">
  
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="1">Overlap</th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="3">
        <span class="gt_column_spanner">False-positive rate</span>
      </th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="3">
        <span class="gt_column_spanner">Number of false-positives</span>
      </th>
    </tr>
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">mean</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">s.e.m.</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">max</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">mean</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">s.e.m.</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">max</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr class="gt_group_heading_row">
      <td colspan="7" class="gt_group_heading">sqrt(100)%</td>
    </tr>
    <tr class="gt_row_group_first"><td class="gt_row gt_right gt_stub">All null-SNPs</td>
<td class="gt_row gt_right">0.0039393</td>
<td class="gt_row gt_right">0.00001604</td>
<td class="gt_row gt_right">0.0043424</td>
<td class="gt_row gt_right">664.95</td>
<td class="gt_row gt_right">2.708</td>
<td class="gt_row gt_right">733</td></tr>
    <tr><td class="gt_row gt_right gt_stub">Null-SNPs not used in PRS</td>
<td class="gt_row gt_right">0.0000001</td>
<td class="gt_row gt_right">0.00000006</td>
<td class="gt_row gt_right">0.0000062</td>
<td class="gt_row gt_right">0.01</td>
<td class="gt_row gt_right">0.01</td>
<td class="gt_row gt_right">1</td></tr>
    <tr><td class="gt_row gt_right gt_stub">Null-SNPs used in PRS</td>
<td class="gt_row gt_right">0.0787060</td>
<td class="gt_row gt_right">0.00035577</td>
<td class="gt_row gt_right">0.0870656</td>
<td class="gt_row gt_right">664.94</td>
<td class="gt_row gt_right">2.708</td>
<td class="gt_row gt_right">733</td></tr>
    <tr class="gt_group_heading_row">
      <td colspan="7" class="gt_group_heading">50%</td>
    </tr>
    <tr class="gt_row_group_first"><td class="gt_row gt_right gt_stub">All null-SNPs</td>
<td class="gt_row gt_right">0.0032088</td>
<td class="gt_row gt_right">0.00001588</td>
<td class="gt_row gt_right">0.0035012</td>
<td class="gt_row gt_right">541.64</td>
<td class="gt_row gt_right">2.68</td>
<td class="gt_row gt_right">591</td></tr>
    <tr><td class="gt_row gt_right gt_stub">Null-SNPs not used in PRS</td>
<td class="gt_row gt_right">0.0000001</td>
<td class="gt_row gt_right">0.00000006</td>
<td class="gt_row gt_right">0.0000062</td>
<td class="gt_row gt_right">0.01</td>
<td class="gt_row gt_right">0.01</td>
<td class="gt_row gt_right">1</td></tr>
    <tr><td class="gt_row gt_right gt_stub">Null-SNPs used in PRS</td>
<td class="gt_row gt_right">0.0641144</td>
<td class="gt_row gt_right">0.00034800</td>
<td class="gt_row gt_right">0.0710935</td>
<td class="gt_row gt_right">541.63</td>
<td class="gt_row gt_right">2.68</td>
<td class="gt_row gt_right">591</td></tr>
    <tr class="gt_group_heading_row">
      <td colspan="7" class="gt_group_heading">0%</td>
    </tr>
    <tr class="gt_row_group_first"><td class="gt_row gt_right gt_stub">All null-SNPs</td>
<td class="gt_row gt_right">0.0024164</td>
<td class="gt_row gt_right">0.00001223</td>
<td class="gt_row gt_right">0.0027133</td>
<td class="gt_row gt_right">407.88</td>
<td class="gt_row gt_right">2.064</td>
<td class="gt_row gt_right">458</td></tr>
    <tr><td class="gt_row gt_right gt_stub">Null-SNPs not used in PRS</td>
<td class="gt_row gt_right">0.0000001</td>
<td class="gt_row gt_right">0.00000006</td>
<td class="gt_row gt_right">0.0000062</td>
<td class="gt_row gt_right">0.01</td>
<td class="gt_row gt_right">0.01</td>
<td class="gt_row gt_right">1</td></tr>
    <tr><td class="gt_row gt_right gt_stub">Null-SNPs used in PRS</td>
<td class="gt_row gt_right">0.0482818</td>
<td class="gt_row gt_right">0.00026984</td>
<td class="gt_row gt_right">0.0546573</td>
<td class="gt_row gt_right">407.87</td>
<td class="gt_row gt_right">2.065</td>
<td class="gt_row gt_right">458</td></tr>
  </tbody>
  
  
</table>
</div>

``` r
gwas_alz_summary_lst <- NULL
for (i in 1:length(sim_results)) {
  gwas_alz_summary_lst[[i]] <- sim_results[[i]][["gwas_alz_summary"]]
}
gwas_alz_summary <- do.call("rbind", gwas_alz_summary_lst)

summary(gwas_alz_summary)
```

    ##    prs_r2liab      ldsc_intercept     ldsc_h2l            K       
    ##  Min.   :0.03700   Min.   :0.990   Min.   :0.0800   Min.   :0.05  
    ##  1st Qu.:0.04400   1st Qu.:1.000   1st Qu.:0.0900   1st Qu.:0.05  
    ##  Median :0.04600   Median :1.000   Median :0.1000   Median :0.05  
    ##  Mean   :0.04687   Mean   :1.001   Mean   :0.0974   Mean   :0.05  
    ##  3rd Qu.:0.05000   3rd Qu.:1.000   3rd Qu.:0.1000   3rd Qu.:0.05  
    ##  Max.   :0.05600   Max.   :1.010   Max.   :0.1200   Max.   :0.05  
    ##      lambda           neff           ncase          ncontrol     
    ##  Min.   :1.000   Min.   :69688   Min.   :18339   Min.   :348432  
    ##  1st Qu.:1.010   1st Qu.:69688   1st Qu.:18339   1st Qu.:348432  
    ##  Median :1.010   Median :69688   Median :18339   Median :348432  
    ##  Mean   :1.011   Mean   :69688   Mean   :18339   Mean   :348432  
    ##  3rd Qu.:1.010   3rd Qu.:69688   3rd Qu.:18339   3rd Qu.:348432  
    ##  Max.   :1.030   Max.   :69688   Max.   :18339   Max.   :348432
