### libraries and functions ----
library(tidyverse)
library(ggforce)
library(magrittr)
library(cowplot)

pca_tidy = function(counts, ntop = 5000, sampletable = NULL) {
  # pick the top n genes by variance
  top_genes <- counts %>%
    rowwise(id) %>% # calculate the variance for each gene
    mutate(var = var(c_across(everything()))) %>%
    ungroup() %>%
    slice_max(var, n = ntop) %>% # pick only the top n
    select(-var)
  
  # perform a PCA on the data
  pca <- top_genes %>%
    select(-id) %>% 
    na.omit() %>% 
    t() %>% 
    prcomp()
  
  # convert it to a table
  pca$x %>%
    as_tibble(rownames = "sample") %>%
    rename_with(~ pca$sdev^2 %>%        # square the sd to get the variance
                  divide_by(sum(.)) %>% # divide it by the sum of variances
                  multiply_by(100) %>%  # *100 to get a percentage of variance
                  round() %>%           # round it
                  str_c(.x, " [", ., "%]"), # paste it to the PCx 
                .cols = -sample) %>% # exclude the sample column 
    {
      if (is.null(sampletable)) .
      else left_join(., sampletable, by = c("sample" = "sampleName"))
    }
}

### Load data ----
cts_nuc <- read_tsv("data/FPKM_imbibed_seeds_7D_abc_clf_noplastids.tsv", show_col_types = FALSE) %>%
  # we need to do log2(counts + 1) or the PC axes go crazy
  mutate(across(-id, ~ log2(.x + 1)))

sample_table <- tibble(sampleName = colnames(cts_nuc)) %>% 
  filter(sampleName != "id") %>% 
  mutate(condition = str_remove(sampleName, "_.$"),
         replica = str_remove(sampleName, ".*_") %>% as_factor(),
         description = case_when(condition == "Col_0" ~ "WT 0 HAI", 
                                 condition == "Col_48" ~ "WT 48 HAI", 
                                 condition == "col0" ~ "WT 7 DAG",
                                 condition == "bmi1abc" ~ "bmi1abc 7 DAG", 
                                 condition == "clfswn" ~ "clfswn 7 DAG", 
                                 .default = condition) %>% 
           factor(levels = c("WT 0 HAI", "WT 48 HAI", "WT 7 DAG", "bmi1abc 7 DAG", "clfswn 7 DAG")))


condition_colors = c("WT 0 HAI" = "#FFEDA0", "WT 48 HAI" = "#FEB24C", "WT 7 DAG" = "#F03B20", 
                     "bmi1abc 7 DAG" = "#3182BD", "clfswn 7 DAG" = "#31A354")
replica_shapes = c("1" = 21, "2" = 24, "3" = 22, "4" = 25)

### PC plot pairwise

# PC1 vs PC2 all genes
pca_tidy(cts_nuc, sampletable = sample_table, ntop = 25000) %>% 
  ggplot(aes(x = `PC1 [61%]`, y = `PC2 [17%]`)) + 
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_mark_ellipse(aes(fill = description, label = description), con.cap = 0, 
                    label.fill = "grey95", label.margin = margin(1.5,1.5,1.5,1.5,"mm")) +
  geom_point(size=3, color="black", aes(fill = description, shape = replica)) +
  #coord_fixed() +
  scale_shape_manual(values = replica_shapes) +
  scale_fill_manual("Line", values = condition_colors) +
  labs(title = "PCA of Nuclear counts log, all genes",
       subtitle = "PC1 separates WT 0 HAI from the rest
PC2 separates PRC mutants and WT 48 HAI from other stages",
       caption = "Percentage of variance explained by each PC in square brackets") +
  guides(fill = guide_legend(override.aes = list(shape = 22, color = NA)), 
         shape = guide_legend(override.aes = list(fill = "darkseagreen2", size = 3))) +
  theme_cowplot() +
  theme(axis.line=element_blank(), panel.border = element_rect(colour = "black", linewidth = 1, fill = NA))

dir.create("figures", showWarnings = FALSE)
ggsave2(filename = "figures/PCA_1-vs-2_all-genes.pdf", width = 7, height = 7)




