if (exists("snakemake")) {
  # input
  fn_summary <- snakemake@input$fn_summary
  fn_phos_unphos <- snakemake@output$fn_phos_unphos
  fn_fbs <- snakemake@output$fn_fbs
  fn_wb_1e <- snakemake@output$fn_wb_1e
  fn_wb_1d <- snakemake@output$fn_wb_1d
  fn_wb_1f <- snakemake@output$fn_wb_1f
  fn_wb_1esup3 <- snakemake@output$fn_wb_1esup3
} else {
  fn_summary <- "resources/raw/20240111_denisiometry_summarized.csv"
  fn_phos_unphos <- "results/phos_unphos_ratio.png"
  fn_fbs <- "results/fbs_unphos.png"
  fn_wb_1e <- "results/wb_1e.png"
  fn_wb_1d <- "results/wb_1d.png"
  fn_wb_1f <- "results/wb_1f.png"
  fn_wb_1esup3 <- "results/wb_1e_sup3.png"
}

library(dplyr)
library(tidyr)
library(ggplot2)

ref_band <- "Actin"
dat_summary <- read.csv(fn_summary, sep = ";")

dat_summary


# Normalize the data
dat_summary_ref <- dat_summary %>%
  filter(band_type == ref_band) %>%
  group_by(experiment, lane_idx) %>%
  summarize(ref_intensity = mean(intensity_lane))

dat_summary <- dat_summary %>%
  left_join(dat_summary_ref, by = c("experiment", "lane_idx")) %>%
  mutate(norm_intensity = intensity_lane / ref_intensity) %>%
  select(-ref_intensity)


dat_summary %>%
  filter(band_type != ref_band) %>%
  ggplot(aes(x = as.factor(methionine), y = norm_intensity, fill = band_type)) +
  facet_grid(experiment ~ mthfr + FBS, scale = "free_y") +
  geom_bar(position = "dodge", stat = "identity") +
  # g(aes(group=experiment+band_type)) +
  # theme_minimal() +
  # scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


dat_summary %>%
  filter(band_type %in% c("Phos", "Non phos")) %>%
  filter(mthfr == "WT") %>%
  group_by(experiment, lane_idx, mthfr, FBS, methionine) %>%
  summarize(norm_intensity = norm_intensity / sum(norm_intensity)) %>%
  ggplot(aes(x = as.factor(methionine), y = norm_intensity)) +
  facet_grid(experiment ~ mthfr + FBS, scale = "free_y") +
  geom_bar(position = "dodge", stat = "identity") +
  # g(aes(group=experiment+band_type)) +
  # theme_minimal() +
  # scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dat_ratio <- dat_summary %>%
  pivot_wider(
    names_from = band_type, values_from = norm_intensity,
    id_cols = c("experiment", "lane_idx", "mthfr", "FBS", "methionine")
  ) %>%
  mutate(log2_phos_ratio = log2(Phos / `Non phos`))


dat_ratio %>%
  ggplot(aes(x = methionine, y = log2_phos_ratio, color = experiment, fill = FBS)) +
  geom_point() +
  geom_line() +
  # g(aes(group=experiment+band_type)) +
  # theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


dat_met_comparison <- dat_ratio %>%
  filter(is.finite(log2_phos_ratio)) %>%
  filter(methionine %in% c(0, 1000)) %>%
  pivot_wider(id_cols = c("experiment", "mthfr", "FBS"), names_from = methionine, values_from = log2_phos_ratio)

res_met <- t.test(dat_met_comparison$`1000`, dat_met_comparison$`0`, paired = TRUE)

print(res_met)

p_ratio <- dat_ratio %>%
  filter(methionine %in% c(0, 1000)) %>%
  filter(mthfr == "WT") %>%
  ggplot(aes(x = as.factor(methionine), y = log2_phos_ratio, color = FBS, shape = experiment)) +
  geom_boxplot(aes(group = methionine), width = 0.5, alpha = 0.5) +
  geom_point(size = 5) +
  geom_line(aes(group = interaction(FBS, experiment))) +
  # g(aes(group=experiment+band_type)) +
  # theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual(values = c("#33a02c", "#6a3d9a", "#e31a1c")) +
  xlab("Methionine [uM]") +
  ylab("log2(Phos / Non phos)") +
  ggtitle("Phosphorylation ratio in WT MTHFR")
p_ratio

ggsave(fn_phos_unphos, p_ratio, width = 5, height = 3)


dat_dial_comparison <- dat_summary %>%
  pivot_wider(
    names_from = c(band_type, FBS), values_from = norm_intensity,
    id_cols = c("experiment", "mthfr", "methionine")
  ) %>%
  filter(is.finite(`Non phos_no`))


t.test(dat_dial_comparison$`Non phos_no`, dat_dial_comparison$`Non phos_dialysed`, paired = TRUE)
res_fbs <- t.test(log2(dat_dial_comparison$`Non phos_dialysed`),
  log2(dat_dial_comparison$`Non phos_no`),
  paired = TRUE
)
print(res_fbs)
print(2**res_fbs$conf.int)
p_fbs <- dat_summary %>%
  filter(FBS != "normal") %>%
  filter(band_type %in% c("Non phos")) %>%
  mutate(FBS = factor(FBS, c("no", "dialysed"))) %>%
  inner_join(dat_dial_comparison, by = c("experiment", "mthfr", "methionine")) %>%
  ggplot(aes(x = FBS, y = log2(norm_intensity), color = as.factor(methionine), shape = experiment)) +
  geom_boxplot(aes(group = FBS), alpha = 0.5, width = 0.3) +
  geom_point(size = 3) +
  # black white theme
  theme_bw() +
  geom_line(aes(group = interaction(experiment, methionine))) +
  ylab("log2(Non phos / Actin)") +
  # change color guide title
  scale_color_manual(values = c("#a6cee3", "#1f78b4"), name = "Methionine [mM]") +
  ggtitle("No FBS vs dialysed FBS in WT")

p_fbs


ggsave(fn_fbs, p_fbs, width = 5, height = 3)


col_band <- c("Phos" = "#7fc97f", "Non phos" = "#beaed4", "Degraded" = "#fdc086")
p_1e <- dat_summary %>%
  filter(experiment == "1e") %>%
  filter(band_type != "Actin") %>%
  ggplot(aes(x = band_type, y = norm_intensity)) +
  facet_grid(. ~ FBS + methionine) +
  geom_point(size = 3) +
  # rotate axis
  geom_segment(aes(yend = 0)) + # change color guide title
  scale_color_discrete(name = "Methionine [uM]")

p_1e

fix_fbs_labels <- function(x) {
  levels(x) <- paste0(levels(x), " FBS")
  return(x)
}



p_1e <- dat_summary %>%
  filter(experiment == "1e") %>%
  filter(band_type != "Actin") %>%
  mutate(FBS = factor(FBS, c("no", "dialysed", "normal"))) %>%
  mutate(FBS = fix_fbs_labels(FBS)) %>%
  mutate(band_type = factor(band_type, rev(names(col_band)))) %>%
  ggplot(aes(x = as.factor(methionine), y = norm_intensity)) +
  facet_grid(. ~ FBS) +
  geom_bar(aes(fill = band_type), position = "stack", stat = "identity") +
  ylab("Actin normalized intensity [a.u.]") +
  xlab("Methionine [mM]") +
  scale_color_discrete(name = "Methionine [uM]") +
  scale_fill_manual(values = col_band, name = "Band type") +
  theme_bw() +
  ggtitle("Quantifications figure 1e")

p_1e

ggsave(fn_wb_1e, p_1e, width = 5, height = 3)


p_1d <- dat_summary %>%
  filter(experiment == "1d") %>%
  filter(band_type != "Actin") %>%
  mutate(band_type = factor(band_type, rev(names(col_band)))) %>%
  ggplot(aes(x = as.factor(methionine), y = norm_intensity)) +
  geom_bar(aes(fill = band_type), position = "stack", stat = "identity") +
  ylab("Actin normalized intensity [a.u.]") +
  xlab("Methionine [mM]") +
  scale_color_discrete(name = "Methionine [uM]") +
  scale_fill_manual(values = col_band, name = "Band type") +
  theme_bw() +
  ggtitle("Quantifications figure 1d")
# rename fill guide

p_1d

ggsave(fn_wb_1d, p_1d, width = 3, height = 3)

p_1f <- dat_summary %>%
  filter(experiment == "1f") %>%
  filter(band_type != "Actin") %>%
  mutate(mthfr = factor(mthfr, c("WT", "T34A", "EV"))) %>%
  mutate(band_type = factor(band_type, rev(names(col_band)))) %>%
  ggplot(aes(x = as.factor(methionine), y = norm_intensity)) +
  facet_grid(. ~ mthfr) +
  geom_bar(aes(fill = band_type), position = "stack", stat = "identity") +
  ylab("Actin normalized intensity [a.u.]") +
  xlab("Methionine [mM]") +
  scale_color_discrete(name = "Methionine [mM]") +
  scale_fill_manual(values = col_band, name = "Band type") +
  theme_bw() +
  ggtitle("Quantifications figure 1f")

p_1f

ggsave(fn_wb_1f, p_1f, width = 5, height = 3)

p_1esup3 <- dat_summary %>%
  filter(experiment %in% c("1e", "sup3")) %>%
  filter(band_type != "Actin") %>%
  mutate(FBS = factor(FBS, c("no", "dialysed", "normal"))) %>%
  mutate(FBS = fix_fbs_labels(FBS)) %>%
  mutate(band_type = factor(band_type, rev(names(col_band)))) %>%
  ggplot(aes(x = as.factor(methionine), y = norm_intensity)) +
  facet_grid(experiment ~ FBS) +
  geom_bar(aes(fill = band_type), position = "stack", stat = "identity") +
  ylab("Actin normalized intensity [a.u.]") +
  xlab("Methionine [mM]") +
  scale_color_discrete(name = "Methionine [uM]") +
  scale_fill_manual(values = col_band, name = "Band type") +
  theme_bw() +
  ggtitle("Quantifications figure 1e and replicate sup. figure 3")

ggsave(fn_wb_1esup3, p_1esup3, width = 5, height = 5)
