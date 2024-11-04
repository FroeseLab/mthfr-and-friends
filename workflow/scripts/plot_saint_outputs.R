if (exists("snakemake")) {
  # input
  fn_saint_annotated <- snakemake@input$saint_annotated
  # output
  fn_plot_sp_fca <- snakemake@output$plot_sp_fca
  fn_plot_sp_fca_kinase <- snakemake@output$plot_sp_fca_kinase
  fn_plot_sp_fca_string <- snakemake@output$plot_sp_fca_string
  fn_plot_sp_fca_crapome_gradient <- snakemake@output$plot_sp_fca_crapome_gradient
  fn_plot_sp_fca_crapome <- snakemake@output$plot_sp_fca_crapome
  fn_plot_heatmap <- snakemake@output$plot_heatmap
} else {
  fn_saint_annotated <- "data/output/saint_output_merged_annotated_MTHFR38to656.rds"
  # output
  fn_plot_sp_fca_check <- "results/tmp/plot_sp_fca_check.png"
  fn_plot_sp_fca_kinase <- "results/tmp/plot_sp_fca_kinase.png"
  fn_plot_sp_fca_string <- "results/tmp/plot_sp_fca_string.png"
  fn_plot_sp_fca_crapome_gradient <- "results/tmp/plot_sp_fca_crapome_gradient.png"
  fn_plot_sp_fca_crapome <- "results/tmp/plot_sp_fca_crapome.png"
  fn_plot_heatmap <- "results/tmp/plot_heatmap.png"
}

# Read the input file
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(grid)
library(gridExtra)
library("ggrepel")
library(STRINGdb)
library(tidyverse)
library(ggpubr)


data_SAINT_merged_modified <- readRDS(fn_saint_annotated)

# Plotting ----------------------------------------------------------------------------

## Main Plot ----

# Temporarely to make crapome points go in the back
plot_main <- data_SAINT_merged_modified %>%
  arrange(is_not_crapome) %>% # Change this to your desired ordering logic


  ggplot(
    data = ,
    aes(
      x = SP,
      y = FCA,
      label = PreyGene,
      fill = is_not_crapome,
      color = merged_activity
    )
  ) +
  geom_point(
    size = 2.5,
    shape = 21,
    stroke = 0.3
  ) +

  # Limits for high confidence interactions
  geom_hline(
    yintercept = 2,
    color = "darkcyan",
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = 0.77,
    y = 2.1,
    label = "FCA = 2",
    color = "darkcyan"
  ) +
  geom_vline(
    xintercept = 0.95,
    color = "darkcyan",
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = 0.945,
    y = 8,
    angle = 90,
    label = "SP = 0.95",
    color = "darkcyan"
  ) +
  scale_color_manual(
    values = c(
      "FALSE" = "white",
      "kinase" = "white",
      "kinase regulator activity" = "white"
    ),
    labels = c("FALSE" = " ")
  ) +
  scale_fill_manual(
    values = c("TRUE" = "black", "FALSE" = "#B38A5E"),
    labels = c("TRUE" = "Not crapome", "FALSE" = "Crapome")
  ) +
  coord_cartesian(xlim = c(0.75, 1), ylim = c(1, 10)) +
  geom_text_repel(
    data = filter(
      data_SAINT_merged_modified,
      FCA > 2 &
        (is_not_crapome) &
        SP > 0.95 & PreyGene != "MTHFR"
    ),
    aes(label = PreyGene),
    color = "black",
    max.overlaps = 140,
    size = 5,
    #
    nudge_y = 0.07,
    nudge_x = -0.12,
    #
    box.padding = 0.4,
    segment.colour = "darkgray",
  ) +
  geom_text_repel(
    data = filter(data_SAINT_merged_modified, PreyGene == "MTHFR"),
    aes(label = PreyGene),
    color = "black",
    max.overlaps = 140,
    size = 5,
    #
    nudge_y = 0.05,
    nudge_x = -0.01,
    #
    box.padding = 0.4,
    segment.colour = "darkgray",
  ) +
  facet_grid(~Bait_name, labeller = label_parsed) + # label_parsed allows for the subscript to show in graph

  scale_y_log10() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 20),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
  ) +
  labs(fill = NULL, color = NULL) + # Remove legend titles
  guides(
    fill = guide_legend(title = NULL), # Remove the fill legend title
    color = "none"
  ) # Remove the color legend title

plot_main

# Have to separate the legend as this affects the size of the images, and the plots dont align with eachother
legend <- as_ggplot(cowplot::get_legend(plot_main))
grid <- grid.arrange(
  plot_main + theme(legend.position = "none"),
  legend,
  nrow = 1,
  widths = c(8, 1)
)


ggsave(
  fn_plot_sp_fca,
  grid,
  height = 10,
  width = 40,
  units = "cm",
  dpi = 600,
)



## Kinases ----
# The protein CAD seems to be incorrectly labled with the kinase GO term and will therfore not be highlighted in the graph
data_SAINT_merged_modified[which(data_SAINT_merged_modified$PreyGene == "CAD"), "kinase"] <-
  FALSE
data_SAINT_merged_modified[which(data_SAINT_merged_modified$PreyGene == "CAD"), "merged_activity"] <-
  FALSE

plot_kinases <- ggplot(
  data = filter(data_SAINT_merged_modified, merged_activity != FALSE),
  aes(
    x = SP,
    y = FCA,
    label = PreyGene,
    fill = is_not_crapome,
    color = merged_activity,
    #     alpha = merged_activity
  )
) +
  geom_point(data = data_SAINT_merged_modified, aes(x = SP, y = FCA), color = "#E5E7E9") +
  geom_point(size = 2.5, shape = 21) +

  # Limits for high confidence interactions
  geom_hline(
    yintercept = 2,
    color = "darkcyan",
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = 0.77,
    y = 2.1,
    label = "FCA = 2",
    color = "darkcyan"
  ) +
  geom_vline(
    xintercept = 0.95,
    color = "darkcyan",
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = 0.945,
    y = 8,
    angle = 90,
    label = "SP = 0.95",
    color = "darkcyan"
  ) +
  scale_color_manual(
    values = c(
      "FALSE" = "white",
      "kinase" = "#D41159",
      "kinase regulator" = "#1A85FF"
    ),
    labels = c("FALSE" = " ")
  ) +
  scale_fill_manual(
    values = c("TRUE" = "black", "FALSE" = "#B38A5E"),
    labels = c("TRUE" = "Not crapome", "FALSE" = "Crapome")
  ) +
  # scale_alpha_manual ( values =c("FALSE" = 0.1), breaks = waiver(), na.value = NA)+
  coord_cartesian(xlim = c(0.75, 1), ylim = c(1, 10)) +
  geom_text_repel(
    data = filter(
      data_SAINT_merged_modified,
      FCA > 2 & SP > 0.95 & merged_activity != FALSE
    ),
    aes(label = PreyGene),
    color = "black",
    max.overlaps = 140,
    size = 5,
    nudge_y = 0.05,
    nudge_x = -0.05,
    box.padding = 0.4,
    segment.colour = "darkgray",
  ) +
  facet_grid(~Bait_name, labeller = label_parsed) + # label_parsed allows for the subscript to show in graph


  scale_y_log10() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 20),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15)
  ) +
  labs(fill = NULL, color = NULL) + # Remove legend titles
  guides(
    fill = guide_legend(title = NULL), # Remove the fill legend title
    color = guide_legend(title = NULL)
  ) # Remove the color legend title

plot_kinases

# Have to separate the legend as this affects the size of the images, and the plots dont align with eachother
legend <- as_ggplot(cowplot::get_legend(plot_kinases))

grid <- grid.arrange(
  plot_kinases + theme(legend.position = "none"),
  legend,
  nrow = 1,
  widths = c(8, 1)
)

ggsave(
  fn_plot_sp_fca_kinase,
  grid,
  height = 10,
  width = 40,
  units = "cm",
  dpi = 600,
)



## STRING  ----

# Temporarely
plot_STRING <- data_SAINT_merged_modified %>%
  arrange(STRING_score) %>% # Change this to your desired ordering logic
  arrange(merged_activity) %>%
  filter(STRING_score != FALSE) %>%
  ggplot(
    data = ,
    aes(
      x = SP,
      y = FCA,
      label = PreyGene,
      fill = STRING_score,
      color = merged_activity
    )
  ) +
  geom_point(data = data_SAINT_merged_modified, aes(x = SP, y = FCA), color = "#E5E7E9") +
  geom_point(
    size = 2.5,
    shape = 21,
    stroke = 0.5
  ) +

  # Limits for high confidence interactions
  geom_hline(
    yintercept = 2,
    color = "darkcyan",
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = 0.77,
    y = 2.1,
    label = "FCA = 2",
    color = "darkcyan"
  ) +
  geom_vline(
    xintercept = 0.95,
    color = "darkcyan",
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = 0.945,
    y = 8,
    angle = 90,
    label = "SP = 0.95",
    color = "darkcyan"
  ) +
  scale_color_manual(
    values = c(
      "FALSE" = "black",
      "kinase" = "black",
      "kinase regulator activity" = "black"
    ),
    labels = c("FALSE" = " ")
  ) +
  scale_fill_gradient(
    low = "gray",
    high = "purple",
    na.value = "black"
  ) +
  coord_cartesian(xlim = c(0.75, 1), ylim = c(1, 10)) +
  geom_text_repel(
    data = filter(
      data_SAINT_merged_modified,
      FCA > 2 &
        SP > 0.95 & STRING_score != FALSE
    ),
    aes(label = PreyGene),
    color = "black",
    max.overlaps = 140,
    size = 5,
    nudge_y = 0.05,
    nudge_x = -0.05,
    box.padding = 0.4,
    segment.colour = "darkgray",
  ) +
  facet_grid(~Bait_name, labeller = label_parsed) + # label_parsed allows for the subscript to show in graph



  scale_y_log10() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 20),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15)
  ) +
  labs(fill = NULL, color = NULL) + # Remove legend titles
  guides(fill = guide_legend(title = "STRING score"), color = "none") # Remove the color legend title

plot_STRING

# Have to separate the legend as this affects the size of the images, and the plots dont align with eachother
legend <- as_ggplot(cowplot::get_legend(plot_STRING))
grid <- grid.arrange(
  plot_STRING + theme(legend.position = "none"),
  legend,
  nrow = 1,
  widths = c(8, 1)
)


ggsave(
  fn_plot_sp_fca_string,
  grid,
  height = 10,
  width = 40,
  units = "cm",
  dpi = 600,
)

# grid arrange
# Arrange the plots using gridExtra
grid_plots <- grid.arrange(plot_main, plot_kinases, plot_STRING, ncol = 1)
grid_plots

## Crapome [NOT IN ARTICLE] ----


ggplot(
  data = data_SAINT_merged_modified,
  aes(
    x = SP,
    y = FCA,
    label = PreyGene,
    fill = crapome_frequency,
    color = merged_activity
  )
) +
  geom_hline(
    yintercept = 2,
    color = "darkcyan",
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = 0.77,
    y = 2.1,
    label = "FCA = 2",
    color = "darkcyan"
  ) +
  geom_vline(
    xintercept = 0.95,
    color = "darkcyan",
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = 0.94,
    y = 9,
    angle = 90,
    label = "SP = 0.95",
    color = "darkcyan"
  ) +
  facet_grid(~Bait) +
  geom_point(size = 2, shape = 21) +
  scale_color_manual(
    values = c(
      "FALSE" = "white",
      "kinase" = "red",
      "kinase regulator activity" = "green"
    ),
    labels = c("FALSE" = " ")
  ) +
  scale_fill_gradient(
    low = "green",
    high = "red",
    na.value = "black"
  ) +
  coord_cartesian(xlim = c(0.75, 1), ylim = c(1, 10)) +
  # geom_text_repel(
  #   data = filter(data_SAINT_merged_modified, kinase==TRUE),
  #   aes(label = PreyGene),
  #   color = 'black',
  #   max.overlaps = 140,
  #   size = 3,
  #   nudge_y = 0.05,
  #   nudge_x = -0.1,
  #   box.padding = 1,

  # ) +
  scale_y_log10() +
  theme_bw() +
  labs(fill = NULL, color = NULL) + # Remove legend titles
  guides(
    fill = guide_legend(title = NULL), # Remove the fill legend title
    color = guide_legend(title = NULL)
  ) # Remove the color legend title

ggsave(
  fn_plot_sp_fca_crapome_gradient,
  height = 2000,
  width = 4000,
  units = "px",
)


ggplot(
  data = filter(data_SAINT_merged_modified),
  aes(
    x = SP,
    y = FCA,
    label = PreyGene,
    fill = is_not_crapome,
    color = merged_activity
  )
) +
  facet_grid(~Bait) +
  geom_point(size = 2, shape = 21) +

  # Limits for high confidence interactions
  geom_hline(
    yintercept = 2,
    color = "darkcyan",
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = 0.77,
    y = 2.1,
    label = "FCA = 2",
    color = "darkcyan"
  ) +
  geom_vline(
    xintercept = 0.95,
    color = "darkcyan",
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = 0.945,
    y = 9,
    angle = 90,
    label = "SP = 0.95",
    color = "darkcyan"
  ) +
  scale_color_manual(
    values = c(
      "FALSE" = "white",
      "kinase" = "red",
      "kinase regulator activity" = "green"
    ),
    labels = c("FALSE" = " ")
  ) +
  scale_fill_manual(
    values = c("TRUE" = "black", "FALSE" = "gray"),
    labels = c("TRUE" = "Not crapome", "FALSE" = "Crapome")
  ) +
  coord_cartesian(xlim = c(0.75, 1), ylim = c(1, 10)) +
  scale_y_log10() +
  theme_bw() +
  labs(fill = NULL, color = NULL) + # Remove legend titles
  guides(
    fill = guide_legend(title = NULL), # Remove the fill legend title
    color = guide_legend(title = NULL)
  ) # Remove the color legend title

ggsave(
  fn_plot_sp_fca_crapome,
  height = 15,
  width = 40,
  units = "cm",
  dpi = 600
)

## Heatmap [NOT IN ARTICLE]----
# filter out prey genes that are significant according to parameters
data_significant_Prey <- data_SAINT_merged_modified %>%
  filter(FCA >= 2, SP > 0.95)
# single out names
data_significant_Prey <- unique(data_significant_Prey$Prey)

# Set up plot settings
pdat <- data_SAINT_merged_modified %>%
  filter(Prey %in% data_significant_Prey)

if (nrow(pdat) == 0) {
  plot_heatmaps <- ggplot() +
    theme_void() +
    ggtitle("No significant interactions found")
} else {
  plot_heatmap <- pdat %>%
    ggplot(
      aes(x = Bait, y = PreyGene)
    ) +
    # Gradient scale
    scale_fill_gradient(
      low = "blue",
      high = "red",
      guide = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        label.position = "bottom",
        ticks = TRUE
      ),
    ) +
    # Theme
    scale_x_discrete(expand = c(0, 0)) + # Remove space on the x-axis
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = -45, hjust = 0),
      # Tilt x-axis title
      legend.position = "top",
      legend.title = element_text(angle = 0, vjust = 1),
      legend.text = element_text(angle = 0),
      panel.background = element_rect(fill = "white"),
      # Gray background
      panel.grid.major = element_blank(),
      # No major grid lines
      panel.grid.minor = element_blank(),
      # No minor grid lines
      panel.border = element_rect(color = "black", fill = NA),
      # Box around the plot
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
    )

  # Generate different plots visualising different SAINT results
  plot_heatmap_FCA <- plot_heatmap +
    geom_tile(aes(fill = FCA))
  plot_heatmap_FCA

  plot_heatmap_SP <- plot_heatmap +
    geom_tile(aes(fill = SP))
  plot_heatmap_SP

  plot_heatmap_Abundance <- plot_heatmap +
    geom_tile(aes(fill = Abundance))
  plot_heatmap_Abundance


  # visualise all plots together
  plot_heatmaps <- grid.arrange(plot_heatmap_FCA, plot_heatmap_SP, plot_heatmap_Abundance, nrow = 1)
}
ggsave(
  fn_plot_heatmap,
  plot = plot_heatmaps,
  height = 2000, width = 3000, units = "px"
)
