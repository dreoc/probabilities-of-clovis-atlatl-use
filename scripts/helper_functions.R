# "Helper" scripts to convert saved analysis results into Table 1 and Figure 1 without refitting the models.
# This code needs code sourcing it to provide data_directory, HomeFolder, and frozen_dir.

library(ggplot2)
library(patchwork)

########################
# Read the chronology and morphology results including the bayes factor output, density output,
# and all pre-prepared input in the settings csv file. 
########################
bayes_factors <- read.csv(file.path(data_directory, "bayes_factors.csv"))
probabilities <- read.csv(file.path(data_directory, "probabilities.csv"))
density <- read.csv(file.path(data_directory, "density.csv"))
settings <- read.csv(file.path(HomeFolder, "data/settings.csv"))
s <- setNames(settings$value, settings$setting)

# get the absence:presence bayes factor value to be called by manuscript body in letter.rmd
bf <- list(absence_presence = bayes_factors$bayes_factor[bayes_factors$technology == "atlatl"])

########################
# For table 1. Arrange the point reference probabilities for comparison.
# the comparison is nested, first is the bow and atlatl, then thrusting, followed by throwing
########################
classLabels <- c(bow = "Bow", atlatl = "Atlatl", contact_thrusting = "Contact/thrusting",
  short_dual_use = "Short-range throwing/thrusting", hand_thrown = "Hand-thrown", Clovis = "Clovis")
classes <- names(classLabels)[names(classLabels) != "Clovis"]
stage_names <- c("stone", "add_thrusting", "add_nonlithic")
values <- matrix(NA_real_, 3, 5, dimnames = list(NULL, as.vector(classLabels[classes])))
for (i in seq_along(stage_names)) {
  z <- probabilities[probabilities$stage == stage_names[i], ]
  values[i, match(z$class, classes)] <- z$probability
}
table1_rawdata <- data.frame(
  Comparison = c("Stone tips", "Add thrusting", "Full comparison"),
  values, check.names = FALSE
)
table1_formatted <- table1_rawdata
for (j in 2:ncol(table1_formatted)) {
  table1_formatted[[j]] <- ifelse(is.na(table1_rawdata[[j]]), intToUtf8(0x2014),
                                sprintf("%.3f", table1_rawdata[[j]]))
}
table1 <- list(raw = table1_rawdata, display = table1_formatted)

########################
# Keep clovis and reference class colors consistent across figure 1's panels 
# Same for survival curves, using their own palette.
########################
technology_colors <- setNames(c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B"),
  as.vector(classLabels[c("Clovis", "bow", "atlatl", "hand_thrown", "short_dual_use", "contact_thrusting")]))
preservation_colors <- c("Wood survival" = "#1F77B4", "Site survival" = "#FF7F0E",
  "Combined survival" = "#111111")

########################
# Compute wood and site survival starting at the older Clovis boundary point from the settings values, 
# then plot it against (calibrated) calendar age.
# The combined curve multiplies wood and site survival under the model's assumptions.
########################
age <- sort(unique(c(seq(0, s[["clovis_old"]], by = 10),
                     s[["clovis_young"]], s[["clovis_old"]])))
wood <- s[["h"]] * exp(-s[["guyette_k"]] * age) +
  (1 - s[["h"]]) * exp(-s[["perez_k"]] * age)
site <- (1 + age / s[["site_offset"]])^(-s[["site_exponent"]])
survival <- rbind(
  data.frame(calendar_ka = (s[["clovis_old"]] - age) / 1000, value = wood,
             curve = "Wood survival"),
  data.frame(calendar_ka = (s[["clovis_old"]] - age) / 1000, value = site,
             curve = "Site survival"),
  data.frame(calendar_ka = (s[["clovis_old"]] - age) / 1000, value = wood * site,
             curve = "Combined survival"))
survival$curve <- factor(survival$curve, levels = names(preservation_colors))
clovis_old_ka <- s[["clovis_old"]] / 1000
clovis_young_ka <- s[["clovis_young"]] / 1000
survival_breaks <- seq(0, 1, by = 0.2)

########################
# The saved atlatl ages turned into SPD provides RCdate context in top panel of figure 1 (a).
# SPD is rescaled for the left axis. The right axis shows the density in original units.
########################
spd <- read.csv(file.path(frozen_dir, "figure_inputs/Figure_1_Atlatl_SPD_Eren_2026.csv"))
spd <- spd[spd$calendar_age_bp >= 0 & spd$calendar_age_bp <= s[["clovis_old"]], ]
spd$calendar_ka <- spd$calendar_age_bp / 1000
density_breaks <- pretty(c(0, max(spd$normalized_spd_per_ka)), n = 4)
density_breaks <- density_breaks[density_breaks >= 0]
density_scale <- max(density_breaks)
spd$plot_height <- spd$normalized_spd_per_ka / density_scale
line_colors <- c(preservation_colors, "Atlatl dates (right axis)" = technology_colors[["Atlatl"]])
line_types <- c("Wood survival" = "solid", "Site survival" = "22",
  "Combined survival" = "solid", "Atlatl dates (right axis)" = "solid")

preservation_plot <- ggplot() +
  annotate("rect", xmin = clovis_young_ka, xmax = clovis_old_ka,
           ymin = -Inf, ymax = Inf, fill = "gray70", alpha = 0.22) +
  geom_area(data = spd, aes(calendar_ka, plot_height),
            fill = technology_colors[["Atlatl"]], alpha = 0.18) +
  geom_line(data = spd, aes(calendar_ka, plot_height,
    color = "Atlatl dates (right axis)", linetype = "Atlatl dates (right axis)"),
    linewidth = 0.4) +
  geom_line(data = survival, aes(calendar_ka, value, color = curve,
    linetype = curve, linewidth = curve, group = curve)) +
  scale_linewidth_manual(values = setNames(c(0.85, 0.85, 1.05), names(preservation_colors)),
                         guide = "none") +
  annotate("text", x = clovis_old_ka, y = 1.065, label = "Clovis interval",
    hjust = 0, vjust = 0.5, family = "Times New Roman", size = 2.8) +
  scale_x_reverse(breaks = c(clovis_old_ka, seq(12, 0, by = -2)),
    labels = c(format(clovis_old_ka, trim = TRUE), as.character(seq(12, 0, by = -2))),
    limits = c(clovis_old_ka, 0), expand = c(0.01, 0)) +
  scale_y_continuous(name = "Survival probability since Clovis",
    breaks = survival_breaks, labels = sprintf("%.2f", survival_breaks),
    limits = c(0, 1.14), expand = c(0, 0),
    sec.axis = sec_axis(~ . * density_scale,
      name = "Atlatl calibrated age density\n(per ka)", breaks = density_breaks,
      labels = sprintf("%.2f", density_breaks))) +
  scale_color_manual(values = line_colors, breaks = names(line_colors), drop = FALSE,
    guide = guide_legend(nrow = 1, byrow = TRUE,
      override.aes = list(linetype = as.vector(line_types), linewidth = c(0.85, 0.85, 1.05, 0.4)))) +
  scale_linetype_manual(values = line_types, breaks = names(line_types),
                        drop = FALSE, guide = "none") +
  labs(title = "a  Preservation and dated atlatl record", x = "Calendar age (ka BP)",
       color = NULL, linetype = NULL) +
  theme_classic(base_family = "Times New Roman", base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0),
    axis.title = element_text(size = 8), axis.text = element_text(color = "black", size = 8),
    axis.line = element_line(color = "black", linewidth = 0.45),
    axis.ticks = element_line(color = "black", linewidth = 0.45),
    axis.title.y.right = element_text(color = technology_colors[["Atlatl"]], size = 8),
    axis.text.y.right = element_text(color = technology_colors[["Atlatl"]], size = 8),
    axis.line.y.right = element_line(color = technology_colors[["Atlatl"]], linewidth = 0.45),
    axis.ticks.y.right = element_line(color = technology_colors[["Atlatl"]], linewidth = 0.45),
    legend.position = "top", legend.title = element_blank(), legend.text = element_text(size = 7),
    legend.key.width = grid::unit(0.23, "in"), legend.background = element_blank(),
    legend.margin = margin(0, 0, 4, 0))

########################
# Use Table 1 class names and a fixed drawing 
# order for the predictive density curves.
########################
rc_density <- density
rc_density$class <- factor(as.vector(classLabels[rc_density$class]), levels = names(technology_colors))
line_types <- setNames(c("22", rep("solid", 5)), names(technology_colors))
panel_titles <- c("b  Stone tips", "c  Add thrusting", "d  Add nonlithic tips")
morphology_breaks <- seq(0, 2.5, by = 0.5)
panels <- list()

########################
# plot Clovis the three reference comparisons, using TCSA as log(mm^2) using ggplot 
# for panels b-d in figure 1.
# Axes are shared, crop as close as possible to tallest density (2.5).
########################
for (i in seq_along(stage_names)) {
  stage <- stage_names[i]
  panel_data <- rc_density[rc_density$stage == stage, ]
  panels[[stage]] <- ggplot(panel_data,
    aes(log_tcsa, density_log, color = class, linetype = class, group = class)) +
    geom_line(linewidth = 0.85) +
    scale_x_continuous(breaks = log(c(1, 10, 100, 1000)),
      labels = c("1", "10", "100", "1,000"), expand = c(0, 0)) +
    scale_y_continuous(breaks = morphology_breaks,
      labels = sprintf("%.2f", morphology_breaks), expand = c(0, 0)) +
    coord_cartesian(xlim = log(c(0.5, 2000)), ylim = c(0, 2.5), expand = FALSE) +
    scale_color_manual(values = technology_colors, breaks = names(technology_colors),
                       drop = FALSE, guide = "none") +
    scale_linetype_manual(values = line_types, breaks = names(technology_colors),
                         drop = FALSE, guide = "none") +
    labs(title = panel_titles[i], x = expression("TCSA (mm"^2*"; log scale)"),
      y = if (i == 1) "Density per unit log TCSA" else NULL, color = NULL, linetype = NULL) +
    theme_classic(base_family = "Times New Roman", base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0),
      legend.position = "none",
      axis.title = element_text(size = 8), axis.text = element_text(color = "black", size = 8),
      axis.line = element_line(color = "black", linewidth = 0.45),
      axis.ticks = element_line(color = "black", linewidth = 0.45))
}

########################
# Include all classes in one legend.
########################
legend_data <- expand.grid(x = c(0, 1), class = names(technology_colors), KEEP.OUT.ATTRS = FALSE)
legend_data$y <- 0
legend_data$class <- factor(legend_data$class, levels = names(technology_colors))
legend_plot <- ggplot(legend_data, aes(x, y, color = class, linetype = class, group = class)) +
  geom_line(linewidth = 0.85, alpha = 0) +
  scale_color_manual(values = technology_colors, breaks = names(technology_colors), drop = FALSE,
    guide = guide_legend(nrow = 1, byrow = TRUE,
                         override.aes = list(linetype = as.vector(line_types), alpha = 1))) +
  scale_linetype_manual(values = line_types, breaks = names(technology_colors),
                       drop = FALSE, guide = "none") +
  labs(color = NULL, linetype = NULL) +
  theme_void(base_family = "Times New Roman", base_size = 7) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal",
    legend.text = element_text(size = 7), legend.key.width = grid::unit(0.25, "in"),
    legend.spacing.x = grid::unit(0.08, "in"), legend.margin = margin(0), plot.margin = margin(0))

Figure_1 <- wrap_plots(preservation_plot, wrap_plots(panels, nrow = 1), legend_plot,
                       ncol = 1, heights = c(1.15, 1.00, 0.13))

########################
# Save figure after adjusting for figure size and resolution.
########################
figure_path <- file.path(HomeFolder, "manuscript/figures/Figure_1.png")
dir.create(dirname(figure_path), showWarnings = FALSE, recursive = TRUE)
ggsave(figure_path, Figure_1, width = 6.9, height = 5.7,
       units = "in", dpi = 600, bg = "white", device = ragg::agg_png)
