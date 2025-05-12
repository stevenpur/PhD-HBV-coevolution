library(tidyverse)
library(gridExtra)

file <- "~/hbv_covar3/analysis/sim_seq/treeso_roc_result.txt"
data <- read.table(file, header=TRUE, sep="\t")

data$sequence_length <- as.numeric(sub(".*l([0-9]+).*", "\\1", rownames(data)))
data$sample_size <- as.numeric(sub(".*n([0-9]+).*", "\\1", rownames(data)))
data$coev_factor <- as.numeric(sub(".*f([0-9.]+).*", "\\1", rownames(data)))
data$mutation_rate <- as.numeric(sub(".*u([0-9.]+).*", "\\1", rownames(data)))

# .x is the treeso result and .y is the dca result

metrics <- 'recall'
metric_vals <- data[[metrics]]
metric_vals2 <- metric_vals[-which(metric_vals == 1)]

metrics_range <- c(min(metric_vals2, na.rm = TRUE), max(metric_vals2, na.rm = TRUE))

# metrics_range <- c(min(data[[metrics]], na.rm = TRUE), max(data[[metrics]], na.rm = TRUE))

# set break points for the color scale
mybreaks <- seq(metrics_range[1], metrics_range[2] + metrics_range[2]/100, length.out = 9)

#Function to return the desired number of colors
mycolors <- function(x) {
   colors <- colorRampPalette(c("darkblue", "green", "yellow"))(8)
}

#Function to create the contour plot
create_contour_plot <- function(data, z, mybreaks, mycolors, title, position) {
   p <- ggplot(data, aes(x=log10(coev_factor), y=log10(mutation_rate), z=z)) +
      geom_contour_filled(breaks = mybreaks) +
      scale_fill_manual(name = title, values = mycolors(10), drop=FALSE) +
      scale_x_continuous(name="Log of Coevolution Factor (f)", 
                         breaks=log10(data$coev_factor), 
                         labels=data$coev_factor) +
      scale_y_continuous(name="Log of Mutation Rate (u)", 
                         breaks=log10(data$mutation_rate), 
                         labels=data$mutation_rate) +
      theme(legend.position = "none") +
      ggtitle(title)
   # if position is not left, then set the text to ""
   if (position != "left") {
      p <- p + theme(axis.title.y = element_text(color = "transparent"))
   }
   # if position is not middle, then hide the plot title (by setting it to empty string)
   if (position != "middle") {
      p <- p + theme(plot.title = element_text(color = "transparent"))
   }
   # if position is not right then hide the legend
   if (position == "legend") {
      p <- p + theme(legend.position = "right")
      p_legend <- ggplot_gtable(ggplot_build(p))
      legend <- p_legend$grobs[[which(p_legend$layout$name == "guide-box")]]
      return(legend)
   }
   return(p)
}

plots <- list()

unique_sample_sizes <- sort(unique(data$sample_size))
for (samp_size in unique_sample_sizes) {
  subset_data <- subset(data, sample_size == samp_size)
   # get the position of the plot
   position <- ""
   if (samp_size == unique_sample_sizes[1]) {
      position <- "left"
   } else if (samp_size == unique_sample_sizes[length(unique_sample_sizes)]) {
      position <- "right"
   } else if (samp_size == unique_sample_sizes[ceiling(length(unique_sample_sizes)/2)]) {
      position <- "middle"
   }
  # Creating contour plots
   p <- create_contour_plot(subset_data, subset_data[[metrics]], mybreaks, mycolors, metrics, position)
   grob <- ggplotGrob(p)
   plots[[length(plots) + 1]] <- grob
}
# adding legend
p_legend <- create_contour_plot(subset_data, subset_data[[metrics]], mybreaks, mycolors, metrics, "legend")
plots[[length(plots) + 1]] <- p_legend

options(repr.plot.width = 30, repr.plot.height = 2)
#grid.arrange(grobs = c(plots_tpr, plots_fpr, plots_indir), ncol = length(unique_sample_sizes) + 1)
g_result <- arrangeGrob(grobs = c(plots), nrow = 1)

ggsave(paste0("~/hbv_covar3/plots/contour_treeso_", metrics, "_plots.pdf"), g_result, width = 40, height = 10)


# Loop through each sample size and generate a contour plot for avg_tpr, avg_fpr, avg_indir
unique_sample_sizes <- sort(unique(data$sample_size))
for (samp_size in unique_sample_sizes) {
  subset_data <- subset(data, sample_size == samp_size)
   # get the position of the plot
   position <- ""
   if (samp_size == unique_sample_sizes[1]) {
      position <- "left"
   } else if (samp_size == unique_sample_sizes[length(unique_sample_sizes)]) {
      position <- "right"
   } else if (samp_size == unique_sample_sizes[ceiling(length(unique_sample_sizes)/2)]) {
      position <- "middle"
   }
  # Creating contour plots
   p_tpr <- create_contour_plot(subset_data, subset_data$avg_tpr, mybreaks_tpr, mycolors, "Average TPR", position)
   grob_tpr <- ggplotGrob(p_tpr)
   plots_tpr[[length(plots_tpr) + 1]] <- grob_tpr

   p_fpr <- create_contour_plot(subset_data, subset_data$avg_fpr, mybreaks_fpr, mycolors, "Average FPR", position)
   grob_fpr <- ggplotGrob(p_fpr)
   plots_fpr[[length(plots_fpr) + 1]] <- grob_fpr

   p_indir <- create_contour_plot(subset_data, subset_data$avg_indir, mybreaks_indir, mycolors, "Average Indirect Coevolution", position)
   grob_indir <- ggplotGrob(p_indir)
   plots_indir[[length(plots_indir) + 1]] <- grob_indir
}

# Adding legend to each row
p_legend_tpr <- create_contour_plot(subset_data, subset_data$avg_tpr, mybreaks_tpr, mycolors, "Average TPR", "legend")
p_legend_fpr <- create_contour_plot(subset_data, subset_data$avg_fpr, mybreaks_fpr, mycolors, "Average FPR", "legend")
p_legend_indir <- create_contour_plot(subset_data, subset_data$avg_indir, mybreaks_indir, mycolors, "Average Indirect Coevolution", "legend")
plots_tpr[[length(plots_tpr) + 1]] <- p_legend_tpr
plots_fpr[[length(plots_fpr) + 1]] <- p_legend_fpr
plots_indir[[length(plots_indir) + 1]] <- p_legend_indir

# plot results
options(repr.plot.width = 12, repr.plot.height = 10)

#grid.arrange(grobs = c(plots_tpr, plots_fpr, plots_indir), ncol = length(unique_sample_sizes) + 1)
g_result <- arrangeGrob(grobs = c(plots_tpr, plots_fpr, plots_indir), ncol = length(unique_sample_sizes) + 1)

ggsave("~/test.pdf", g_result)

