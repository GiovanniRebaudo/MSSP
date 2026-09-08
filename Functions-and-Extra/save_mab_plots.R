# Export the three existing MAB comparison panels, without running any sampler.

mab_result_plots <- function(experiment = c("simul", "trees"), envir = parent.frame()){
  experiment = match.arg(experiment)
  objects = c(Uniform = "results_random", `Ind DP` = "results_indepDP",
              `Ind PY` = "results_indepPY", `+DP` = "results_plusDP",
              `+PY` = "results_plusPY", HDP = "results_HDP", HPY = "results_HPY")
  if(experiment == "simul") objects = c(objects, Oracle = "results_oracle")
  if(experiment == "trees") objects = stats::setNames(paste0(objects, "_real"),
                                                       names(objects))

  missing = objects[!vapply(objects, exists, logical(1), envir = envir, inherits = FALSE)]
  if(length(missing)) stop("Missing results: ", paste(missing, collapse = ", "),
                          ". Wait for the run to finish or load its saved workspace.")
  results = lapply(objects, get, envir = envir, inherits = FALSE)
  valid = vapply(results, function(x)
    is.matrix(x) && is.numeric(x) && nrow(x) > 0L && ncol(x) > 0L &&
      all(is.finite(x)) && identical(dim(x), dim(results[[1]])), logical(1))
  if(!all(valid)) stop("Results must be complete, finite matrices of matching dimensions.")
  curves = lapply(results, rowMeans, na.rm = TRUE)
  n_steps = nrow(results[[1]])

  panels = list(independent = c("Uniform", "Ind DP", "Ind PY"),
                additive = c("Uniform", "+DP", "+PY"),
                hierarchical = c("Uniform", "HDP", "HPY"))
  if(experiment == "simul") panels = lapply(panels, c, "Oracle")
  titles = c(independent = "Independent Processes", additive = "Additive processes",
             hierarchical = "Hierarchical Processes")
  plots = lapply(names(panels), function(panel){
    models = panels[[panel]]
    data_plot = data.frame(time = rep(seq_len(n_steps), length(models)),
                           model = rep(models, each = n_steps),
                           value = unlist(curves[models], use.names = FALSE))
    plot = ggplot2::ggplot(data_plot, ggplot2::aes(x = time, y = value,
                                                 color = as.factor(model))) +
      ggplot2::geom_line(ggplot2::aes(linetype = as.factor(model)), linewidth = 1.2) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Additional Samples", y = "Discoveries") +
      ggplot2::scale_color_brewer(palette = "Dark2") +
      ggplot2::theme(text = ggplot2::element_text(size = 20),
                     legend.position = "right", legend.title = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::ggtitle(titles[[panel]])
    if(experiment == "trees") plot = plot + ggplot2::scale_y_continuous(limits = c(0, 100))
    plot
  })
  names(plots) = names(panels)
  plots
}

save_mab_plots <- function(experiment = c("simul", "trees"),
                            output_dir = file.path("Data-and-Results", "MAB-plots"),
                            width = 8, height = 6, overwrite = FALSE,
                            envir = parent.frame()){
  experiment = match.arg(experiment)
  plots = mab_result_plots(experiment, envir = envir)
  stopifnot(length(width) == 1L, is.finite(width), width > 0,
            length(height) == 1L, is.finite(height), height > 0,
            is.logical(overwrite), length(overwrite) == 1L, !is.na(overwrite))
  paths = file.path(output_dir, paste0("MAB_", experiment, "_", names(plots), ".pdf"))
  if(!overwrite && any(file.exists(paths)))
    stop("Some PDFs already exist. Use overwrite = TRUE to replace them, or choose another output_dir.")
  if(!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE))
    stop("Cannot create output directory: ", output_dir)
  for(i in seq_along(plots)){
    ggplot2::ggsave(filename = paths[i], plot = plots[[i]], device = "pdf",
                    width = width, height = height, units = "in",
                    bg = "white", useDingbats = FALSE)
  }
  paths = stats::setNames(normalizePath(paths, mustWork = TRUE), names(plots))
  message("Saved PDFs:\n", paste(paths, collapse = "\n"))
  invisible(paths)
}
