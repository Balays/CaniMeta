suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

n50 <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & x > 0]
  if (!length(x)) return(NA_real_)
  x <- sort(x, decreasing = TRUE)
  x[which(cumsum(x) >= sum(x)/2)[1]]
}

theme_pub <- function(base_size = 10) {
  theme_ipsum(base_size = base_size) +
    theme(plot.title = element_text(face = "bold"))
}

save_plot <- function(p, out_dir, name, w = 170, h = 120, dpi = 300,  ...) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", limitsize = FALSE, device = cairo_pdf)
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, units = "mm", limitsize = FALSE, dpi = dpi, ...)
}

summarize_ont_run <- function(dt,
                              out_dir = "ONT_QC",
                              sample_name = "Dog_M0 MinION WGS (LSK114)",
                              length_cap = 200000,
                              q_bins = 80,
                              len_bins = 100) {

  if (!is.data.table(dt)) dt <- as.data.table(dt)

  # required columns
  req <- c("sequence_length_template", "mean_qscore_template", "passes_filtering")
  stopifnot(all(req %in% names(dt)))

  # coerce types
  dt[, passes_filtering := as.logical(passes_filtering)]
  dt[, sequence_length_template := as.numeric(sequence_length_template)]
  dt[, mean_qscore_template := as.numeric(mean_qscore_template)]

  # optional columns
  if (!("alias" %in% names(dt))) dt[, alias := "unknown"]
  if (!("start_time" %in% names(dt))) dt[, start_time := NA_real_]
  dt[, alias := fifelse(is.na(alias) | alias == "", "unknown", as.character(alias))]
  dt[, start_time := as.numeric(start_time)]

  # derived
  dt[, pass_fail := fifelse(passes_filtering, "pass", "fail")]

  # ---- overall
  overall <- dt[, .(
    n_reads = .N,
    n_pass = sum(passes_filtering, na.rm = TRUE),
    n_fail = sum(!passes_filtering, na.rm = TRUE),
    pass_rate = sum(passes_filtering, na.rm = TRUE) / .N,
    total_bases = sum(sequence_length_template, na.rm = TRUE),
    mean_len = mean(sequence_length_template, na.rm = TRUE),
    median_len = median(sequence_length_template, na.rm = TRUE),
    n50_len = n50(sequence_length_template),
    p90_len = as.numeric(quantile(sequence_length_template, 0.90, na.rm = TRUE)),
    mean_q = mean(mean_qscore_template, na.rm = TRUE),
    median_q = median(mean_qscore_template, na.rm = TRUE)
  ), by = .(sample)]
  #overall[, sample := sample_name]
  setcolorder(overall, c("sample", setdiff(names(overall), "sample")))

  # ---- by pass/fail
  by_pf <- dt[, .(
    n_reads = .N,
    total_bases = sum(sequence_length_template, na.rm = TRUE),
    mean_len = mean(sequence_length_template, na.rm = TRUE),
    median_len = median(sequence_length_template, na.rm = TRUE),
    n50_len = n50(sequence_length_template),
    mean_q = mean(mean_qscore_template, na.rm = TRUE),
    median_q = median(mean_qscore_template, na.rm = TRUE)
  ), by = .(sample, pass_fail)][order(pass_fail, sample)]
  #by_pf[, sample := sample_name]
  setcolorder(by_pf, c("sample","pass_fail", setdiff(names(by_pf), c("sample","pass_fail"))))

  # ---- by alias (barcode)
  by_alias <- dt[, .(
    n_reads = .N,
    pass_rate = mean(passes_filtering, na.rm = TRUE),
    total_bases = sum(sequence_length_template, na.rm = TRUE),
    mean_len = mean(sequence_length_template, na.rm = TRUE),
    median_len = median(sequence_length_template, na.rm = TRUE),
    n50_len = n50(sequence_length_template),
    mean_q = mean(mean_qscore_template, na.rm = TRUE)
  ), by = .(sample, alias)][order(-n_reads)]
  #by_alias[, sample := sample_name]
  setcolorder(by_alias, c("sample","alias", setdiff(names(by_alias), c("sample","alias"))))

  # write tables
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(overall,  paste0(out_dir, '/', sample_name, "_overall_summary.tsv"),     sep = '\t')
  fwrite(by_pf,    paste0(out_dir, '/', sample_name, "_summary_by_passfail.tsv"), sep = '\t')
  fwrite(by_alias, paste0(out_dir, '/', sample_name, "_summary_by_alias.tsv"),    sep = '\t')

  # ---- plots
  dt_plot <- copy(dt)
  dt_plot[, len_plot := pmin(sequence_length_template, length_cap)]

  pf_counts <- dt_plot[, .N, by = .(sample, pass_fail)]
  p_pass <- ggplot(pf_counts, aes(x = pass_fail, y = N)) +
    geom_col(width = 0.7) +
    scale_y_continuous(labels = label_comma()) +
    labs(title = paste0(sample_name, ": pass/fail reads"), x = NULL, y = "Reads") +
    theme_pub() +
    facet_nested_wrap(~sample, nrow = 1)

  save_plot(p_pass, out_dir, paste0(sample_name, "_pass_fail_counts"), w = 140, h = 110)

  p_len <- ggplot(dt_plot[is.finite(len_plot) & len_plot > 0], aes(x = len_plot, fill = pass_fail)) +
    geom_histogram(bins = len_bins, alpha = 0.6, position = "identity") +
    scale_x_continuous(labels = label_comma()) +
    scale_y_continuous(labels = label_comma()) +
    labs(title = paste0(sample_name, ": read length distribution (capped)"),
         x = "Read length (bp)", y = "Reads", fill = NULL) +
    theme_pub()+
    facet_nested_wrap(~sample, nrow = 1)

  save_plot(p_len, out_dir, paste0(sample_name, "_read_length_hist"), w = 170, h = 120)

  p_q <- ggplot(dt_plot[is.finite(mean_qscore_template)], aes(x = mean_qscore_template, fill = pass_fail)) +
    geom_histogram(bins = q_bins, alpha = 0.6, position = "identity") +
    scale_y_continuous(labels = label_comma()) +
    labs(title = paste0(sample_name, ": mean Q-score distribution"),
         x = "Mean Q-score", y = "Reads", fill = NULL) +
    theme_pub()+
    facet_nested_wrap(~sample, nrow = 1)

  save_plot(p_q, out_dir, paste0(sample_name, "_qscore_hist"), w = 170, h = 120)

  # yield over time (if start_time mostly present)
  if (sum(is.finite(dt_plot$start_time)) > 0.9 * nrow(dt_plot)) {
    t0 <- min(dt_plot$start_time, na.rm = TRUE)
    dt_plot[, minutes_from_start := (start_time - t0) / 60]
    dt_plot[, time_bin := floor(minutes_from_start)]
    yield <- dt_plot[, .(bases = sum(sequence_length_template, na.rm = TRUE)), by = .(sample, time_bin)][order(time_bin)]

    p_yield <- ggplot(yield, aes(x = time_bin, y = bases)) +
      geom_line(linewidth = 0.6) +
      scale_y_continuous(labels = label_number(scale_cut = cut_si(""), accuracy = 1)) +
      labs(title = paste0(sample_name, ": yield over time"),
           x = "Minutes from run start", y = "Bases per minute") +
      theme_pub()+
      facet_nested_wrap(~sample, nrow = 1)

    save_plot(p_yield, out_dir, paste0(sample_name, "_yield_over_time"), w = 180, h = 120)

    fwrite(yield, paste0(out_dir, '/', sample_name, "_yield_over_time_by_minute.tsv"), sep = '\t')
  }

  invisible(list(overall = overall, by_passfail = by_pf, by_alias = by_alias, dt_plot = dt_plot))
}
