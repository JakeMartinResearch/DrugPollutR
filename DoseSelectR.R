select_doses <- function(data,
                         variable,
                         central_tendency = "median",
                         min_spacing = 3.2,
                         max_doses = 5,
                         LOQ = 0.1,
                         limit_to_observed = FALSE) {
  
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  
  # Extract values 
  values <- data[[variable]]
  values <- values[!is.na(values) & values > 0]
  if (length(values) < 3) stop("Not enough data points.")
  
  # Central tendency
  center <- switch(tolower(central_tendency),
                   "mean" = mean(values),
                   "median" = median(values),
                   "mode" = {
                     d <- density(values)
                     d$x[which.max(d$y)]
                   },
                   stop("Invalid central_tendency"))
  
  upper <- quantile(values, 0.95)
  lower <- quantile(values, 0.05)
  max_val <- max(values)
  min_val <- min(values)
  
  # Generate doses above and below
  above <- c(); below <- c()
  
  val <- center
  while(length(above) < max_doses){
    val <- val * min_spacing
    if(limit_to_observed && val > max_val) break
    above <- c(above, val)
  }
  
  val <- center
  while(length(below) < max_doses){
    val <- val / min_spacing
    if(val <= LOQ) break
    if(limit_to_observed && val < min_val) break
    below <- c(below, val)
  }
  
  dose_list <- c(center)
  i <- 1
  while(length(dose_list) < max_doses){
    if(i <= length(above)) dose_list <- c(dose_list, above[i])
    if(length(dose_list) == max_doses) break
    if(i <= length(below)) dose_list <- c(dose_list, below[i])
    if(length(dose_list) == max_doses) break
    i <- i + 1
    if(i > 20) break
  }
  
  if(length(dose_list) < max_doses){
    warning("Unable to select enough doses within constraints.")
    return(NA)
  }
  
  dose_list <- sort(dose_list)
  dose_df <- data.frame(dose = paste0("D", seq_along(dose_list)),
                        value = dose_list)
  
  summary_tbl <- data.frame(
    central = center,
    upper_95 = upper,
    lower_05 = lower,
    max = max_val,
    min = min_val,
    LOQ = LOQ,
    row.names = NULL
  )
  summary_tbl$distribution_covered <- mean(values >= min(dose_list) & values <= max(dose_list))
  
  # ---------------- Density plots ----------------
  # Raw scale
  dens <- density(values, from = min(values), to = max(values), n = 512)
  dens_df <- data.frame(x = dens$x, y = dens$y) %>%
    mutate(in_cri = x >= lower & x <= upper)
  
  base_plot <- ggplot(dens_df, aes(x = x, y = y)) +
    geom_line(size = 1.2) +
    geom_ribbon(aes(ymin = 0, ymax = ifelse(in_cri, y, 0)), fill = "grey80", alpha = 0.5) +
    geom_vline(xintercept = center, color = "blue", linetype = "dashed", size = 1) +
    geom_vline(xintercept = LOQ, color = "red", linetype = "dotdash") +
    geom_point(data = dose_df, aes(x = value, y = 0), size = 3, color = "black") +
    geom_text(data = dose_df, aes(x = value, y = max(dens_df$y)*0.02, label = dose), vjust = -0.5) +
    ggtitle("Dose Selection (Raw Scale)") +
    theme_minimal()
  
  # Log scale
  log_values <- log10(values)
  dens_log <- density(log_values, n = 512)
  dens_df_log <- data.frame(x = 10^dens_log$x, y = dens_log$y)
  
  log_plot <- ggplot(dens_df_log, aes(x = x, y = y)) +
    geom_line(size = 1.2) +
    annotate("rect", xmin = lower, xmax = upper, ymin = 0, ymax = max(dens_df_log$y),
             fill = "grey80", alpha = 0.4) +
    geom_vline(xintercept = center, color = "blue", linetype = "dashed", size = 1) +
    geom_vline(xintercept = LOQ, color = "red", linetype = "dotdash") +
    geom_point(data = dose_df, aes(x = value, y = 0), size = 3, color = "black") +
    geom_text(data = dose_df, aes(x = value, y = max(dens_df_log$y)*0.02, label = dose), vjust = -0.5) +
    scale_x_log10(
      breaks = scales::log_breaks(n = 6),
      labels = scales::label_number(accuracy = 0.01)
    ) +
    ggtitle("Dose Selection (Log Scale)") +
    theme_minimal()
  
  # Patchwork combine
  combined_plot <- base_plot + log_plot
  
  return(list(doses = dose_df,
              summary = summary_tbl,
              plot = combined_plot))
}