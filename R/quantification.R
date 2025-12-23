# quantification.R — GC-IMS Quantification and Calibration Utilities
#
# This script implements utilities for the quantification of GC-IMS signals based on
# external calibration data. The workflow is designed to estimate analyte concentrations
# from measured peak intensities using calibrated reference curves and optional scaling
# and offset correction.
#
# The quantification strategy includes:
# - Interpolation of concentrations from calibration curves within the valid intensity range.
# - Estimation of scale and offset parameters to align target datasets to a reference.
# - Polynomial calibration modeling to map concentration–intensity relationships.
# - Robust handling of missing analytes and flexible analyte-to-peak mapping.
#
# The main entry point is `quantify_gcims()`, which applies the full quantification pipeline
# starting from a reference calibration dataset and a GC-IMS peak table.
#
# Additional utilities are provided to support error analysis and variance estimation.
#
# Typical usage:
#   quantified_table <- quantify_gcims(ref, analytes, scales, peak_table)



interpolate_concentration <- function(intensity, calibration_intensity, calibration_concentration) {
  # extreure els límits del calibratge
  min_int <- min(calibration_intensity, na.rm = TRUE)
  max_int <- max(calibration_intensity, na.rm = TRUE)
  
  # aproximem només dins del rang (rule = 1 → NA fora de [min, max])
  est <- approx(
    x    = calibration_intensity,
    y    = calibration_concentration,
    xout = intensity,
    rule = 1
  )$y
  
  # per sota del mínim: posem 0
  est[intensity < min_int] <- 0
  
  return(est)
}


# Finding Scales

find_offset_scale <- function(ref, target, analytes, fix_offset = FALSE) {
  results <- data.frame()
  
  for (analyte in analytes) {
    c_col <- paste0("c_", analyte)
    
    # Comprovació que les columnes existeixen
    if (!(analyte %in% names(ref)) || !(c_col %in% names(ref)) || !(analyte %in% names(target))) {
      stop(paste("Missing columns for analyte:", analyte))
    }
    
    # Crear dataframe temporal amb la concentració
    ref_tmp <- ref %>%
      dplyr::select(!!c_col, !!analyte) %>%
      dplyr::rename(concentration = !!rlang::sym(c_col))
    
    # Optimització per aquest analit
    res <- optimize_offset_scale_2(
      df_fixed = ref_tmp,
      df_to_align = target,
      analyte = analyte,
      align_to_origin = fix_offset
    )
    
    results <- rbind(results, data.frame(
      analyte = analyte,
      scale = round(res$scale, 3),
      offset = round(res$offset, 3)
    ))
  }
  
  return(results)
}

# Full Quantification

quantify_gcims <- function(
    ref,
    analytes,
    scales = rep(1, length(analytes)),
    peak_table,
    mapping = NULL
) {

  
  # Safety checks
  if (!is.data.frame(ref) || !is.data.frame(peak_table))
    stop("Both 'ref' and 'peak_table' must be data frames.")
  
  if (length(scales) != length(analytes))
    stop("Length of 'scales' must match length of 'analytes'.")
  
  # Normalize analyte names
  analytes_norm <- normalize_name(analytes)
  ref_names_norm <- normalize_name(names(ref))
  peak_names_norm <- normalize_name(names(peak_table))
  
  # Identify analytes that exist in both ref and peak_table
  valid_analytes <- c()
  missing_analytes <- c()
  
  for (i in seq_along(analytes_norm)) {
    analyte_norm <- analytes_norm[i]
    analyte_orig <- analytes[i]
    
    c_col_norm <- normalize_name(paste0("c_", analyte_orig))
    has_ref <- analyte_norm %in% ref_names_norm
    has_conc <- c_col_norm %in% ref_names_norm
    
    has_peak <- analyte_norm %in% peak_names_norm || 
      (!is.null(mapping) && analyte_orig %in% names(mapping))
    
    if (has_ref && has_conc && has_peak) {
      valid_analytes <- c(valid_analytes, analyte_orig)
    } else {
      missing_analytes <- c(missing_analytes, analyte_orig)
    }
  }
  
  # Report what will be quantified
  if (length(valid_analytes) > 0) {
    message("Quantifying the following analytes: ", paste(valid_analytes, collapse = ", "))
  } else {
    stop("No valid analytes found in both reference and peak_table.")
  }
  if (length(missing_analytes) > 0) {
    message("Skipping missing analytes: ", paste(missing_analytes, collapse = ", "))
  }
  
  # Prepare output container
  # Identify sample ID column (common cases)
  id_col <- names(peak_table)[tolower(names(peak_table)) %in% c("sampleid", "sample_id", "id")]
  if (length(id_col) == 0) id_col <- NULL
  
  # Initialize results with sample IDs if available
  if (!is.null(id_col)) {
    quantified <- data.frame(SampleID = peak_table[[id_col[1]]])
  } else {
    quantified <- data.frame(SampleID = seq_len(nrow(peak_table)))
  }
  
  
  for (i in seq_along(valid_analytes)) {
    analyte <- valid_analytes[i]
    scale_factor <- scales[i]
    
    analyte_norm <- normalize_name(analyte)
    c_col <- paste0("c_", analyte)
    c_col_norm <- normalize_name(c_col)
    
    # Match actual column names by normalized versions
    ref_col <- names(ref)[which(ref_names_norm == analyte_norm)]
    conc_col <- names(ref)[which(ref_names_norm == c_col_norm)]
    intensity_col <- NULL
    
    # Match in peak_table or use mapping
    if (!is.null(mapping) && analyte %in% names(mapping)) {
      intensity_col <- mapping[[analyte]]
    } else {
      match_idx <- which(peak_names_norm == analyte_norm)
      if (length(match_idx) > 0) intensity_col <- names(peak_table)[match_idx[1]]
    }
    
    if (is.null(intensity_col))
      next  # skip silently (shouldn't happen since pre-filtered)
    
    # Prepare reference calibration data
    df_ref <- ref %>%
      select(!!sym(conc_col), !!sym(ref_col)) %>%
      rename(concentration = !!sym(conc_col), intensity = !!sym(ref_col))
    
    # Fit polynomial model
    fit_poly <- lm(intensity ~ poly(concentration, 3, raw = TRUE) - 1, data = df_ref)
    
    # Generate scaled calibration curve
    conc_seq <- seq(0, max(df_ref$concentration), length.out = 300)
    predicted_intensity <- predict(fit_poly, newdata = data.frame(concentration = conc_seq))
    
    fit_df <- data.frame(
      intensity = predicted_intensity,
      conc_scaled = conc_seq / scale_factor
    )
    
    # Interpolate concentrations for measured intensities
    measured_intensities <- peak_table[[intensity_col]]
    
    conc_est <- interpolate_concentration(
      measured_intensities,
      fit_df$intensity,
      fit_df$conc_scaled
    )
    
    # Append results
    quantified[[paste0("conc_", analyte)]] <- conc_est
  }
  
  return(quantified)
}




# FOR ERRORS

weighted_within_var <- function(data, group_var, value_var) {
  data %>%
    group_by({{group_var}}) %>%
    summarise(
      n           = n(),
      var_within  = var({{value_var}}, na.rm = TRUE),
      .groups     = "drop"
    ) %>%
    summarise(weighted_var = sum(n / sum(n) * var_within)) %>%
    pull(weighted_var)
}