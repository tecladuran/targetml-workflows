#' Identify target VOCs from cluster stats with reliability warning
#'
#' @param cluster_stats Dataframe with detected clusters. Must contain columns: cluster, dt_apex_ms, rt_apex_s.
#' @param target_vocs Dataframe with expected VOCs. Must contain columns: compound, dt, rt.
#' @param max_distance Threshold to consider a match as reliable. Matches with a distance below this threshold are considered reliable.
#' @param ratio A value to adjust the importance of drift time (dt) relative to retention time (rt). The higher the ratio, the more important the drift time becomes in the distance calculation.
#'
#' @return A dataframe with target VOCs matched to clusters, including reliability flag.
#'
#' @export
identify_target_vocs <- function(cluster_stats, target_vocs, max_distance = 2, ratio = 1) {
  
  # Match all VOCs to nearest cluster
  results <- purrr::map_dfr(seq_len(nrow(target_vocs)), function(i) {
    voc <- target_vocs[i, ]
    
    # Calculate distance based on ratio
    cluster_stats <- cluster_stats %>%
      dplyr::mutate(
        # Calculate distance, scaling RT difference by the ratio
        distance = sqrt((dt_apex_ms - voc$dt)^2 + (ratio*(rt_apex_s - voc$rt))^2)
      )
    
    closest <- cluster_stats %>%
      dplyr::arrange(distance) %>%
      dplyr::slice(1)
    
    match <- closest %>%
      dplyr::select(cluster, dt_apex_ms, rt_apex_s, distance) %>%
      dplyr::mutate(reliable_match = distance <= max_distance)
    
    dplyr::bind_cols(voc, match)
  })
  
  # Warn if some matches are not reliable
  unreliable <- results %>% dplyr::filter(!reliable_match)
  if (nrow(unreliable) > 0) {
    warning(
      "Unreliable matches found for the following VOCs:\n- ",
      paste(unreliable$compound, collapse = "\n- "),
      "\n\nThese VOCs were matched to the closest cluster despite exceeding the `max_distance` threshold.\n",
      "You may want to:\n",
      "- Use `modify_voc_assignment(results, voc_name, new_cluster, cluster_stats)` to manually correct a match.\n",
      "- Use `remove_voc(results, voc_name)` to remove a VOC from the results.\n"
    )
  }
  
  return(results)
}

modify_voc_assignment <- function(results, voc_name, new_cluster, cluster_stats, ratio = 1) {
  row_idx <- which(results$compound == voc_name)
  cluster_row <- cluster_stats %>% dplyr::filter(cluster == new_cluster)
  
  if (length(row_idx) == 1 && nrow(cluster_row) == 1) {
    results$cluster[row_idx] <- new_cluster
    results$dt_apex_ms[row_idx] <- cluster_row$dt_apex_ms
    results$rt_apex_s[row_idx] <- cluster_row$rt_apex_s
    results$distance[row_idx] <- sqrt(
      (results$dt[row_idx] - cluster_row$dt_apex_ms)^2 +
        (ratio*(results$rt[row_idx] - cluster_row$rt_apex_s))^2
    )
    results$reliable_match[row_idx] <- TRUE
  } else {
    message("VOC name or cluster not found. No modifications made.")
  }
  
  return(results)
}



remove_voc <- function(results, voc_name) {
  results <- results %>% dplyr::filter(compound != voc_name)
  return(results)
}