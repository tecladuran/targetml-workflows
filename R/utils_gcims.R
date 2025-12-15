plot_sample <- function(sample_obj, colormap, minmax = NULL, minmax_for_legend = NULL, dt_min = NULL, dt_max = NULL, rt_min = NULL, rt_max = NULL){
  intmat <- intensity(sample_obj)
  
  if (is.null(dt_min)) {
    dt_min <- as.numeric(rownames(intmat)[1L])
  }
  if (is.null(dt_max)) {
    dt_max <- as.numeric(rownames(intmat)[nrow(intmat)])
  }
  if (is.null(rt_min)) {
    rt_min <- as.numeric(colnames(intmat)[1L])
  }
  if (is.null(rt_max)) {
    rt_max <- as.numeric(colnames(intmat)[ncol(intmat)])
  }
  
  trans <- cubic_root_trans()
  intmat_trans <- trans$transform(intmat)
  if (is.null(minmax)){
    minmax <- range(intmat_trans)
  }
  if (is.null(minmax_for_legend)){
    minmax_for_legend <- range(intmat)
  }
  nr <- build_nr(intmat_trans, minmax = c(minmax[1], minmax[2]), colormap)
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      xmin = dt_min, xmax = dt_min,
      ymin = rt_min, ymax = rt_min,
      ggplot2::aes(fill = .data$x),
      data = data.frame(
        x = NA_real_,
        dt_ms_min = dt_min, dt_ms_max = dt_max,
        rt_s_min = rt_min, rt_s_max = rt_max
      )
    ) +
    ggplot2::annotation_raster(
      nr,
      xmin = dt_min, xmax = dt_max,
      ymin = rt_min, ymax = rt_max
    ) +
    scale_fill_gradientn(colors = colormap,   # Custom Spectral colormap
                         limits = minmax_for_legend,              # Use your predefined limits
                         na.value = "#00000000",       # Transparent for NA values
                         trans = trans) +
    ggplot2::lims(
      x = c(dt_min, dt_max),
      y = c(rt_min, rt_max)
    ) +
    ggplot2::labs(
      x = "Drift time (ms)",
      y = "Retention time (s)",
      fill = "Intensity (a.u.)"
    ) +
    ggplot2::theme_minimal()
  print(p)
}

# Create colormap for plot_sample()
colormap <- c(
  colorRampPalette(c("black", "#26456E"))(70),          # black → dark blue
  colorRampPalette(c("#26456E", "#1C65A3"))(30),        # dark blue → blue
  colorRampPalette(c("#1C65A3", "#4993C0"))(20),        # blue → light blue
  colorRampPalette(c("#4993C0", "orange"))(20),         # light blue → orange
  colorRampPalette(c("orange", "red"))(60),             # orange → red
  colorRampPalette(c("red", "#CB1618"))(70)             # red → dark red
)


build_nr <- function(x, minmax = NULL, colormap){
  if (is.null(minmax)) {
    minmax <- range(x)
  }
  breaks <- seq(from = minmax[1], to = minmax[2], length.out = length(colormap))
  xdim <- dim(x)
  rev_cols <- seq.int(ncol(x), 1L, by = -1L)
  x <- x[, rev_cols]
  x <- findInterval(x, breaks, rightmost.closed = TRUE)
  x <- colormap[x]
  x <- matrix(x, nrow = xdim[1], ncol = xdim[2], byrow = FALSE)
  
  nr <- structure(
    x,
    dim = c(xdim[2], xdim[1]),
    class = "nativeRaster",
    channels = 4L
  )
  return(nr)
}
