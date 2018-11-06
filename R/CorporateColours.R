#' @importFrom grDevices colorRampPalette
# llsr corporate colors
llsr_colors <- c(
  `red`        = "#d11141",
  `green`      = "#00b159",
  `blue`       = "#00aedb",
  `orange`     = "#f37735",
  `yellow`     = "#ffc425",
  `light grey` = "#cccccc",
  `dark grey`  = "#8c8c8c",
  #
  `dark blue`        = "#396AB1",
  `darf orange`      = "#DA7C30",
  `dark green`       = "#3E9651",
  `gray`     = "#535154",
  `purple` = "#6B4C9A",
  `dark mustard`  = "#948B3D"
)
#
#' Function to extract llsr colors as hex codes
#'
#' @param ... Character names of llsr_colors 
#'
llsr_cols <- function(...) {
  cols <- c(...)
  
  if (is.null(cols))
    return (llsr_colors)
  
  llsr_colors[cols]
}
#
llsr_palettes <- list(
  `main`  = llsr_cols("blue", "green", "yellow"),
  
  `cool`  = llsr_cols("blue", "green"),
  
  `hot`   = llsr_cols("yellow", "orange", "red"),
  
  `mixed` = llsr_cols("blue", "green", "yellow", "orange", "red"),
  
  `grey`  = llsr_cols("light grey", "dark grey")
)
#
llsr_pal <- function(palette = "main", reverse = FALSE, ...) {
  pal <- llsr_palettes[[palette]]
  
  if (reverse) pal <- rev(pal)
  
  colorRampPalette(pal, ...)
}
#
#' Color scale constructor for llsr colors
#'
#' @param palette Character name of palette in llsr_palettes
#' @param discrete Boolean indicating whether color aesthetic is discrete or not
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_color_gradientn(), used respectively when discrete is TRUE or FALSE
#'
scale_color_llsr <-
  function(palette = "main",
           discrete = TRUE,
           reverse = FALSE,
           ...) {
    pal <- llsr_pal(palette = palette, reverse = reverse)
    
    if (discrete) {
      discrete_scale("colour", paste0("llsr_", palette), palette = pal, ...)
    } else {
      scale_color_gradientn(colours = pal(256), ...)
    }
  }

#' Fill scale constructor for llsr colors
#'
#' @param palette Character name of palette in llsr_palettes
#' @param discrete Boolean indicating whether color aesthetic is discrete or not
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_fill_gradientn(), used respectively when discrete is TRUE or FALSE
#'
scale_fill_llsr <-
  function(palette = "main",
           discrete = TRUE,
           reverse = FALSE,
           ...) {
    pal <- llsr_pal(palette = palette, reverse = reverse)
    
    if (discrete) {
      discrete_scale("fill", paste0("llsr_", palette), palette = pal, ...)
    } else {
      scale_fill_gradientn(colours = pal(256), ...)
    }
  }
