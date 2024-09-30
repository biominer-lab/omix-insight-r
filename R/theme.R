#' Custom Theme
#'
#' A custom theme for ggplot2
#' @return custom theme
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_line
#' @export
custom_theme <- function() {
  custom <- theme_minimal() + theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    axis.title.x = element_text(color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1, size = 16),
    axis.text.x = element_text(color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 45, size = 14),
    axis.title.y = element_text(color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 450, size = 16),
    axis.text.y = element_text(color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 450, size = 14),
    legend.title = element_text(colour = "black", angle = 1, face = "bold", size = 16),
    legend.text = element_text(colour = "black", angle = 1, face = "bold", size = 14),
    legend.position = "right", legend.justification = c(1, 1),
    strip.text.x = element_text(face = "bold", size = 16),
    strip.text.y = element_text(face = "bold", size = 16),
    panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1),
    panel.grid.minor = element_blank(), axis.ticks = element_line(color = "black", size = 0.5)
  ) + theme(
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    # Set the whole plot area to white
    plot.background = element_rect(fill = "white", color = NA),
  )

  return(custom)
}
