#' Extract guide box legend
#' @description extract_guide_box_lgd() is an Extract function that extracts data from an object. Specifically, this function implements an algorithm to extract guide box legend. The function returns Legend (a character vector of length one).
#' @param plot_plt Plot (a plot)
#' @return Legend (a character vector of length one)
#' @rdname extract_guide_box_lgd
#' @export 
#' @importFrom ggplot2 ggplot_gtable ggplot_build
#' @keywords internal
extract_guide_box_lgd <- function (plot_plt) 
{
    tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plot_plt))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend_1L_chr <- tmp$grobs[[leg]]
    return(legend_1L_chr)
}
