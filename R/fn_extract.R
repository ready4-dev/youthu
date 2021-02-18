#' Extract guide box lgd
#' @description extract_guide_box_lgd() is an Extract function that extracts data from an object. Specifically, this function implements an algorithm to extract guide box lgd. The function is called for its side effects and does not return a value.
#' @param a.gplot PARAM_DESCRIPTION
#' @return NA ()
#' @rdname extract_guide_box_lgd
#' @export 
#' @importFrom ggplot2 ggplot_gtable ggplot_build
#' @keywords internal
extract_guide_box_lgd <- function (a.gplot) 
{
    tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}
