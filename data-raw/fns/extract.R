extract_guide_box_lgd <- function (plot_plt)
{
    tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plot_plt))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend_1L_chr <- tmp$grobs[[leg]]
    return(legend_1L_chr)
}
