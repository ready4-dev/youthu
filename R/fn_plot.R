#' Plot auto lm
#' @description plot_auto_lm() is a Plot function that plots data Specifically, this function implements an algorithm to plot auto lm. The function is called for its side effects and does not return a value.
#' @param ... Additional arguments
#' @param which_dbl Which (a double vector), Default: 1:6
#' @param ncol_1L_int Ncol (an integer vector of length one), Default: 3
#' @param label_size_1L_int Label size (an integer vector of length one), Default: 3
#' @return NULL
#' @rdname plot_auto_lm
#' @export 
#' @importFrom pacman p_load
#' @importFrom ggplot2 autoplot
plot_auto_lm <- function (mdl, which_dbl = 1:6, ncol_1L_int = 3L, label_size_1L_int = 3) 
{
    pacman::p_load(char = "ggfortify")
    plt <- ggplot2::autoplot(mdl, which = which_dbl, ncol = ncol_1L_int, 
        label.size = label_size_1L_int)
    if (6 %in% which_dbl) 
        plt[which(which_dbl == 6)] <- plt[which(which_dbl == 
            6)] + ggtitle("Cook's vs Leverage")
    plt
}
#' Plot lnr cmprsn
#' @description plot_lnr_cmprsn() is a Plot function that plots data Specifically, this function implements an algorithm to plot lnr cmprsn. The function is called for its side effects and does not return a value.
#' @param data_tb Data (a tibble)
#' @param predn_ds_tb Predn dataset (a tibble)
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'aqol6d_total_w'
#' @param predr_var_nm_1L_chr Predr var name (a character vector of length one)
#' @param dep_var_desc_1L_chr Dep var description (a character vector of length one), Default: 'AQoL-6D utility score'
#' @param predr_var_desc_1L_chr Predr var description (a character vector of length one)
#' @return NULL
#' @rdname plot_lnr_cmprsn
#' @export 
#' @importFrom dplyr filter
#' @importFrom rlang sym
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth geom_line theme_bw labs
plot_lnr_cmprsn <- function (data_tb, predn_ds_tb, dep_var_nm_1L_chr = "aqol6d_total_w", 
    predr_var_nm_1L_chr, dep_var_desc_1L_chr = "AQoL-6D utility score", 
    predr_var_desc_1L_chr) 
{
    data_tb <- data_tb %>% dplyr::filter(!is.na(!!rlang::sym(predr_var_nm_1L_chr)))
    ggplot2::ggplot(data_tb, ggplot2::aes(x = !!rlang::sym(predr_var_nm_1L_chr), 
        y = !!rlang::sym(dep_var_nm_1L_chr))) + ggplot2::geom_point() + 
        ggplot2::geom_smooth(method = "loess", size = 1.5) + 
        ggplot2::geom_line(data = predn_ds_tb, ggplot2::aes(x = !!rlang::sym(predr_var_nm_1L_chr), 
            y = !!rlang::sym(dep_var_nm_1L_chr)), col = "red") + 
        ggplot2::theme_bw() + ggplot2::labs(x = predr_var_desc_1L_chr, 
        y = dep_var_desc_1L_chr)
}
#' Plot lnr cmprsn sctr plt
#' @description plot_lnr_cmprsn_sctr_plt() is a Plot function that plots data Specifically, this function implements an algorithm to plot lnr cmprsn sctr plt. The function is called for its side effects and does not return a value.
#' @param tfd_data_tb Transformed data (a tibble)
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'aqol6d_total_w'
#' @param predd_val_var_nm_1L_chr Predd value var name (a character vector of length one), Default: 'Predicted'
#' @return NULL
#' @rdname plot_lnr_cmprsn_sctr_plt
#' @export 
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth theme_bw geom_abline xlim ylim
#' @importFrom rlang sym
plot_lnr_cmprsn_sctr_plt <- function (tfd_data_tb, dep_var_nm_1L_chr = "aqol6d_total_w", 
    predd_val_var_nm_1L_chr = "Predicted") 
{
    tfd_data_tb %>% ggplot2::ggplot(ggplot2::aes(x = !!rlang::sym(dep_var_nm_1L_chr), 
        y = !!rlang::sym(predd_val_var_nm_1L_chr))) + ggplot2::geom_point() + 
        ggplot2::geom_smooth(method = "loess", size = 1.5) + 
        ggplot2::theme_bw() + ggplot2::geom_abline(intercept = 0, 
        slope = 1, color = "red", linetype = "dashed", size = 1.5) + 
        ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1)
}
#' Plot obsd predd dnst
#' @description plot_obsd_predd_dnst() is a Plot function that plots data Specifically, this function implements an algorithm to plot obsd predd dnst. The function is called for its side effects and does not return a value.
#' @param tfd_data_tb Transformed data (a tibble)
#' @param dep_var_desc_1L_chr Dep var description (a character vector of length one), Default: 'AQoL-6D utility score'
#' @param predd_val_var_nm_1L_chr Predd value var name (a character vector of length one), Default: 'Predicted'
#' @return NULL
#' @rdname plot_obsd_predd_dnst
#' @export 
#' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @importFrom rlang sym
#' @importFrom ggplot2 ggplot aes geom_rug theme_bw theme labs
#' @importFrom ggalt geom_bkde
#' @importFrom viridis scale_fill_viridis
plot_obsd_predd_dnst <- function (tfd_data_tb, dep_var_desc_1L_chr = "AQoL-6D utility score", 
    predd_val_var_nm_1L_chr = "Predicted") 
{
    tfd_data_tb %>% dplyr::mutate(Observed = aqol6d_total_w) %>% 
        tidyr::gather(variable, value, !!rlang::sym(predd_val_var_nm_1L_chr), 
            Observed) %>% ggplot2::ggplot(ggplot2::aes(x = value, 
        fill = variable)) + ggalt::geom_bkde(alpha = 0.5) + ggplot2::geom_rug() + 
        viridis::scale_fill_viridis(discrete = TRUE) + ggplot2::theme_bw() + 
        ggplot2::theme(legend.position = "bottom") + ggplot2::labs(x = dep_var_desc_1L_chr, 
        y = "Density", fill = "")
}
#' Plot obsd predd sctr
#' @description plot_obsd_predd_sctr() is a Plot function that plots data Specifically, this function implements an algorithm to plot obsd predd sctr. The function is called for its side effects and does not return a value.
#' @param tfd_data_tb Transformed data (a tibble)
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'aqol6d_total_w'
#' @param dep_var_desc_1L_chr Dep var description (a character vector of length one), Default: 'AQoL-6D utility score'
#' @param round_var_nm_1L_chr Round var name (a character vector of length one), Default: 'round'
#' @param args_ls Arguments (a list), Default: NULL
#' @param predd_val_var_nm_1L_chr Predd value var name (a character vector of length one), Default: 'Predicted'
#' @return NULL
#' @rdname plot_obsd_predd_sctr
#' @export 
#' @importFrom ggplot2 ggplot geom_point aes theme_bw xlim ylim scale_color_manual labs theme
#' @importFrom rlang exec sym
plot_obsd_predd_sctr <- function (tfd_data_tb, dep_var_nm_1L_chr = "aqol6d_total_w", 
    dep_var_desc_1L_chr = "AQoL-6D utility score", round_var_nm_1L_chr = "round", 
    args_ls = NULL, predd_val_var_nm_1L_chr = "Predicted") 
{
    ggplot2::ggplot(tfd_data_tb) + rlang::exec(ggplot2::geom_point, 
        ggplot2::aes(x = !!rlang::sym(dep_var_nm_1L_chr), y = !!rlang::sym(predd_val_var_nm_1L_chr), 
            col = !!rlang::sym(round_var_nm_1L_chr)), size = 1, 
        !!!args_ls) + ggplot2::theme_bw() + ggplot2::xlim(0, 
        1) + ggplot2::ylim(0, 1) + ggplot2::scale_color_manual(values = c("#D55E00", 
        "#56B4E9")) + ggplot2::labs(x = paste0("Observed ", dep_var_desc_1L_chr), 
        y = paste0(predd_val_var_nm_1L_chr, " ", dep_var_desc_1L_chr), 
        col = "") + ggplot2::theme(legend.position = "bottom")
}
