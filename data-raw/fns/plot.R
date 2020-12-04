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
