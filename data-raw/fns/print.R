print_all_plts_for_mdl_set <- function (output_ls, start_from_1L_int = 0L) 
{
    label_refs_mat <- paste0("lab", 1:(output_ls %>% purrr::flatten() %>% 
        length()) + start_from_1L_int) %>% matrix(nrow = length(output_ls), 
        byrow = T)
    purrr::walk(1:length(output_ls), ~print_ts_mdl_plts(output_ls[[.x]], 
        title_1L_chr = names(output_ls)[.x], label_refs_chr = label_refs_mat[.x, 
            ]))
}
print_table_xx <- function (data_tb, output_type_1L_chr = "PDF", caption_1L_chr = NA_character_, 
    footnotes_chr = NA_character_, merge_row_idx_int = NA_integer_, 
    digits_dbl = NULL, big_mark_1L_chr = " ", label, hline.after, 
    addtorow, sanitize_fn = getOption("xtable.sanitize.text.function", 
        NULL)) 
{
    if (output_type_1L_chr == "PDF") {
        data_x_tb <- data_tb %>% xtable::xtable(caption = caption_1L_chr, 
            label = label, digits = digits_dbl)
        data_x_tb %>% print(comment = F, floating = TRUE, hline.after = hline.after, 
            caption.placement = "top", add.to.row = addtorow, 
            sanitize.text.function = sanitize_fn, format.args = list(big.mark = big_mark_1L_chr), 
            include.colnames = FALSE, include.rownames = FALSE)
    }
    else {
        if (output_type_1L_chr == "HTML") {
            data_kb <- data_tb %>% kableExtra::kbl(format = "html", 
                escape = F, caption = caption_1L_chr, align = c("l", 
                  rep("r", nrow(data_tb) - 1))) %>% kableExtra::kable_paper("hover", 
                full_width = F)
            if (!all(is.na(merge_row_idx_int))) {
                data_kb <- data_kb %>% kableExtra::row_spec(merge_row_idx_int, 
                  bold = T)
            }
            data_kb <- data_kb %>% kableExtra::footnote(general = footnotes_chr, 
                general_title = "", footnote_as_chunk = T, threeparttable = T)
            data_xx <- data_kb
        }
        if (output_type_1L_chr == "Word") {
            j2_1L_dbl <- ncol(data_tb)
            data_fx <- flextable::flextable(data_tb) %>% flextable::set_caption(caption_1L_chr, 
                autonum = officer::run_autonum())
            if (!all(is.na(merge_row_idx_int))) {
                data_fx <- data_fx %>% flextable::merge_h_range(i = merge_row_idx_int, 
                  j1 = 1, j2 = j2_1L_dbl) %>% flextable::bold(i = merge_row_idx_int, 
                  j = 1)
            }
            if (!is.null(digits_dbl)) {
                numeric_lgl <- unlist(lapply(data_tb, is.numeric)) %>% 
                  unname()
                data_fx <- unique(digits_dbl) %>% purrr::reduce(.init = data_fx, 
                  ~.x %>% flextable::colformat_num(j = names(data_tb)[which(digits_dbl == 
                    .y & numeric_lgl)], digits = .y))
            }
            data_fx <- data_fx %>% flextable::autofit(add_w = 0, 
                add_h = 0)
            if (!all(is.na(footnotes_chr))) {
                data_fx <- data_fx %>% flextable::add_footer_lines(values = footnotes_chr)
            }
            data_xx <- data_fx
        }
        return(data_xx)
    }
}
print_ts_mdl_plts <- function (paths_to_plts_chr, title_1L_chr, label_refs_chr, mdl_smry_ls) 
{
    cat("\n")
    cat("## ", title_1L_chr, "\n")
    cat("\n")
    purrr::pwalk(list(paths_to_plts_chr, paste0("Caption", 1:length(paths_to_plts_chr)), 
        label_refs_chr), ~{
        cat(paste0("\n![", "](", ..1, ")\n\n"))
        cat(fig_nums(..3, ..2))
        cat("\n\n")
    })
    cat("\n\n")
}
