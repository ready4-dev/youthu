knit_mdl_rprt <- function (knit_pars_ls, path_to_mdl_rprt_tmpl_1L_chr = system.file("_Model_Report_Template.Rmd", 
    package = "ready4u")) 
{
    src <- purrr::pmap(knit_pars_ls, ~knitr::knit_expand(path_to_mdl_rprt_tmpl_1L_chr, 
        plt_nms_chr = ..1 %>% deparse(), path_to_mdl_1L_chr = ..2 %>% 
            deparse(), caption_1L_chr = ..3 %>% deparse(), label_stub_1L_chr = ..4 %>% 
            deparse(), output_type_1L_chr = ..5 %>% deparse(), 
        section_ttl_1L_chr = ..6))
    res <- knitr::knit_child(text = unlist(src), quiet = TRUE)
    cat(res, sep = "\n")
}
