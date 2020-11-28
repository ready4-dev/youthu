randomise_changes_in_fct_levs <- function (vector_fct, prob_unchanged_dbl) 
{
    labels_chr <- levels(vector_fct)
    levels_dbl <- 1:length(labels_chr)
    unchanged_lgl <- runif(length(vector_fct)) > prob_unchanged_dbl
    vector_fct %>% purrr::map2_dbl(unchanged_lgl, ~{
        ifelse(.y, .x, sample(levels_dbl[levels_dbl != .x], 1))
    }) %>% factor(levels = levels_dbl, labels = labels_chr)
}
