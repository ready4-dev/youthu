replace_with_missing_vals <- function (data_tbl, synth_data_spine_ls, idx_int) 
{
    missing_def_tbl <- purrr::reduce(1:length(synth_data_spine_ls$var_names_chr), 
        .init = NULL, ~simstudy::defMiss(.x, varname = synth_data_spine_ls$var_names_chr[.y], 
            formula = synth_data_spine_ls$missing_ls[[idx_int]][.y]/synth_data_spine_ls$nbr_obs_dbl[idx_int], 
            logit.link = FALSE))
    missing_mat <- simstudy::genMiss(data_tbl, missing_def_tbl, 
        idvars = "id")
    synth_tbl <- simstudy::genObs(data_tbl, missing_mat, idvars = "id")
    return(synth_tbl)
}
