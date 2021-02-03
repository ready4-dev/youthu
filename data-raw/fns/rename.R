rename_from_nmd_vec <- function(data_tb,
                                nmd_vec_chr,
                                vec_nms_as_new_1L_lgl = T){
  if(vec_nms_as_new_1L_lgl){
    renamed_data_tb <- purrr::reduce(1:length(nmd_vec_chr),
                                     .init = data_tb,
                                     ~ dplyr::rename(.x,
                                                     !!rlang::sym(names(nmd_vec_chr)[.y]) := nmd_vec_chr[.y] %>% as.vector()))
  }else{
    renamed_data_tb <- purrr::reduce(1:length(nmd_vec_chr),
                                     .init = data_tb,
                                     ~ dplyr::rename(.x,
                                                     !!rlang::sym(nmd_vec_chr[.y]) := names(nmd_vec_chr)[.y] %>% as.vector()))
  }
  return(renamed_data_tb)
}