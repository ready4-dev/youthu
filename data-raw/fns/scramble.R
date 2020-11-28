scramble_xx <- function (vector_xx) 
{
    scrambled_vec_xx <- vector_xx[sample(1:length(vector_xx))]
    return(scrambled_vec_xx)
}
