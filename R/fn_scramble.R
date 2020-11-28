#' Scramble
#' @description scramble_xx() is a Scramble function that randomly reorders an object. Specifically, this function implements an algorithm to scramble output object of multiple potential types. The function returns Scrambled vec (an output object of multiple potential types).
#' @param vector_xx Vector (an output object of multiple potential types)
#' @return Scrambled vec (an output object of multiple potential types)
#' @rdname scramble_xx
#' @export 

scramble_xx <- function (vector_xx) 
{
    scrambled_vec_xx <- vector_xx[sample(1:length(vector_xx))]
    return(scrambled_vec_xx)
}
