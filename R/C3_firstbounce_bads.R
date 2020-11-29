
#' First Bounce S3 class for Behavioural Activation for Depression Scale (BADS) scores
#' @description Create a new valid instance of the First Bounce S3 class for Behavioural Activation for Depression Scale (BADS) scores
#' @param x A prototype for the First Bounce S3 class for Behavioural Activation for Depression Scale (BADS) scores, Default: make_prototype_firstbounce_bads()
#' @return A validated instance of the First Bounce S3 class for Behavioural Activation for Depression Scale (BADS) scores
#' @details First Bounce S3 class for Behavioural Activation for Depression Scale (BADS) scores
#' @rdname firstbounce_bads
#' @export 

firstbounce_bads <- function(x = make_prototype_firstbounce_bads()){ 
validate_firstbounce_bads(make_new_firstbounce_bads(x))
}
#' Make new First Bounce S3 class for Behavioural Activation for Depression Scale (BADS) scores
#' @description Create a new unvalidated instance of the First Bounce S3 class for Behavioural Activation for Depression Scale (BADS) scores
#' @param x A prototype for the First Bounce S3 class for Behavioural Activation for Depression Scale (BADS) scores
#' @return An unvalidated instance of the First Bounce S3 class for Behavioural Activation for Depression Scale (BADS) scores
#' @details First Bounce S3 class for Behavioural Activation for Depression Scale (BADS) scores
#' @rdname make_new_firstbounce_bads
#' @export 

make_new_firstbounce_bads <- function(x){ 
stopifnot(is.integer(x))
class(x) <- append(c("firstbounce_bads",setdiff(make_prototype_firstbounce_bads() %>% class(),class(x))),
class(x))
x
}
#' Make prototype First Bounce S3 class for Behavioural Activation for Depression Scale (BADS) scores
#' @description Create a new prototype for the First Bounce S3 class for Behavioural Activation for Depression Scale (BADS) scores

#' @return A prototype for First Bounce S3 class for Behavioural Activation for Depression Scale (BADS) scores
#' @details First Bounce S3 class for Behavioural Activation for Depression Scale (BADS) scores
#' @rdname make_prototype_firstbounce_bads
#' @export 

make_prototype_firstbounce_bads <- function(){ 
integer(0)
}
#' Validate First Bounce S3 class for Behavioural Activation for Depression Scale (BADS) scores
#' @description Validate an instance of the First Bounce S3 class for Behavioural Activation for Depression Scale (BADS) scores
#' @param x An unvalidated instance of the First Bounce S3 class for Behavioural Activation for Depression Scale (BADS) scores
#' @return A prototpe for First Bounce S3 class for Behavioural Activation for Depression Scale (BADS) scores
#' @details First Bounce S3 class for Behavioural Activation for Depression Scale (BADS) scores
#' @rdname validate_firstbounce_bads
#' @export 

validate_firstbounce_bads <- function(x){
if(any(x < 0)){
stop("All values in valid firstbounce_bads object must be greater than or equal to 0.",
call. = FALSE)
}
 if(any(x > 150)){
stop("All values in valid firstbounce_bads object must be less than or equal to 150.",
call. = FALSE)
}
x}
#' Is First Bounce S3 class for Behavioural Activation for Depression Scale (BADS) scores
#' @description Check whether an object is a valid instance of the First Bounce S3 class for Behavioural Activation for Depression Scale (BADS) scores
#' @param x An object of any type
#' @return A logical value, TRUE if a valid instance of the First Bounce S3 class for Behavioural Activation for Depression Scale (BADS) scores
#' @details First Bounce S3 class for Behavioural Activation for Depression Scale (BADS) scores
#' @rdname is_firstbounce_bads
#' @export 

is_firstbounce_bads <- function(x) inherits(validate_firstbounce_bads(x), "firstbounce_bads")
