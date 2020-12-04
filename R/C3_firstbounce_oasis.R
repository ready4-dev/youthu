
#' First Bounce S3 class for Overall Anxiety Severity and Impairment Scale (OASIS) scores
#' @description Create a new valid instance of the First Bounce S3 class for Overall Anxiety Severity and Impairment Scale (OASIS) scores
#' @param x A prototype for the First Bounce S3 class for Overall Anxiety Severity and Impairment Scale (OASIS) scores, Default: make_pt_firstbounce_oasis()
#' @return A validated instance of the First Bounce S3 class for Overall Anxiety Severity and Impairment Scale (OASIS) scores
#' @details First Bounce S3 class for Overall Anxiety Severity and Impairment Scale (OASIS) scores
#' @rdname firstbounce_oasis
#' @export 

firstbounce_oasis <- function(x = make_pt_firstbounce_oasis()){ 
validate_firstbounce_oasis(make_new_firstbounce_oasis(x))
}
#' Make new First Bounce S3 class for Overall Anxiety Severity and Impairment Scale (OASIS) scores
#' @description Create a new unvalidated instance of the First Bounce S3 class for Overall Anxiety Severity and Impairment Scale (OASIS) scores
#' @param x A prototype for the First Bounce S3 class for Overall Anxiety Severity and Impairment Scale (OASIS) scores
#' @return An unvalidated instance of the First Bounce S3 class for Overall Anxiety Severity and Impairment Scale (OASIS) scores
#' @details First Bounce S3 class for Overall Anxiety Severity and Impairment Scale (OASIS) scores
#' @rdname make_new_firstbounce_oasis
#' @export 

make_new_firstbounce_oasis <- function(x){ 
stopifnot(is.integer(x))
class(x) <- append(c("firstbounce_oasis",setdiff(make_pt_firstbounce_oasis() %>% class(),class(x))),
class(x))
x
}
#' Make prototype First Bounce S3 class for Overall Anxiety Severity and Impairment Scale (OASIS) scores
#' @description Create a new prototype for the First Bounce S3 class for Overall Anxiety Severity and Impairment Scale (OASIS) scores

#' @return A prototype for First Bounce S3 class for Overall Anxiety Severity and Impairment Scale (OASIS) scores
#' @details First Bounce S3 class for Overall Anxiety Severity and Impairment Scale (OASIS) scores
#' @rdname make_pt_firstbounce_oasis
#' @export 

make_pt_firstbounce_oasis <- function(){ 
integer(0)
}
#' Validate First Bounce S3 class for Overall Anxiety Severity and Impairment Scale (OASIS) scores
#' @description Validate an instance of the First Bounce S3 class for Overall Anxiety Severity and Impairment Scale (OASIS) scores
#' @param x An unvalidated instance of the First Bounce S3 class for Overall Anxiety Severity and Impairment Scale (OASIS) scores
#' @return A prototpe for First Bounce S3 class for Overall Anxiety Severity and Impairment Scale (OASIS) scores
#' @details First Bounce S3 class for Overall Anxiety Severity and Impairment Scale (OASIS) scores
#' @rdname validate_firstbounce_oasis
#' @export 

validate_firstbounce_oasis <- function(x){
if(any(x < 0)){
stop("All values in valid firstbounce_oasis object must be greater than or equal to 0.",
call. = FALSE)
}
 if(any(x > 20)){
stop("All values in valid firstbounce_oasis object must be less than or equal to 20.",
call. = FALSE)
}
x}
#' Is First Bounce S3 class for Overall Anxiety Severity and Impairment Scale (OASIS) scores
#' @description Check whether an object is a valid instance of the First Bounce S3 class for Overall Anxiety Severity and Impairment Scale (OASIS) scores
#' @param x An object of any type
#' @return A logical value, TRUE if a valid instance of the First Bounce S3 class for Overall Anxiety Severity and Impairment Scale (OASIS) scores
#' @details First Bounce S3 class for Overall Anxiety Severity and Impairment Scale (OASIS) scores
#' @rdname is_firstbounce_oasis
#' @export 

is_firstbounce_oasis <- function(x) inherits(validate_firstbounce_oasis(x), "firstbounce_oasis")
