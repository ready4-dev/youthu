
#' First Bounce S3 class for Screen for Child Anxiety Related Disorders (SCARED) scores
#' @description Create a new valid instance of the First Bounce S3 class for Screen for Child Anxiety Related Disorders (SCARED) scores
#' @param x A prototype for the First Bounce S3 class for Screen for Child Anxiety Related Disorders (SCARED) scores, Default: make_prototype_firstbounce_scared()
#' @return A validated instance of the First Bounce S3 class for Screen for Child Anxiety Related Disorders (SCARED) scores
#' @details First Bounce S3 class for Screen for Child Anxiety Related Disorders (SCARED) scores
#' @rdname firstbounce_scared
#' @export 

firstbounce_scared <- function(x = make_prototype_firstbounce_scared()){ 
validate_firstbounce_scared(make_new_firstbounce_scared(x))
}
#' Make new First Bounce S3 class for Screen for Child Anxiety Related Disorders (SCARED) scores
#' @description Create a new unvalidated instance of the First Bounce S3 class for Screen for Child Anxiety Related Disorders (SCARED) scores
#' @param x A prototype for the First Bounce S3 class for Screen for Child Anxiety Related Disorders (SCARED) scores
#' @return An unvalidated instance of the First Bounce S3 class for Screen for Child Anxiety Related Disorders (SCARED) scores
#' @details First Bounce S3 class for Screen for Child Anxiety Related Disorders (SCARED) scores
#' @rdname make_new_firstbounce_scared
#' @export 

make_new_firstbounce_scared <- function(x){ 
stopifnot(is.integer(x))
class(x) <- append(c("firstbounce_scared",setdiff(make_prototype_firstbounce_scared() %>% class(),class(x))),
class(x))
x
}
#' Make prototype First Bounce S3 class for Screen for Child Anxiety Related Disorders (SCARED) scores
#' @description Create a new prototype for the First Bounce S3 class for Screen for Child Anxiety Related Disorders (SCARED) scores

#' @return A prototype for First Bounce S3 class for Screen for Child Anxiety Related Disorders (SCARED) scores
#' @details First Bounce S3 class for Screen for Child Anxiety Related Disorders (SCARED) scores
#' @rdname make_prototype_firstbounce_scared
#' @export 

make_prototype_firstbounce_scared <- function(){ 
integer(0)
}
#' Validate First Bounce S3 class for Screen for Child Anxiety Related Disorders (SCARED) scores
#' @description Validate an instance of the First Bounce S3 class for Screen for Child Anxiety Related Disorders (SCARED) scores
#' @param x An unvalidated instance of the First Bounce S3 class for Screen for Child Anxiety Related Disorders (SCARED) scores
#' @return A prototpe for First Bounce S3 class for Screen for Child Anxiety Related Disorders (SCARED) scores
#' @details First Bounce S3 class for Screen for Child Anxiety Related Disorders (SCARED) scores
#' @rdname validate_firstbounce_scared
#' @export 

validate_firstbounce_scared <- function(x){
if(any(x < 0)){
stop("All values in valid firstbounce_scared object must be greater than or equal to 0.",
call. = FALSE)
}
 if(any(x > 82)){
stop("All values in valid firstbounce_scared object must be less than or equal to 82.",
call. = FALSE)
}
x}
#' Is First Bounce S3 class for Screen for Child Anxiety Related Disorders (SCARED) scores
#' @description Check whether an object is a valid instance of the First Bounce S3 class for Screen for Child Anxiety Related Disorders (SCARED) scores
#' @param x An object of any type
#' @return A logical value, TRUE if a valid instance of the First Bounce S3 class for Screen for Child Anxiety Related Disorders (SCARED) scores
#' @details First Bounce S3 class for Screen for Child Anxiety Related Disorders (SCARED) scores
#' @rdname is_firstbounce_scared
#' @export 

is_firstbounce_scared <- function(x) inherits(validate_firstbounce_scared(x), "firstbounce_scared")
