
#' First Bounce S3 class for Kessler Psychological Distress Scale (K6) - US Scoring System scores
#' @description Create a new valid instance of the First Bounce S3 class for Kessler Psychological Distress Scale (K6) - US Scoring System scores
#' @param x A prototype for the First Bounce S3 class for Kessler Psychological Distress Scale (K6) - US Scoring System scores, Default: make_prototype_firstbounce_k6()
#' @return A validated instance of the First Bounce S3 class for Kessler Psychological Distress Scale (K6) - US Scoring System scores
#' @details First Bounce S3 class for Kessler Psychological Distress Scale (K6) - US Scoring System scores
#' @rdname firstbounce_k6
#' @export 

firstbounce_k6 <- function(x = make_prototype_firstbounce_k6()){ 
validate_firstbounce_k6(make_new_firstbounce_k6(x))
}
#' Make new First Bounce S3 class for Kessler Psychological Distress Scale (K6) - US Scoring System scores
#' @description Create a new unvalidated instance of the First Bounce S3 class for Kessler Psychological Distress Scale (K6) - US Scoring System scores
#' @param x A prototype for the First Bounce S3 class for Kessler Psychological Distress Scale (K6) - US Scoring System scores
#' @return An unvalidated instance of the First Bounce S3 class for Kessler Psychological Distress Scale (K6) - US Scoring System scores
#' @details First Bounce S3 class for Kessler Psychological Distress Scale (K6) - US Scoring System scores
#' @rdname make_new_firstbounce_k6
#' @export 

make_new_firstbounce_k6 <- function(x){ 
stopifnot(is.integer(x))
class(x) <- append(c("firstbounce_k6",setdiff(make_prototype_firstbounce_k6() %>% class(),class(x))),
class(x))
x
}
#' Make prototype First Bounce S3 class for Kessler Psychological Distress Scale (K6) - US Scoring System scores
#' @description Create a new prototype for the First Bounce S3 class for Kessler Psychological Distress Scale (K6) - US Scoring System scores

#' @return A prototype for First Bounce S3 class for Kessler Psychological Distress Scale (K6) - US Scoring System scores
#' @details First Bounce S3 class for Kessler Psychological Distress Scale (K6) - US Scoring System scores
#' @rdname make_prototype_firstbounce_k6
#' @export 

make_prototype_firstbounce_k6 <- function(){ 
integer(0)
}
#' Validate First Bounce S3 class for Kessler Psychological Distress Scale (K6) - US Scoring System scores
#' @description Validate an instance of the First Bounce S3 class for Kessler Psychological Distress Scale (K6) - US Scoring System scores
#' @param x An unvalidated instance of the First Bounce S3 class for Kessler Psychological Distress Scale (K6) - US Scoring System scores
#' @return A prototpe for First Bounce S3 class for Kessler Psychological Distress Scale (K6) - US Scoring System scores
#' @details First Bounce S3 class for Kessler Psychological Distress Scale (K6) - US Scoring System scores
#' @rdname validate_firstbounce_k6
#' @export 

validate_firstbounce_k6 <- function(x){
if(any(x < 0)){
stop("All values in valid firstbounce_k6 object must be greater than or equal to 0.",
call. = FALSE)
}
 if(any(x > 24)){
stop("All values in valid firstbounce_k6 object must be less than or equal to 24.",
call. = FALSE)
}
x}
#' Is First Bounce S3 class for Kessler Psychological Distress Scale (K6) - US Scoring System scores
#' @description Check whether an object is a valid instance of the First Bounce S3 class for Kessler Psychological Distress Scale (K6) - US Scoring System scores
#' @param x An object of any type
#' @return A logical value, TRUE if a valid instance of the First Bounce S3 class for Kessler Psychological Distress Scale (K6) - US Scoring System scores
#' @details First Bounce S3 class for Kessler Psychological Distress Scale (K6) - US Scoring System scores
#' @rdname is_firstbounce_k6
#' @export 

is_firstbounce_k6 <- function(x) inherits(validate_firstbounce_k6(x), "firstbounce_k6")
