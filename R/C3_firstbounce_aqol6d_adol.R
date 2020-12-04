
#' First Bounce S3 class for Assessment of Quality of Life Six Dimension Health Utility - Adolescent Version (AQoL6d Adolescent))
#' @description Create a new valid instance of the First Bounce S3 class for Assessment of Quality of Life Six Dimension Health Utility - Adolescent Version (AQoL6d Adolescent))
#' @param x A prototype for the First Bounce S3 class for Assessment of Quality of Life Six Dimension Health Utility - Adolescent Version (AQoL6d Adolescent)), Default: make_pt_firstbounce_aqol6d_adol()
#' @return A validated instance of the First Bounce S3 class for Assessment of Quality of Life Six Dimension Health Utility - Adolescent Version (AQoL6d Adolescent))
#' @details First Bounce S3 class for Assessment of Quality of Life Six Dimension Health Utility - Adolescent Version (AQoL6d Adolescent))
#' @rdname firstbounce_aqol6d_adol
#' @export 

firstbounce_aqol6d_adol <- function(x = make_pt_firstbounce_aqol6d_adol()){ 
validate_firstbounce_aqol6d_adol(make_new_firstbounce_aqol6d_adol(x))
}
#' Make new First Bounce S3 class for Assessment of Quality of Life Six Dimension Health Utility - Adolescent Version (AQoL6d Adolescent))
#' @description Create a new unvalidated instance of the First Bounce S3 class for Assessment of Quality of Life Six Dimension Health Utility - Adolescent Version (AQoL6d Adolescent))
#' @param x A prototype for the First Bounce S3 class for Assessment of Quality of Life Six Dimension Health Utility - Adolescent Version (AQoL6d Adolescent))
#' @return An unvalidated instance of the First Bounce S3 class for Assessment of Quality of Life Six Dimension Health Utility - Adolescent Version (AQoL6d Adolescent))
#' @details First Bounce S3 class for Assessment of Quality of Life Six Dimension Health Utility - Adolescent Version (AQoL6d Adolescent))
#' @rdname make_new_firstbounce_aqol6d_adol
#' @export 

make_new_firstbounce_aqol6d_adol <- function(x){ 
stopifnot(is.numeric(x))
class(x) <- append(c("firstbounce_aqol6d_adol",setdiff(make_pt_firstbounce_aqol6d_adol() %>% class(),class(x))),
class(x))
x
}
#' Make prototype First Bounce S3 class for Assessment of Quality of Life Six Dimension Health Utility - Adolescent Version (AQoL6d Adolescent))
#' @description Create a new prototype for the First Bounce S3 class for Assessment of Quality of Life Six Dimension Health Utility - Adolescent Version (AQoL6d Adolescent))

#' @return A prototype for First Bounce S3 class for Assessment of Quality of Life Six Dimension Health Utility - Adolescent Version (AQoL6d Adolescent))
#' @details First Bounce S3 class for Assessment of Quality of Life Six Dimension Health Utility - Adolescent Version (AQoL6d Adolescent))
#' @rdname make_pt_firstbounce_aqol6d_adol
#' @export 

make_pt_firstbounce_aqol6d_adol <- function(){ 
numeric(0)
}
#' Validate First Bounce S3 class for Assessment of Quality of Life Six Dimension Health Utility - Adolescent Version (AQoL6d Adolescent))
#' @description Validate an instance of the First Bounce S3 class for Assessment of Quality of Life Six Dimension Health Utility - Adolescent Version (AQoL6d Adolescent))
#' @param x An unvalidated instance of the First Bounce S3 class for Assessment of Quality of Life Six Dimension Health Utility - Adolescent Version (AQoL6d Adolescent))
#' @return A prototpe for First Bounce S3 class for Assessment of Quality of Life Six Dimension Health Utility - Adolescent Version (AQoL6d Adolescent))
#' @details First Bounce S3 class for Assessment of Quality of Life Six Dimension Health Utility - Adolescent Version (AQoL6d Adolescent))
#' @rdname validate_firstbounce_aqol6d_adol
#' @export 

validate_firstbounce_aqol6d_adol <- function(x){
if(any(x < 0.03)){
stop("All values in valid firstbounce_aqol6d_adol object must be greater than or equal to 0.03.",
call. = FALSE)
}
 if(any(x > 1)){
stop("All values in valid firstbounce_aqol6d_adol object must be less than or equal to 1.",
call. = FALSE)
}
x}
#' Is First Bounce S3 class for Assessment of Quality of Life Six Dimension Health Utility - Adolescent Version (AQoL6d Adolescent))
#' @description Check whether an object is a valid instance of the First Bounce S3 class for Assessment of Quality of Life Six Dimension Health Utility - Adolescent Version (AQoL6d Adolescent))
#' @param x An object of any type
#' @return A logical value, TRUE if a valid instance of the First Bounce S3 class for Assessment of Quality of Life Six Dimension Health Utility - Adolescent Version (AQoL6d Adolescent))
#' @details First Bounce S3 class for Assessment of Quality of Life Six Dimension Health Utility - Adolescent Version (AQoL6d Adolescent))
#' @rdname is_firstbounce_aqol6d_adol
#' @export 

is_firstbounce_aqol6d_adol <- function(x) inherits(validate_firstbounce_aqol6d_adol(x), "firstbounce_aqol6d_adol")
