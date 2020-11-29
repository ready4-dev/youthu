
#' First Bounce S3 class for Patient Health Questionnaire (PHQ-9) scores
#' @description Create a new valid instance of the First Bounce S3 class for Patient Health Questionnaire (PHQ-9) scores
#' @param x A prototype for the First Bounce S3 class for Patient Health Questionnaire (PHQ-9) scores, Default: make_prototype_firstbounce_phq9()
#' @return A validated instance of the First Bounce S3 class for Patient Health Questionnaire (PHQ-9) scores
#' @details First Bounce S3 class for Patient Health Questionnaire (PHQ-9) scores
#' @rdname firstbounce_phq9
#' @export 

firstbounce_phq9 <- function(x = make_prototype_firstbounce_phq9()){ 
validate_firstbounce_phq9(make_new_firstbounce_phq9(x))
}
#' Make new First Bounce S3 class for Patient Health Questionnaire (PHQ-9) scores
#' @description Create a new unvalidated instance of the First Bounce S3 class for Patient Health Questionnaire (PHQ-9) scores
#' @param x A prototype for the First Bounce S3 class for Patient Health Questionnaire (PHQ-9) scores
#' @return An unvalidated instance of the First Bounce S3 class for Patient Health Questionnaire (PHQ-9) scores
#' @details First Bounce S3 class for Patient Health Questionnaire (PHQ-9) scores
#' @rdname make_new_firstbounce_phq9
#' @export 

make_new_firstbounce_phq9 <- function(x){ 
stopifnot(is.integer(x))
class(x) <- append(c("firstbounce_phq9",setdiff(make_prototype_firstbounce_phq9() %>% class(),class(x))),
class(x))
x
}
#' Make prototype First Bounce S3 class for Patient Health Questionnaire (PHQ-9) scores
#' @description Create a new prototype for the First Bounce S3 class for Patient Health Questionnaire (PHQ-9) scores

#' @return A prototype for First Bounce S3 class for Patient Health Questionnaire (PHQ-9) scores
#' @details First Bounce S3 class for Patient Health Questionnaire (PHQ-9) scores
#' @rdname make_prototype_firstbounce_phq9
#' @export 

make_prototype_firstbounce_phq9 <- function(){ 
integer(0)
}
#' Validate First Bounce S3 class for Patient Health Questionnaire (PHQ-9) scores
#' @description Validate an instance of the First Bounce S3 class for Patient Health Questionnaire (PHQ-9) scores
#' @param x An unvalidated instance of the First Bounce S3 class for Patient Health Questionnaire (PHQ-9) scores
#' @return A prototpe for First Bounce S3 class for Patient Health Questionnaire (PHQ-9) scores
#' @details First Bounce S3 class for Patient Health Questionnaire (PHQ-9) scores
#' @rdname validate_firstbounce_phq9
#' @export 

validate_firstbounce_phq9 <- function(x){
if(any(x < 0)){
stop("All values in valid firstbounce_phq9 object must be greater than or equal to 0.",
call. = FALSE)
}
 if(any(x > 27)){
stop("All values in valid firstbounce_phq9 object must be less than or equal to 27.",
call. = FALSE)
}
x}
#' Is First Bounce S3 class for Patient Health Questionnaire (PHQ-9) scores
#' @description Check whether an object is a valid instance of the First Bounce S3 class for Patient Health Questionnaire (PHQ-9) scores
#' @param x An object of any type
#' @return A logical value, TRUE if a valid instance of the First Bounce S3 class for Patient Health Questionnaire (PHQ-9) scores
#' @details First Bounce S3 class for Patient Health Questionnaire (PHQ-9) scores
#' @rdname is_firstbounce_phq9
#' @export 

is_firstbounce_phq9 <- function(x) inherits(validate_firstbounce_phq9(x), "firstbounce_phq9")
