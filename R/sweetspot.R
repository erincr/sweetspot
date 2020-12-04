#' Perform a sweet spot analysis on clinical trial data. 
#' 
#' @description Identifying heterogeneous treatment effects (HTEs) in randomized controlled trials is an important step toward understanding and acting on trial results. However, HTEs are often small and difficult to identify, and HTE modeling methods which are very general can suffer from low power. This method exploits any existing relationship between illness severity and treatment effect, and identifies the "sweet spot", the contiguous range of illness severity where the estimated treatment benefit is maximized. We further compute a bias-corrected estimate of the conditional average treatment effect (CATE) in the sweet spot, and a p-value. More information here: \url{https://arxiv.org/abs/2011.10157}. 
#' 
#' @param treated A binary vector with one entry for every individual in the trial. The value 1 indicates that the individual was treated, 
#'  and the value 0 indicates that they were a control.
#' @param covariates A numeric matrix containing rows of covariates for every individual in the trial (in the same order as in `treated`).
#' @param negative.outcome The response we use to compute our predilection score. This will determine the model we use.
#' This is a vector with values that may be binary (logistic regression) or continuous (linear regression). 
#' In the binary case, this may be 1 if the patient died, and 0 else.
#' @param positive.outcome The response we use to compute the treatment effect. In the binary case, this may be 1 if the patient lived, and 0 else.
#' @param family A string indicating the response type. Options are "binomial" (for logistic regression) or "gaussian" (for linear regression).
#' @param regularized Boolean indicating whether to use regularization in the predilection score model. If TRUE, we use the optimal model returned by cv.glmnet.
#' @param control.treated.ratio The (positive integer) number of controls for each treated individual.
#' @param predilection.score.nfolds The (positive integer) number of folds we use in pre-validation for the predilection score. The default is 10 folds.
#' @param ntrials.significance The number of bootstraps to run when computing the p-value. Default is 1000.
#' @param ntrials.bias The number of bootstraps to run when debiasing the estimate of the treatment effect. Default is 1000.
#' @param ncores The number of cores to use when finding the sweet spot. If NULL, defaults to parallel::detectCores()-1.
#' @return A list of results, containing: \itemize{
#'         \item{"matches", the matrix of matched treated and control patients, their mean predilection score and treatment effect estimate}
#'         \item{"predilection.score.model", the fitted glmnet object that models the predilection score}
#'         \item{"predilection.scores", the vector of prevalidated predilection scores for each patient}
#'         \item{"model", the sweet spot model, a list containing:\itemize{
#'             \item{the start and end index of the sweet spot}
#'             \item{the bootstrapped p-value}
#'             \item{the mean inside and outside the sweet spot, before and after debiasing}
#'             \item{the distribution around the start and end indices, from the debiasing bootstrap}
#'             \item{the maximum statistic used to compute the p-value}
#'             }}
#' }
#' 
#' @include match_individuals.R
#' @include predilection_scores.R
#' @include find_sweetspot.R
#' @include plot_result.R
#' @include aquamat-data.R
#' 
#' @importFrom glmnet glmnet 
#' @importFrom optmatch pairmatch
#' @import doParallel
#' @importFrom foreach %dopar% foreach
#' @importFrom stats smooth.spline predict
#' @importFrom parallel detectCores
#' 
#' @examples 
#' # Example data with a sweet spot
#' # We generate data with a sweetspot in the middle 20% of predilection scores.
#' # Outside the sweet spot, the treatment effect is 5%.
#' # Inside the sweet spot, the treatment effect is 25%.
#' set.seed(1234)
#' n <- 500; p <- 10;
#' treated    <- sample(c(0,1), n, replace=TRUE)
#' covariates <- matrix(rnorm(n * p), nrow=n, ncol=p)
#' beta       <- rnorm(p)
#' outcome.prob <- 1/(1+exp(-(covariates %*% beta)))
#' in.sweet.spot <- !is.na(cut(outcome.prob, c(.4, .6)))
#' outcome.prob[treated==1] <- outcome.prob[treated==1] + .05
#' outcome.prob[treated==1 & in.sweet.spot] <- outcome.prob[treated==1 & in.sweet.spot] + .2
#' outcome.prob     <- pmin(outcome.prob, 1)
#' positive.outcome <- rbinom(n, 1, prob=outcome.prob)
#' negative.outcome <- 1-positive.outcome
#' 
#' # We limit the number of trials run only to illustrate how to use this function.
#' # In practice, we recommend more bootstrap trials.
#' result <- sweetspot(treated, covariates, 
#'                     negative.outcome, positive.outcome, 
#'                     "binomial", ntrials.bias=250, ntrials.significance=250)
#' plot_sweetspot(result, title="Sweet spot on simulated data")
#' 
#' 
#' # Example data without a sweet spot. 
#' # The overall treatment effect is 5%.
#' set.seed(1234)
#' n <- 500; p <- 10;
#' treated    <- sample(c(0,1), n, replace=TRUE)
#' covariates <- matrix(rnorm(n * p), nrow=n, ncol=p)
#' beta       <- rnorm(p)
#' outcome.probability <- 1/(1+exp(-(covariates %*% beta)))
#' outcome.probability[treated == 1] <- pmin(outcome.probability[treated == 1] + .05, 1) 
#' positive.outcome <- rbinom(n, 1, prob=outcome.probability)
#' negative.outcome <- 1-positive.outcome
#' 
#' result <- sweetspot(treated, covariates, 
#'                     negative.outcome, positive.outcome, 
#'                     "binomial", ntrials.bias=250, ntrials.significance=250)
#' plot_sweetspot(result, title="Sweet spot on simulated data")
#' 
#' @export sweetspot

sweetspot <- function(treated, covariates, negative.outcome, positive.outcome, family, 
                      regularized = F, control.treated.ratio = 1, predilection.score.nfolds = 10,
                      ntrials.significance=1000, ntrials.bias=1000, ncores=NULL){
  
  # check types
  if(!((class(covariates)[1] == "matrix") & (typeof(covariates) == "double"))) { 
    message("The covariates are not a numeric matrix."); stop()
  } 
  
  if(!(family %in% c("binomial", "gaussian"))){ message("family should be 'binomial' or 'gaussian'."); stop()}
   
  if(!is.null(ncores)){
    if(!is.numeric(ncores)) { message("Number of cores should be an integer or NULL."); stop() }
    if(floor(ncores) != ncores)  { message("Number of cores should be an integer or NULL."); stop() }
  }
  if(!is.numeric(c(ntrials.significance, ntrials.bias))){ message("Number of trials should be integers."); stop() }
  if(floor(ntrials.significance) != ntrials.significance) { message("ntrials.significance should be an integer."); stop() }
  if(floor(ntrials.bias) != ntrials.bias) { message("ntrials.bias should be an integer."); stop() }
  
  # compute the predilection scores using prevalidation:
  predilection.scores <- predilection_scores(treated, 
                                             covariates, 
                                             negative.outcome,
                                             family = family, 
                                             regularized = regularized, 
                                             nfolds  = predilection.score.nfolds)
  predilection.score.model = predilection.scores$model
  predilection.scores      = predilection.scores$scores
  
  # match individuals using predilection.scores:
  matches <- match.individuals(treated, predilection.scores, control.treated.ratio, 
                                           positive.outcome, function(treated, control) treated - mean(control))
  
  # find the sweet spot, estimate the p-value, and debias the CATE.
  model <- find.sweetspot(matches[, "apparent_benefit"], 
                          ntrials.significance = ntrials.significance, 
                          ntrials.bias = ntrials.bias,
                          ncores = ncores)

  return(list(
              matches = matches,
              predilection.score.model = predilection.score.model,
              predilection.scores = predilection.scores,
              model = model
  ))
}
