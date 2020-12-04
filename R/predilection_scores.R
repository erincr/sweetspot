library(glmnet)

####################################################################################################################
# Compute predilection scores.
#
# Compute predilection scores on clinical trial data, using pre-validation. 
# To compute the predilection score, we use Cox regression (survival response), linera regerssion (continuous response),
# or logistic regression (binary response).
# 1. To compute predilection scores for the controls, we use pre-validation: we partition the data into `nfolds` folds, we train a model on `nfolds-1` folds, and compute the score on the remaining fold. 
#    We do this to avoid introducing bias: each individual's score is derived from other individuals in the data.
# 2. To compute predilections cores for the the treated individuals, we use a model trained a model on all controls.
# 
# 
# param treated A binary vector with one entry for every individual in the trial. The value 1 indicates that the individual was treated, 
#  and the value 0 indicates that they were a control.
# param covariates A numeric matrix containing a row of covariates for every individual in the trial (in the same order as in `treated`).
# param response The response we use to compute our predilection score. This will determine the model we use.
# This may be a binary (logistic regression) vector, a continuous (linear regression) vector, or a `Surv` object (Cox proportional hazards).
# param family The response type: one of "Cox", "binomial" or "gaussian".
# param regularized Logical flag for whether to use regularization when computing scores. Default is `FALSE`. If `TRUE`, we predict using 
# param nfolds The (positive integer) number of folds we use in pre-validation. The default is 10 folds.
# returns A vector of predilection scores, one for every individual in the data.
####################################################################################################################

predilection_scores <- function(treated, covariates, response, family = "binomial", regularized = F, nfolds = NULL){
  # Argument checking:
  # if((class(response) == "Surv") & (family != "cox")) { message("The response type (Surv) does not match the family. Use `family = \"cox\" for survival outcomes."); stop()} 
  if((class(response) == "Surv") & (family != "cox")) { message("The response type (Surv) is not binary or continuous."); stop()} 
  if((family == "binomial") & (length(unique(response)) > 2)) { message("The response type has more than two unique values. Use `family = \"gaussian\" for continuous outcomes."); stop()} 
  if(!is.null(nfolds)){
    if(nfolds > sum(treated == 0)) { message("Too many folds for prevalidating predilection scores. Make sure the number of folds is less than or equal to the number of controls."); stop() } 
    if(nfolds < 0) { message("Number of folds for prevalidation should be greater than zero."); stop() }
  }
  
  nf <- if(is.null(nfolds)) {sum(treated == 0)} else {nfolds}
  
  # Fit and predict:
  fit.model <- glmnet::cv.glmnet(covariates[treated == 0, ], 
                         response[treated == 0], 
                         family=family, 
                         nfolds=nf,
                         grouped=F,
                         # intercept=T,
                         keep=T)
  
  s         <- if(regularized) { fit.model$lambda.1se } else {  min(fit.model$lambda)  }  
  
  scores               <- rep(NA, nrow(covariates))
  scores[treated == 0] <- fit.model$fit.preval[, which(fit.model$lambda == s)]
  scores[treated == 1] <- predict(fit.model, s=s, newx=covariates[treated == 1, ])
  
  return(list(model = fit.model, scores = scores))
}
