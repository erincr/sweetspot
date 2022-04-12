library(optmatch)
# library(RcppRoll)

###################################################################################################
# Matching individuals
#
# We use optmatch::pairmatch to perform optimal matching to pair controls with treated individuals.
# We take the following arguments:
# 1. treated: a binary vector with one entry for every individual in the trial. The value 1 indicates that the individual was treated, 
#  and the value 0 indicates that they were a control.
# 2. risk.scores: a numeric vector containing risk scores.
# 3. control.treated.ratio: the integer number of controls for each treated individual.
# 4. response: the response for each individual.
###################################################################################################

match.individuals <- function(treated, risk.scores, control.treated.ratio, response){
  
  df <- data.frame(trt=treated, score=risk.scores, name=1:length(treated), ix=1:length(treated))
  df[df$trt == 1, "name"] = paste0("treated_", df[df$trt == 1, "name"])
  df[df$trt == 0, "name"] = paste0("control_", df[df$trt == 0, "name"])
  

  pm <- pairmatch(trt ~ score, data=df, controls=control.treated.ratio, remove.unmatchables = T)
  pm <- pm[!is.na(pm)]
  
  # Compute information for each match (mean risk score, treatment effect)
  get_group_info <- function(group){
    gp <- df[df[, "ix"] %in% as.numeric(names(pm[pm == group])), "name"]
    
    # Index of controls and treated:
    ctrls <- as.numeric(gsub("control_", "", gp[grepl("control_", gp)]))
    trtds <- as.numeric(gsub("treated_", "", gp[grepl("treated_", gp)]))
    
    treatment.effect <- mean(response[trtds] - response[ctrls])
    
    c(ctrls, trtds,                                    # id of each control/treated
      risk.scores[ctrls],                              # risk scores
      risk.scores[trtds],
      mean(c(risk.scores[ctrls], risk.scores[trtds])), # mean score for the group
      treatment.effect
    )
  }
  all.matches <- t(sapply(unique(pm), get_group_info))

  ctrl.names <- if(control.treated.ratio == 1) {"control"} else {paste0("control_", 1:control.treated.ratio)}
  colnames(all.matches) <- c(ctrl.names, "treated", paste0(ctrl.names, "_score"), "treated_score", "mean_score", "treatment_effect")

  # return the matches in order of increasing risk score
  all.matches[order(all.matches[, "mean_score"]), ]
}
