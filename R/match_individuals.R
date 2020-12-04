# library(optmatch)

###################################################################################################
# Matching individuals
#
# We use optmatch::pairmatch to perform optimal matching to pair controls with treated individuals.
# We take the following arguments:
# 1. treated: a binary vector with one entry for every individual in the trial. The value 1 indicates that the individual was treated, 
#  and the value 0 indicates that they were a control.
# 2. predilection.scores: a numeric vector containing predilection scores.
# 3. control.treated.ratio: the integer number of controls for each treated individual.
# 4. response: the response for each individual.
# 5. apparent.benefit.function The function to use to compute apparent benefit for each group.
###################################################################################################
match.individuals <- function(treated, predilection.scores, control.treated.ratio, response, apparent.benefit.function){
  df <- data.frame(trt=treated, score=predilection.scores, name=1:length(treated))
  df[df$trt == 1, "name"] = paste0("treated_", df[df$trt == 1, "name"])
  df[df$trt == 0, "name"] = paste0("control_", df[df$trt == 0, "name"])

  # Match treated individuals with controls, and remove those we could not match
  pm <- pairmatch(trt ~ score, data=df, controls=control.treated.ratio, remove.unmatchables=T)
  pm <- pm[!is.na(pm)]

  # Compute information for each match (mean predilection score, apparent benefit)
  get_group_info <- function(group){
    gp <- df$name[as.numeric(names(pm[pm == group]))]

    # Index of controls and treated:
    ctrls <- as.numeric(gsub("control_", "", gp[grepl("control_", gp)]))
    trtds <- as.numeric(gsub("treated_", "", gp[grepl("treated_", gp)]))
    
    c(ctrls, trtds,                                                          # id of each control/treated
      predilection.scores[ctrls],                                            # predilection score
      predilection.scores[trtds],
      mean(c(predilection.scores[ctrls], predilection.scores[trtds])),                # mean score for the group
      apparent.benefit.function(treated = response[trtds], control = response[ctrls]) # apparent benefit for the group
    )
  }
  all.matches <- t(sapply(unique(pm), get_group_info))

  ctrl.names <- if(control.treated.ratio == 1) {"control"} else {paste0("control_", 1:control.treated.ratio)}
  colnames(all.matches) <- c(ctrl.names, "treated", paste0(ctrl.names, "_score"), "treated_score", "mean_score", "apparent_benefit")

  # return the matches in order of increasing predilection score
  all.matches[order(all.matches[, "mean_score"]), ]
}