#' Visualize the sweet spot analysis on clinical trial data. 
#' 
#' @name plot_sweetspot
#' @param sweetspot_result The result from a sweetspot analysis.
#' @param title The title for your plot.
#' @return A ggplot.
#' 
#' @import ggplot2
#' 
#' @export plot_sweetspot

library(ggplot2)

GREY   <- "#E5E5E5"
BLUE   <- "#62CBD7" 
GREEN1 <- "#9ED393"
GREEN2 <- "#55763C"
SPLINEGREY <- "#333333"
BASE_SIZE  <- 15

plot_sweetspot <- function(sweetspot_result, title=""){
  model    <- sweetspot_result$model
  matches  <- sweetspot_result$matches
  
  smoothed <- smooth.spline(matches[, "mean_score"], matches[, "apparent_benefit"])
    
  n     <- nrow(matches)
  
  zero  <- matches[1, "mean_score"]
  lrge  <- matches[n, "mean_score"]
  start <- matches[model$start.index, "mean_score"]
  end   <- matches[model$end.index,   "mean_score"]
  height.base <- abs(0.1*(max(model$mean.outside, model$mean.outside.debiased) - min(model$mean.inside, model$mean.inside.debiased)))
  base <- min(0, model$mean.outside-height.base, model$mean.outside.debiased-height.base)
  
  p <- ggplot() 
  p <- p + geom_ribbon(aes(x=c(zero, zero, start, start, end, end, lrge, lrge), 
                           ymin=rep(base, 8),
                           ymax=c(base, model$mean.outside, model$mean.outside, model$mean.inside, model$mean.inside, model$mean.outside, model$mean.outside, base),
                           fill="original"), alpha=.8)
  p <- p + geom_ribbon(aes(x=c(zero, zero, start, start, end, end, lrge, lrge), 
                           ymin=rep(base, 8),
                           ymax=c(base, model$mean.outside.debiased, model$mean.outside.debiased, model$mean.inside.debiased, model$mean.inside.debiased, model$mean.outside.debiased, model$mean.outside.debiased, base),
                          fill="debiased"), alpha=.4)
  p <- p + geom_line(aes(x=smoothed$x, y=smoothed$y, color="smoothed\ntreatment effect"), size=.6)
  p <- p + xlim(zero, lrge)
  p <- p + scale_color_manual(breaks=c("debiased", "original", "start", "end", "smoothed\ntreatment effect"), values=c(BLUE, GREY, GREEN1, GREEN2, SPLINEGREY))
  p <- p + scale_fill_manual( breaks=c("debiased", "original", "start", "end"), values=c(BLUE, GREY, GREEN1, GREEN2))
  p <- p + geom_rug(aes(x=matches[sample(model$start.indices, min(n, 1000)), "mean_score"], color="start"), size=1, length = unit(0.05, "npc"), alpha=.1)
  p <- p + geom_rug(aes(x=matches[sample(model$end.indices, min(n, 1000)), "mean_score"],   color="end"),   size=1, length = unit(0.05, "npc"), alpha=.1)
  p <- p + labs(x = "Predilection score", y = "Treatment effect estimate", x="", y="", color="", fill="", 
                subtitle=paste0("p-value: ", model$p.value), title=title) 
  p <- p + theme_bw(base_size=BASE_SIZE) 
  p <- p + theme(
              axis.line = element_line(colour = "#d5d5d5"),
              panel.grid.major.y = element_line(linetype = "dotted", color="#d5d5d5"),
              panel.grid.major.x = element_line(linetype = "dotted", color="#d5d5d5"),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank()) 
  p
}
