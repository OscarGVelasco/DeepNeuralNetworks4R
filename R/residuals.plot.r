#' Plot of the residuals Results
#'
#' This function plots the residuals from the results once the training has finish.
#'
#' @import ggplot2
#' @param observed Vector. The observed values from the data training.
#' @param predicted Vector. The predicted values from the model's result.
#' @param title Character. The title to be shown in the plot.
#' @export
residuals_plot <- function(observed = NULL, predicted = NULL, title = "Residuals plot"){
  library(ggplot2)
  if(!is.vector(observed)){
    if(ncol(observed)==1)observed <- c(observed)
    else stop("observed values need to be a vector or a df with a single column")
  }
  if(!is.vector(predicted)){
    if(ncol(predicted)==1)predicted <- c(predicted)
    else stop("predicted values need to be a vector or a df with a single column")
  }
  if(length(observed)!=length(predicted))stop("The lengths of the vectors does not match.")
  residuals <- observed - predicted
  # r.test <- runs.test(factor(residuals>0), alternative = c("two.sided"))
  labels <- paste(
    paste("runs test for randomness of residuals:"),
    # paste("- p.value =", round(r.test$p.value,digits = 2)),
    # paste("- alternative hypothesis ",r.test$alternative),
    sep="\n")

  residuals.studentized <- residuals/sd(residuals)
  df <- data.frame(residuals = residuals.studentized,
                   predicted = predicted,
                   deviation = abs(residuals.studentized),
                   sample.name = names(residuals.studentized),
                   stringsAsFactors = F)

  g <- ggplot(df, aes(x = predicted, y = residuals,color=deviation)) +
    geom_point(size = 2, shape = 16) +
    geom_hline(yintercept = 0) +
    # We print the outliers that are > 1.5 * IQR
    geom_text(aes(label=ifelse(sample.name%in%names(boxplot.stats(residuals.studentized)$out),as.character(sample.name),''),colour=0),hjust=0,vjust=0) +
    labs(title = title,
         x = "predicted values", y = "studentized residuals") +
    theme(legend.justification = c(0, 1), legend.position = c(0, 1),
          axis.text = element_text(size = 22), axis.title=element_text(size=20),
          plot.title = element_text(size=24,face="bold",hjust = 0.5))
  return(g)
}
