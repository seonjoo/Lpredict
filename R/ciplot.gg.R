#' Draw 95\% Confidence Interval
#'
#' @param d Bootstrapped coefficients matrix.
#' @param exp Whether take exponential for the estimated CIs.
#' @param yrange set up range for the 95\% CIs.
#' @examples
#' ciplot.gg(d)
#' @author Liwen Wu, Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references
#' TBA
#' @keywords confidence interval
#' @export

ciplot.gg<-function (d, exp = FALSE, yrange = c(-5, 5),
                     group=c('Right','Left'),
                     group.label=c("Right Hemisphere", "Left Hemisphere"))
{
  require(matrixStats)
  require(ggplot2)

  if (is.null(d)) {
    print("d is missing.")
    break
  }
  namelist = colnames(d)
  if (is.null(namelist)){namelist=1:ncol(d)}

  citable <- data.frame(x = namelist, y = colMedians(d),
                        ylo = colQuantiles(d, probs = 0.025),
                        yhi = colQuantiles(d, probs = 0.975))
  if (exp == TRUE) {
    citable <- transform(citable, y = exp(y), ylo = exp(ylo),
                         yhi = exp(yhi))
  }
  citable <- transform(citable, right = grepl(group[1], citable[,
                                                               1]))
  citable <- transform(citable, right = ifelse(right == TRUE,
                                               group.label[1],group.label[2]), x = gsub(group[1],
                                                                                                "", x))
  citable <- transform(citable, x = gsub(group[2], "", x))
  if (exp == TRUE) {
    p <- ggplot(citable, aes(x = x, y = y, ymin = ylo, ymax = yhi)) +
      geom_pointrange() + facet_grid(. ~ right) + coord_flip() +
      geom_hline(aes(yintercept = 1), lty = 2, col = 2) +
      xlab("") + ylab("") + theme(axis.text = element_text(size = 14),
                                  strip.text = element_text(size = 18, face = "bold")) +
      ylim(yrange)
  }
  if (exp == FALSE) {
    p <- ggplot(citable, aes(x = x, y = y, ymin = ylo, ymax = yhi)) +
      geom_pointrange() + facet_grid(. ~ right) + coord_flip() +
      geom_hline(aes(yintercept = 0), lty = 2, col = 2) +
      xlab("") + ylab("") + theme(axis.text = element_text(size = 14),
                                  strip.text = element_text(size = 18, face = "bold")) +
      ylim(yrange)
  }
  return(p)
}
