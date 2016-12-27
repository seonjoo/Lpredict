#' Draw 95\% Confidence Interval
#'
#' @param d Bootstrapped coefficients
#' @param lfpcrobj lfpcrobj
#' @examples

#' @author Liwen Wu, Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references
#' TBA
#' @keywords
#' hdlfpca glmnet
#' @export

ciplot.gg <- function(lfpcrobj, d){
  require(matrixStats)
  require(ggplot2)
  if (is.null(lfpcrobj)){print('lfpcr object is missing.');break}
  if (is.null(d)){print('bootstrapped object is missing.');break}

  namelist = rownames(lfpcrobj$phix0)
  if (is.null(namelist)) namelist<-paste('X',1:nrow(lfpcrobj$phix0))

  citable <- data.frame(x=rownames(lfpcrobj$phix0),y=colMedians(d),ylo=colQuantiles(d,probs=0.025),yhi=colQuantiles(d,probs=0.975))
  citable <- transform(citable,y=exp(y),ylo=exp(ylo),yhi=exp(yhi))
  citable <- transform(citable,right=grepl('Right',citable[,1]))
  citable <- transform(citable, right =ifelse(right == TRUE,'Right Hemisphere','Left Hemisphere'), x=gsub('Right','',x))
  citable <- transform(citable, x=gsub('Left','',x))
  p <- ggplot(citable, aes(x=x, y=y, ymin=ylo, ymax=yhi))+geom_pointrange() + facet_grid(.~right) +
    coord_flip() + geom_hline(aes(yintercept=1), lty=2,col=2) + xlab('') + ylab('') +
    theme(axis.text=element_text(size=14), strip.text=element_text(size=18,face="bold"))
  return(p)
}
