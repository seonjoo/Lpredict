#' Conduct longitudinal principal component regression analysis
#'
#' @param Y A number
#' @param T A number
#' @param J A number
#' @param I A number
#' @param visit A number
#' @param cov Covariate
#' @param varthresh A number
#' @param projectthresh A number
#' @param timeadjust A number
#' @return xi
#' @return phix0
#' @return phix1
#' @return zeta
#' @return phiw
#' @examples
#' re<-hd_lfpca(Ydat,T,J,I,visit, varthresh=0.85, timeadjust=FALSE)
#' lpcr_lasso(Y,lfpca=re)
#' @author Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references
#' @keywords hdlfpca glmnet
#' @export

lpcr_lasso <- function(Y,T=NULL,J=NULL,I=NULL,visit=NULL, lfpca=NULL, cov=NULL, Nx,nfold=10,lambdalist = 2^(c(-10:10)/2), penalty.factor=NULL, M=NULL){
  library(glmnet)
  ## If lfpca is has not been conducted, it will run.
  if (is.null(lfpca)){
    if (is.null(T) |  is.null(J) | is.null(I) | is.null(visit)){
      print('LFPCA has not been run, and reuiqred parameters were not specified.')
      break
    }else{
        system.time(lfpc<-hd_lfpca(Y,T,J,I,visit,varthresh = 0.85, timeadjust=FALSE, figure=TRUE))
    }
  }

  X = lfpc$xi

  if(is.null(penalty.factor)){
    if (is.null(cov)){
      penalty.factor=rep(1,ncol(X))
    }else{
      if(is.null(ncol(cov))){
        ncov=1
      }else{
        ncov=ncol(cov)
      }
      penalty.factor=c(rep(0,ncov),rep(1,ncol(X)))
    }
  }

  ## lasso
  fitcv=cv.glmnet(t(X),Y, family='binomial', alpha = 1, lambda = lambdalist, nfold=nfold, penalty.factor=penalty.factor)
  fit = glmnet(t(X),Y, family='binomial', alpha = 1, lambda =fitcv$lambda.min, penalty.factor=penalty.factor)
  coef <- coef(fit)

  colnames(lfpca$phix0) <- paste('X', 1:lfpca$Nx, sep='')
  colnames(lfpca$phix1) <- paste('X', 1:lfpca$Nx, sep='')
  if (is.null(M)) {
    rownames(X) <- paste('X', 1:lfpca$Nx, sep='')
  } else {rownames(X)[-c(1:dim(M)[2])] <- paste('X', 1:lfpca$Nx, sep='')}

  ## extract selected components
  if (is.null(M)) {
    logic <- (coef!=0)[-1]
  } else{
    logic <- (coef!=0)[-c(0:dim(M)[2]+1)]
  }
  phix0 = lfpca$phix0[,which(logic==1)]
  phix1 = lfpca$phix1[,which(logic==1)]
  Nx = sum(logic)
  if (!is.null(M)) {X = X[-c(1:dim(M)[2]),]}
  xi = X[which(logic==1),]


  return(list(  xi = xi,
                phix0 = phix0,
                phix1 = phix1,
                Nx = dim(xi)[1] ))
}
