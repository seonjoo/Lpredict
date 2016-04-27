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
#' @references TBA
#' @keywords hdlfpca glmnet
#' @export

lpcr_lasso <- function(Y,T=NULL,J=NULL,I=NULL,visit=NULL, lfpca=NULL, cov=NULL,
                       nfold=10,lambdalist = 2^(c(-10:10)/2), penalty.factor=NULL, M=NULL,
                       seednum=1234){
  library(glmnet)
  ## If lfpca is has not been conducted, it will run.
  if (is.null(lfpca)){
    if (is.null(T) |  is.null(J) | is.null(I) | is.null(visit)){
      print('LFPCA has not been run, and the required parameters for LFPCA were not specified.')
      break
    }else{
        system.time(lfpc<-hd_lfpca(Y,T,J,I,visit,varthresh = 0.85, timeadjust=FALSE, figure=TRUE))
    }
  }

  if(is.null(ncol(cov))){
    if(is.null(cov)) ncov=0
    else{
      ncov=1
      }
    }else{
    ncov=ncol(cov)
  }

  rownames(lfpca$xi)<-paste('LPC',1:lfpca$Nx,sep='')
  X = cbind(cov,t(lfpca$xi))
#  print(head(X))
  Nx = lfpca$Nx

  if(is.null(penalty.factor)){
    if (is.null(cov)){
      penalty.factor=rep(1,ncol(X))
    }else{
      penalty.factor=c(rep(0,ncov),rep(1,ncol(X)))
    }
  }

  ## lasso
  set.seed(seednum)
  fitcv=cv.glmnet(X,Y, family='binomial', alpha = 1, lambda = lambdalist, nfold=nfold, penalty.factor=penalty.factor)
  fit = glmnet(X,Y, family='binomial', alpha = 1, lambda =fitcv$lambda.min, penalty.factor=penalty.factor)
  coefs <- coef(fit)

  colnames(lfpca$phix0) <- paste('LPC', 1:lfpca$Nx, sep='')
  colnames(lfpca$phix1) <- paste('LPC', 1:lfpca$Nx, sep='')

  nonzeroindx =  which(coefs[-1]!=0)

  nonzeroindx2 =  (coefs[-c(0:ncov +1)]!=0)
## extract selected components
#  if (is.null(M)) {
#    logic <- (coef!=0)[-1]
#  } else{
#    logic <- (coef!=0)[-c(0:dim(M)[2]+1)]
#  }

  refit = glmnet(X[,nonzeroindx],Y, family='binomial', alpha = 1, lambda =0)


  hattheta = as.matrix(coef(refit))

#  print(dim(lfpca$phix0[,nonzeroindx2]))
#  print(hattheta)
#  print(ncov)
  beta0 = lfpca$phix0[,nonzeroindx2] %*% hattheta[-c(0:ncov +1)]
  beta1 = lfpca$phix0[,nonzeroindx2] %*% hattheta[-c(0:ncov +1)]
  return(list(nonzeroindx = nonzeroindx, nonzeroindx2 = nonzeroindx2,refit=refit, beta0=beta0,beta1=beta1,hattheta=hattheta))

}
