#' Run boot strapping to compute 95\% CI
#' @param idlist List of sample unit
#' @param seednum random seed number
#' @param multicore whether conducting multicore computation (default=FALSE)
#' @param mc.cores number of cores (default=2)
#' @param idvar ID variable name from the dataset
#' @param timevar time variable name from the dataset in days from the baseline
#' @param inYrs whether the time variable needs to be converted in years, default FALSE
#' @param predvarlist columns of the predictors
#' @param outcomevar outcome variable
#' @param varthresh threshold for LFPCA
#' @example
#' @references
#' TBA
#' @author
#' Seonjoo lee \email{sl2670@cumc.columbia.edu}
#' @export

lfpcr_boot<-function(datct=NULL,idlist=NULL,B=100,
                             idvar=NULL,timevar=NULL,inYrs=FALSE,predvarlist=NULL,outcomevar=NULL,covariates=NULL,
                             varthresh=0.85,
                             lambdalist = 2^(c(-8:8)), penalty.factor=NULL,nfold=10,
                     multicore=FALSE,mc.cores=2){
require(parallel)
  if(.Platform$OS.type=='windows'){
    print("multicore computing is under construction for windows environment." )
    z<-lapply(1:B,function(j){try(lfpcr_boot_oneunit(datct=datct,idlist=idlist,seednum=j,
                                                   idvar=idvar,timevar=timevar,inYrs=inYrs,predvarlist=predvarlist,outcomevar=outcomevar,covariates=covariates,
                                                   varthresh=varthresh,lambdalist = lambdalist, penalty.factor=penalty.factor,nfold=nfold))})
  }
  if(.Platform$OS.type=='unix'){
    print("parallel package is used for multicore computation." )
    z<-mclapply(1:B,function(j){try(lfpcr_boot_oneunit(datct=datct,idlist=idlist,seednum=j,
                                                   idvar=idvar,timevar=timevar,inYrs=inYrs,predvarlist=predvarlist,outcomevar=outcomevar,covariates=covariates,
                                                   varthresh=varthresh,lambdalist = lambdalist, penalty.factor=penalty.factor,nfold=nfold))},mc.cores=mc.cores)
  }
return(z)
#  return(list(beta0=lapply(z,function(x){x$lpcrre$beta0}), beta1=lapply(z,function(x){x$lpcrre$beta1})))

}
