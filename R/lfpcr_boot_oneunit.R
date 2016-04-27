#' Run one unit boot strapping
#' @param idlist List of sample unit
#' @param seednum random seed number
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
lfpcr_boot_oneunit<-function(datct=NULL,idlist=NULL,seednum=1000,
                             idvar=NULL,timevar=NULL,inYrs=FALSE,predvarlist=NULL,outcomevar=NULL,covariates=NULL,
                             varthresh=0.85,
                             lambdalist = 2^(c(-8:8)), penalty.factor=NULL,nfold=10){
  if(is.null(idvar)){print("ID variable has to be specified.");break}
  if(is.null(datct)){print("dataset (datct) has to be specified.");break}
  if(is.null(idlist)){print("Unique ID (idlist) list is not specified. We are computing from the id variable.");idlist <- unique(datct[,idvar])}
  if(is.null(timevar)){print("time variable (timevar) has to be specified.");break}
  if(is.null(predvarlist)){print("predictor list (predvarlist) have to be specified.");break}
  if(is.null(outcomevar)){print("An outcome variable (outcomevar) has to be specified.");break}

  set.seed(seednum)
  samplelist =sort(sample(idlist, size=length(idlist),replace=TRUE))
  indx=unlist(lapply(samplelist, function(x){which(datct[,idvar]==x) } ))
  visit=unlist(lapply(samplelist, function(x){length(which(datct[,idvar]==x))}))

  ##lfpca
  Y <- as.matrix(datct[indx,predvarlist])
  if (inYrs==TRUE){
    T=as.vector(datct[indx,timevar]/365)
  }else{
    T=as.vector(datct[indx,timevar])
  }
  J=length(indx)
  I=length(unique(datct[,idvar]))

  lfpcre<-hd_lfpca(t(Y),T,J,I,visit,varthresh=varthresh)

  ##principle component selection
  tempdata = datct[indx,]

  outcome = tempdata[tempdata[,timevar]==0, outcomevar]
  if(is.null(covariates)){
    lpcrre=lpcr_lasso(Y=outcome, lfpca=lfpcre, nfold=nfold,lambdalist = lambdalist, penalty.factor=penalty.factor,seednum=seednum)
  }else{
    lpcrre=lpcr_lasso(Y=outcome,cov=as.matrix(tempdata[tempdata[,timevar]==0, covariates]), lfpca=lfpcre, nfold=nfold,lambdalist = lambdalist, penalty.factor=penalty.factor,seednum=seednum)
  }

  return(list(lpcrre=lpcrre, lfpcre=lfpcre))
}
