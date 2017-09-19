#' Conduct longitudinal principal component regression with model selection
#'
#' @param Y Outcome in PCR
#' @param Xmat Longitudinal Predictor Matrix
#' @param cov Covariate
#' @param lfpca (default=NULL) LFPC results. If it is not specified, the following parameters have to be specified to run LPCA
#' @param T Time of theimage collection
#' @param J Total number of observations
#' @param I Total number of subjects
#' @param visit Vector of number of visits per subjects
#' @param verbose (default=FALSE)
#' @param method (default=AUC) backward selection criteria. c('AUC','bic','aic')
#' @param Nx Dimension of the subject-specific components
#' @param Nw Dimension of the subject-visit specific components
#' @param varthresh (default=0.99) Threshold for variance explained for both subject-specific and subject-visit specific compoents for dimension selection
#' @param projectthresh Threshold for variance explain in the first step of SVD
#' @param timeadjust (default=TRUE) Scale time per subject
#' @return xi
#' @return phix0
#' @return phix1
#' @return zeta
#' @return phiw
#' @examples
#' I=200
#' visit=rpois(I,1)+3
#' time = unlist(lapply(visit, function(x) scale(c(0,cumsum(rpois(x-1,1)+1)))))
#' J = sum(visit)
#' V=1000
#' phix0 = matrix(0,V,3);phix0[1:50,1]<-.1;phix0[1:50 + 50,2]<-.1;phix0[1:50 + 100,3]<-.1
#' phix1 = matrix(0,V,3);phix1[1:50+150,1]<-.1;phix1[1:50 + 200,2]<-.1;phix1[1:50 + 250,3]<-.1
#' phiw = matrix(0,V,3);phiw[1:50+300,1]<-.1;phiw[1:50 + 350,2]<-.1;phiw[1:50 + 400,3]<-.1
#' xi = t(matrix(rnorm(I*3),ncol=I)*c(8,4,2))*3
#' zeta = t(matrix(rnorm(J*3),ncol=J)*c(8,4,2))*2
#' Xmat = phix0%*% t(xi[rep(1:I, visit),]) + phix1%*% t(time * xi[rep(1:I, visit),]) + phiw %*% t(zeta) + matrix(rnorm(V*J,0,.1),V,J)
#' beta=c(1,-1,1)/10
#' p=exp( - xi %*% beta)/(1 + exp( - xi %*% beta))
#' Y = unlist(lapply(p,function(x)rbinom(1,1,prob=x)))
#' re<-hd_lfpca(Ydat,Xmat,T,J,I,visit, varthresh=0.85, timeadjust=FALSE)
#' lpcr.backward(Y,lfpca=re)
#' @author Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references TBA
#' @keywords hdlfpca glmnet
#' @impot pROC
#' @export

lpcr.backward <- function(Y,Xmat=NULL,lfpca=NULL, cov=NULL,
                       T=NULL,J=NULL,I=NULL,visit=NULL,
                       varthresh=0.85, timeadjust=FALSE,
                       nfold=10,lambdalist = 2^(c(-10:10)/2), penalty.factor=NULL, M=NULL,
                       method='AUC',
                       seednum=1234){
  library(glmnet)
  ## If lfpca is has not been conducted, it will run.
  if (is.null(lfpca)){
    if (is.null(Xmat) | is.null(T) |  is.null(J) | is.null(I) | is.null(visit)){
      print('LFPCA has not been run, and the required parameters for LFPCA were not specified.')
      break
    }else{
      system.time(lfpca<-hd_lfpca(Y=Xmat,T=T,J=J,I=I,visit=visit,varthresh = varthresh, timeadjust=timeadjust))
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
  rownames(lfpca$xi)<-paste('xi',1:lfpca$Nx,sep='')

  ## lasso
  library(pROC)
  aucs=runaucs(lfpcobj=lfpca, Y =Y, M=M,method='fast')
  aucvarlists = which(paste('xi',1:lfpca$Nx,sep='') %in% aucs$droppedPC[1:which.max(aucs$AUC)] ==FALSE)

  tmpdata0=data.frame(Y=Y,t(lfpca$xi)[,aucvarlists])
  fitfinal=glm(Y~. ,data=tmpdata0, family=binomial(link = "logit"))
  beta0 = lfpca$phix0[,aucvarlists] %*% as.matrix(coef(fitfinal)[-1])
  beta1 = lfpca$phix1[,aucvarlists] %*% as.matrix(coef(fitfinal)[-1])

  return(list(nonzeroindx = aucvarlists, beta0=beta0,beta1=beta1))
}




prediction<-function(lfpcre=NULL,X=NULL,M=NULL,DX.final=NULL,subset=1:20){
  if (is.null(lfpcre) & is.null(X) ){stop('Either LFPC object or predictor X or subset have to be specified.')}
  if (is.null(lfpcre)==FALSE & is.null(X) & length(subset)>0){X = t(lfpcre$xi[subset,])}
  if (is.null(lfpcre)==FALSE & is.null(X) & length(subset)==1){X = (lfpcre$xi[subset,])}
  if (is.null(lfpcre)==FALSE & is.null(X) & length(subset)==0){X = t(lfpcre$xi)}

  N=length(DX.final)
  pred=rep(NA, N)
  ### leave-one-out cross-validation

  #  print(data.frame(Y=DX.final[1],M=t(as.matrix(M[1,])), x=t(as.matrix(X[1,]))))
  for (j in 1:N){
    if (is.null(M)==FALSE){
      if(length(subset)>1){tmpdata=data.frame(Y=DX.final[-j],M=as.matrix(M)[-j,],x=X[-j,])}
      if(length(subset)==1){tmpdata=data.frame(Y=DX.final[-j],M=as.matrix(M)[-j,],x=X[-j])}
      if(length(subset)==0){tmpdata=data.frame(Y=DX.final[-j], M=as.matrix(M)[-j,])}
    }
    if(is.null(M)==TRUE){
      if(length(subset)>1){tmpdata=data.frame(Y=DX.final[-j],x=X[-j,])}
      if(length(subset)==1){tmpdata=data.frame(Y=DX.final[-j],x=X[-j])}
      if(length(subset)==0){tmpdata=data.frame(Y=DX.final[-j])}
    }

    fit=glm(Y~., family='binomial',tmpdata)
    #  print(fit)
    if (is.null(M)==FALSE){
      if(length(subset)>1){newdata=data.frame(Y=DX.final[j],M=(as.matrix(M[j,])), x=t(as.matrix(X[j,])))}
      if(length(subset)==1){newdata=data.frame(Y=DX.final[j],M=(as.matrix(M[j,])), x=X[j])}
      if(length(subset)==0){newdata=data.frame(Y=DX.final[j],M=(as.matrix(M[j,])))}
    }
    if(is.null(M)==TRUE){
      if(length(subset)>1){newdata=data.frame(Y=DX.final[j], x=t(as.matrix(X[j,])))}
      if(length(subset)==1){newdata=data.frame(Y=DX.final[j], x=X[j])}
      if(length(subset)==0){newdata=data.frame(Y=DX.final[j])}
    }
    # print(length(subset))
    #print(newdata)
    pred[j]<-predict( fit, newdata=newdata, type="response" )
  }
  return(pred)
}

run_lpredict_backward2<-function (lfpcobs, seednum=1234,idlist, Y, cov=NULL,method='fast',selmethod='BIC', start=NULL)
{
  if(is.null(cov)==FALSE){numcov=ncol(cov);if(is.null(numcov)==TRUE){numcov=1}}
  if (is.null(cov)==TRUE){numcov=0}
  set.seed(seednum)
  samplelist = sort(sample(1:length(idlist), size = length(idlist), replace = TRUE))

  tmpY = Y[samplelist]
  tmpcov=cov
  if(is.null(cov)==FALSE){tmpcov=(cov[samplelist,])}

  tmpaucs=runaucs(lfpcobj=lfpcobs, Y =tmpY, M=tmpcov,method=method,sampleorder=samplelist)

  xi = t(lfpcobs$xi)[samplelist,]


  if (selmethod=='AUC'){ maxindx = which.max(tmpaucs$AUC) }
  if (selmethod=='BIC'){ maxindx = which.min(tmpaucs$BIC) }
  namelist = paste('xi',1:lfpcobs$Nx,sep='')
  varlists = which(namelist %in% tmpaucs$droppedPC[1:maxindx] ==FALSE)

  tmpdata0=data.frame(Y=tmpY,M=tmpcov,xi[,varlists])
  run<-tryCatch({fitfinal<-glm(Y~. ,data=tmpdata0, family='binomial',maxit=150)}, warning=function(x){NA} )
  if(is.na(run)==FALSE){
    beta0 = lfpcobs$phix0[,varlists] %*% as.matrix(coef(fitfinal)[-c(1:(1+numcov))])
    beta1 = lfpcobs$phix1[,varlists] %*% as.matrix(coef(fitfinal)[-c(1:(1+numcov))])
    return(list(fit=fitfinal,beta0=beta0,beta1=beta1, varlists=varlists,aucs=tmpaucs))
  }else{return(NULL)}
}


runaucs<-function(lfpcobj=lfpcobj, Y =NULL, M=NULL,method='bootstrap',sampleorder=1:length(Y),start=NULL){
  rownames(lfpcobj$xi)<-paste('xi',1:lfpcobj$Nx,sep='')
  xi=t(lfpcobj$xi)[sampleorder,]
  namelist = rownames(lfpcobj$xi)
  N = length(Y)
  numcov = 0;if(is.null(M)==FALSE){numcov=ncol(M)}
  droplist=c();aucs=c();bic<-c()
  nulval=100;minval=0;

  if (is.null(M)==TRUE){niter = lfpcobj$Nx-2}
  if (is.null(M)==FALSE){niter = lfpcobj$Nx-1}

  for (k in 1:niter){
    if (is.null(M)==FALSE){tmpdata0=data.frame(Y=Y,M=as.matrix(M),xi[,namelist %in% droplist ==FALSE])}
    if (is.null(M)==TRUE){tmpdata0=data.frame(Y=Y,xi[,namelist %in% droplist ==FALSE])}
    fit0=glm(Y~. ,data=tmpdata0, family='binomial',start=start, maxit=150)
    a=drop1(fit0,k=log(N))
    #    print(a)

    dropvar=which.min(a$AIC[-c(1:(numcov+1))])
    minval=min(a$AIC[-c(1:(numcov+1))])
    nulval=a$AIC[1]
    tmpname =colnames(xi[,namelist %in% droplist ==FALSE])
    print(c(k,tmpname[dropvar]))

    if (is.null(tmpname)==TRUE){droplist=c(droplist,namelist[includevarlist])}
    if (is.null(tmpname)==FALSE){droplist=c(droplist,tmpname[dropvar])}

    includevarlist = which(namelist %in% droplist ==FALSE)
    print(includevarlist)
    if( is.null(includevarlist)==FALSE){pred=prediction(lfpcre=lfpcobj,M=M,DX.final=Y,subset=includevarlist)}
    if( is.null(includevarlist)==TRUE){pred=prediction(X=M,DX.final=Y)}

    print(c(length(DX.final), length(pred)))
    rocobj <- try(roc(DX.final, pred))
    set.seed(123456)
    if (length(rocobj)>1){
      if (method=='bootstrap'){aucs<-cbind(aucs,ci.auc(rocobj, method='bootstrap',boot.n=5000))}
      if (method=='fast'){aucs<-cbind(aucs,ci.auc(rocobj))}
    }
    bic<-c(bic, minval)
  }

  aucdat=data.frame(t(aucs))
  names(aucdat)=c('LCI','AUC','UCI')
  aucdat$droppedPC=droplist
  aucdat$x=1:niter
  aucdat$BIC=bic

  return(aucdat)

}
