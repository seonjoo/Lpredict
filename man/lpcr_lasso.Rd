% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lpcr_bic.R, R/lpcr_lasso.R
\name{lpcr_lasso}
\alias{lpcr_lasso}
\alias{lpcr_lasso}
\title{Conduct longitudinal principal component regression analysis}
\usage{
lpcr_lasso(Y, Xmat = NULL, lfpca = NULL, cov = NULL,
  type.measure = "auc", T = NULL, J = NULL, I = NULL, visit = NULL,
  varthresh = 0.85, timeadjust = FALSE, nfold = 10,
  lambdalist = 2^(c(-10:10)/2), penalty.factor = NULL, M = NULL,
  seednum = 1234)

lpcr_lasso(Y, Xmat = NULL, lfpca = NULL, cov = NULL,
  type.measure = "auc", T = NULL, J = NULL, I = NULL, visit = NULL,
  varthresh = 0.85, timeadjust = FALSE, nfold = 10,
  lambdalist = 2^(c(-10:10)/2), penalty.factor = NULL, M = NULL,
  seednum = 1234)
}
\arguments{
\item{Y}{Outcome in PCR}

\item{Xmat}{Longitudinal Predictor Matrix}

\item{lfpca}{(default=NULL) LFPC results. If it is not specified, the following parameters have to be specified to run LPCA}

\item{cov}{Covariate}

\item{type.measure}{(default= "auc")}

\item{T}{Time of theimage collection}

\item{J}{Total number of observations}

\item{I}{Total number of subjects}

\item{visit}{Vector of number of visits per subjects}

\item{varthresh}{(default=0.99) Threshold for variance explained for both subject-specific and subject-visit specific compoents for dimension selection}

\item{timeadjust}{(default=TRUE) Scale time per subject}

\item{verbose}{(default=FALSE)}

\item{Nx}{Dimension of the subject-specific components}

\item{Nw}{Dimension of the subject-visit specific components}

\item{projectthresh}{Threshold for variance explain in the first step of SVD}

\item{Y}{Outcome in PCR}

\item{Xmat}{Longitudinal Predictor Matrix}

\item{cov}{Covariate}

\item{lfpca}{(default=NULL) LFPC results. If it is not specified, the following parameters have to be specified to run LPCA}

\item{T}{Time of theimage collection}

\item{J}{Total number of observations}

\item{I}{Total number of subjects}

\item{visit}{Vector of number of visits per subjects}

\item{verbose}{(default=FALSE)}

\item{Nx}{Dimension of the subject-specific components}

\item{Nw}{Dimension of the subject-visit specific components}

\item{varthresh}{(default=0.99) Threshold for variance explained for both subject-specific and subject-visit specific compoents for dimension selection}

\item{projectthresh}{Threshold for variance explain in the first step of SVD}

\item{timeadjust}{(default=TRUE) Scale time per subject}
}
\value{
xi

phix0

phix1

zeta

phiw

nonzeroindx Non-zero index including covariates

nonzeroindx2 Non-zero index excluding covariates

refit refitted regression object

beta0 Estimated parameter for intercept

beta1 Estimated parameter for slope

auc : if type.measure = "auc", cross-validated
}
\description{
Conduct longitudinal principal component regression analysis

Conduct longitudinal principal component regression analysis
}
\examples{
re<-hd_lfpca(Ydat,Xmat,T,J,I,visit, varthresh=0.85, timeadjust=FALSE)
lpcr_lasso(Y,lfpca=re)
re<-hd_lfpca(Ydat,Xmat,T,J,I,visit, varthresh=0.85, timeadjust=FALSE)
lpcr_lasso(Y,lfpca=re)
}
\references{
TBA

TBA
}
\author{
Seonjoo Lee, \email{sl3670@cumc.columbia.edu}

Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
}
\keyword{glmnet}
\keyword{hdlfpca}
