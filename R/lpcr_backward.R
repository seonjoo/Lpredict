## backward stepwise selection
library(MASS)
lpcr_backword <- function(DX=DX.final, xi=lfpcre$xi, phix0=lfpcre$phix0, phix1=lfpcre$phix1, Nx=lfpcre$Nx){
  ##stepwise selection
  data.pcs <- data.frame(DX, t(xi))
  model.full <- glm(DX ~., family=binomial, data=data.pcs)
  model.step <- stepAIC(model.full, direction="backward", trace=FALSE)

  ##extract selected components
  colnames(phix0) <- paste('X', 1:Nx, sep='')
  colnames(phix1) <- paste('X', 1:Nx, sep='')
  rownames(xi) <- paste('X', 1:Nx, sep='')
  pcs.phix0 = phix0[,which(colnames(phix0) %in% attr(model.step$terms,"term.labels"))]
  pcs.phix1 = phix1[,which(colnames(phix1) %in% attr(model.step$terms,"term.labels"))]
  pcs.xi = xi[which(rownames(xi) %in% attr(model.step$terms,"term.labels")),]

  result = new.env()
  result$xi = pcs.xi
  result$phix0 = pcs.phix0
  result$phix1 = pcs.phix1
  result$Nx = dim(pcs.xi)[1]
  result=as.list(result)
  return(result)
}

## component selection with lasso
