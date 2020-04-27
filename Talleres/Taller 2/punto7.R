#  *** Análisis numérico 2020 - 1 ***
# Trabajo realizado por:
#  Mónica A. Álvarez C.
#  Santiago Palacios L.
#  Paula C. Piñeros P.
# Implementacio vectorial de polinomio de Lagrange.
f <- function( x ) { exp(x) }
lagrange = function(x,y,a){ 
  n = length(x) 
  if(a < min(x) || max(x) < a) stop("No está interpolando")
  X = matrix(rep(x, times=n), n, n, byrow=T)
  mN = a - X; diag(mN) = 1
  mD = X - t(X); diag(mD) = 1
  Lnk = apply(mN, 1, prod)/apply(mD, 2, prod)
  sum(y*Lnk)
}


x<-c(0,0.333,0.666,0.999)
y<-f(x)
lagrange(x,y,0.5)

x<-c(0,0.25,0.5,0.75,1)
y<-f(x)
lagrange(x,y,0.5)
