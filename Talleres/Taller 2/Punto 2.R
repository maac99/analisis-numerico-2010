#     * Analisis Numerico 2020-1 *
#--Trabajo realizado por:
# Monica A. Alvarez C.
# Santiago Palacios L.
# Paula C. Piñeros P.

x = c(0,  1,  2)
y = c(10, 15, 5)
a = 2

MatrizDD = function(x, y){
  n = length(x)
  A = matrix(rep(NA, times = n^2), nrow = n, ncol = n)  #A es la matriz de diferenciass divididas
  A[,1] = y
  for(i in 2:n){
    A[i:n, i] = ((A[i:n, i - 1] - A[(i - 1):(n-1), i - 1]) / (x[i:n] - x[1:(n - i + 1)]))
  }
  print(A)
  return(A)
}

A = MatrizDD(x,y)

w = ncol(A)
m = length(x)
xx = c()
for(i in 1:m){
  xx = c(xx, "(x - ", x[i], ")")
}

cat("El polinomio es: ")
for(i in 2:w){
  if(i == 2){
    cat(A[1,1], " + ")
  }
  cat(A[i,i])
  j = 1
  while(j <=(i-1)*3){
    cat(xx[j],xx[j+1],xx[j+2])
    j = j + 3
  }
  if(i != w){
    cat(" + ")
  }
}
