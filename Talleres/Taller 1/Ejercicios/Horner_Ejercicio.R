#     *** Analisis Numerico 2020-1 *** 
#--Trabajo realizado por: 
# Monica A. Alvarez C.
# Santiago Palacios L.
#   Ejercicio 5). Evaluacion de un polinomio con el teorema de Horner  

Horner <- function(coeficientes, x0){
  valor = coeficientes[1] # Ajuste de las variables 
  multi =0
  sumas=0
  for(i in coeficientes[2:length(coeficientes)]){ # Ciclo que se repite por la cantidad de 
    valor <- x0*valor + i                         # terminos del polinomio propuesto.
    multi = multi + 1
    sumas = sumas + 1
    
  }
  return(cat("El valor del polinomio es: ", round(valor,4), "\nEl total de operaciones es de: ",
  multi+sumas,"; Multiplicaciones:",multi,"; Sumas:  ",sumas))
  
}
x0 = -2
coeficiente <- c(2,0,-3,3,4)
Horner(coeficiente,x0)