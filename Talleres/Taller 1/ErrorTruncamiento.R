#     * Analisis Numerico 2020-1 * 
#--Trabajo realizado por: 
# Monica A. Alvarez C.
# Santiago Palacios L.
#   Ejercicio 1). Error de truncamiento.

error <- function(num){
  m = 4 #numero de cifras que se pueden almacenar
  n = 0 
  numtemp = num
  while(numtemp>1){
    numtemp = numtemp/10
    n = n+1 #n representa la cantidad de enteros del numero normalizado en potencias de 10
  }
  
  numtemptrunc <- trunc(numtemp*10^m)/10^m  #se trunca el numero de acuerdo al numero de cifras
  errort <- (numtemp-numtemptrunc)*10^(n-1)  #se calcula el error con el valor original y el valor truncado
  cat("El error de truncamiento es de: ", errort)
}