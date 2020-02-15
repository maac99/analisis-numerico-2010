#     *** Analisis Numerico 2020-1 *** 
#--Trabajo realizado por: 
# Monica A. Alvarez C.
# Santiago Palacios L.
#   Algoritmo para sumar los n^2.

binario_a_entero <- function (parte_entera) {
  
  numero_potencia = 0
  binario = 0
  total = 0
  
  while(parte_entera > 0) {
    
    binario = parte_entera %% 10
    total = total + (binario*(2 ^ numero_potencia))
    numero_potencia = numero_potencia + 1
    parte_entera = parte_entera %/% 10
  }
  
  return (total)
  
}

binario_a_entero(101)


AlCuadrado <- function (x)
{
  i=0
  resultados<-c()
  while(i<x){
    resultados = c(resultados,(i^2))
    i=i+1
  }
  return (resultados)
}

eficiencia <- function(x){
  n = x
  tot = 0
  while(n>0){
    d = n %% 2
    n = n%/%2
    tot = tot + 2
    cat(d)
  }
  cat("\n")
  #cat(" ",x," ",tot," \n")
}

eficiencia(2)
eficiencia(10)

resultados<-c()
resultados=AlCuadrado(7)
eAnterior<-c()
eActual<-c()

eficiencia(10)
binario_a_entero(eficiencia(10))

binario_a_entero(eficiencia(resultados[2]))


eActual<-c(eActual,binario_a_entero(eficiencia(resultados[2])))
eAnterior<-c(eAnterior,binario_a_entero(eficiencia(resultados[1])))

eActual<-c(eActual,binario_a_entero(eficiencia(resultados[3])))
eAnterior<-c(eAnterior,binario_a_entero(eficiencia(resultados[2])))

eActual<-c(eActual,binario_a_entero(eficiencia(resultados[4])))
eAnterior<-c(eAnterior,binario_a_entero(eficiencia(resultados[3])))

eActual<-c(eActual,binario_a_entero(eficiencia(resultados[5])))
eAnterior<-c(eAnterior,binario_a_entero(eficiencia(resultados[4])))


eActual<-c(eActual,binario_a_entero(eficiencia(resultados[6])))
eAnterior<-c(eAnterior,binario_a_entero(eficiencia(resultados[5])))

eActual<-c(eActual,binario_a_entero(eficiencia(resultados[7])))
eAnterior<-c(eAnterior,binario_a_entero(eficiencia(resultados[6])))



plot(eActual, # Grafica del ejercicio.
     eAnterior,
     main = "Error Actual vs Anterior",
     xlab = "Error Anterior",
     ylab = "Error Actual")
points (eActual, eAnterior, col = "red")
lines(eActual, eAnterior, col = "blue")
