#     *** Analisis Numerico 2020-1 *** 
#--Trabajo realizado por: 
  # Monica A. Alvarez C.
  # Santiago Palacios L.
#   Ejercicio 2). Calculo de la raiz. 

Raiz <- function (n, E, x)
{
  eActual <- c() # Creacion de vectores para graficar. 
  eAnterior <- c()
  y = 0.5 * (x + (n / x)) # Primer calculo de la raiz.
  k = abs(x - y)
  while (k > E) { # Ciclo para iterar el calculo hasta alcanzar la tolerancia dada por el error. 
    x <- y
    y = 0.5 * (x + (n / x))
    eAnterior <- c(eAnterior, k)
    x =  k = abs(x - y)
    eActual <- c(eActual, k)
  }
  
  plot(eActual, # Grafica del ejercicio.
       eAnterior,
       main = "Error Actual vs Anterior",
       xlab = "Error Anterior",
       ylab = "Error Actual")
  points (eActual, eAnterior, col = "red")
  lines(eActual, eAnterior, col = "blue")
  return (cat("El resultado es: ", y, "con error de ", E))
}
Raiz(7, 0.00000001, 100)