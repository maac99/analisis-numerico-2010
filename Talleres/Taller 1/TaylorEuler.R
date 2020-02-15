#     *** Analisis Numerico 2020-1 *** 
#--Trabajo realizado por: 
# Monica A. Alvarez C.
# Santiago Palacios L.
#   Ejercicio 3). Calculo de e^n mediante el teorema de Taylor. 

TaylorEuler <- function(n)
{
  eActual <- c() # Creacion de vectores para graficar la aproximacion al valor final.
  eAnterior <- c()
  i = 0
  y = 0
  E=80000
  Actual=0
  Previa=0
  y = y + ((n ^ i) / factorial(i)) # Primer iteracion del polinomio de Taylor
  i = 1
  while (E > 0.00005) { #Ciclos para seguir el polinomio hasta la tolerancia definida. 
    
    eAnterior <- c(eAnterior, y)
    y = y + ((n ^ i) / factorial(i))
    eActual <- c(eActual, y)
    i = i + 1
    E=((tail(eActual,n=1)-tail(eAnterior,n=1))/tail(eActual,n=1))*100
  }
  
  plot(eActual,
       eAnterior,
       main = "Aproximaciones de Euler",
       xlab = "X",
       ylab = "Y")
  points (eActual, eAnterior, col = "red")
  lines(eActual, eAnterior, col = "blue")
  return (cat("El resultado es: ", round(y,4), "con error de ", E))
}
TaylorEuler(0.5)