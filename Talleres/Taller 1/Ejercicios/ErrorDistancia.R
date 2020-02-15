#     * Analisis Numerico 2020-1 * 
#--Trabajo realizado por: 
# Monica A. Alvarez C.
# Santiago Palacios L.
#   Ejercicio 1). Calculo del error en la distancia.

distancia <- function(v,ev,t,et){
  d <- v*t
  ea<-(v*ev)+(t*et) #calculo del error absoluto
  erel<-(ev/v)+(et/t) #calculo del error relativo
  er<- erel*100
  rangosup<-d-ea
  rangoinf<-d+ea
  cat("Distancia:",d,"Rango: ", rangoinf, "-", rangosup, "\n")
  cat("Error relativo de:", er, "%")
}