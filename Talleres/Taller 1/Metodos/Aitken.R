#     *** Analisis Numerico 2020-1 *** 
#--Trabajo realizado por: 
# Monica A. Alvarez C.
# Santiago Palacios L.
#   Metodo de Aitken con 20 valores.

Iniciales  <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
Datos <- c()

for(j in Iniciales){
  print(j)
  Datos<-c(Datos, cos(1/j))
  
}

x =Datos[1]
x1 =Datos[2]
x2 =Datos[3]
valor=c()

i=1
for(k in Datos[4:(length(Datos)-2)]){ # Ciclo que se repite por la cantidad de
  valor = c(valor ,( x2-((x2-x1)^2/x2-(2*x1)+x) ))   # terminos del polinomio propuesto.
  x =Datos[i]
  x1 =Datos[i+1]
  x2=Datos[i+2]
  i=i+1
  print(valor)
}

options(scipen=6)
