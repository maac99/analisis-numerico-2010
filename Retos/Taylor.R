#     *** Analisis Numerico 2020-1 ***
#--Trabajo realizado por:
# Monica A. Alvarez C.
# Santiago Palacios L.
# Aproximacion a una funcion en un intervalo con el Teorema de Taylor.


require(pracma)
t=seq(-pi/64,pi/64, ((2*(pi/64))/5)) # Solo tienes que cambiar aqui, la cantidad que coloques pues van a ser la de filas de la tabla, pon mas de 10 pero no te pases de verga
y=sin(t)
plot(t,y,type="l", xlab="time", ylab="Sine wave")
print(t)
segmentos=seq((-pi/64),(pi/64),((2*(pi/64))/5)) # Y aqui, tienen que se el mismo valor, evidentemente

Taylor<-function(segmentos){
resultados=c()
aux=0
for(i in segmentos){
  #aux=i+(-(i^3)/6)+((i^5)/120) #Tercer Grado
  aux=i+(-(i^3)/6) #Segundo Grado
  #aux=i  #Primer Grado
  resultados=c(resultados, aux)
}
print (resultados)

plot(segmentos,resultados,type="l", xlab="time", ylab="Taylor")
return (resultados)
}

# Los resultados los sacas de aqui
print(y) # Estos son los valores reales
valorr=Taylor(segmentos) # Estos los calculados


#Calculo de los errores
original=y
errorabs = abs(original-(valorr))  
print(round(errorabs,9))
errorrela = (errorabs/valorr)*100
print(round(errorrela,9))

