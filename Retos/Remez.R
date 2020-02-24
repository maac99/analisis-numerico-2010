#     *** Analisis Numerico 2020-1 ***
#--Trabajo realizado por:
# Monica A. Alvarez C.
# Santiago Palacios L.
#   Metodo de Remez para la aproximacion a una funcion en un intervalo.

rm(list=ls())
require(PolynomF)
Puntos<-function (a,b,n){
  resultados =c()
  i=1
  while(i<n+1){
    aux=(0.5*(a+b)+(0.5*(b-a)*cos(((2*i-1)/(2*n))*pi)))
    resultados =c(resultados,round(aux,9))
    i=i+1
  }
  return(resultados)
}

x=Puntos((-pi/64),(pi/64),5) 
print(x)


1+X1+a2*(0^2)+a3*(0^3)=10
1+X2+a2*(1^2)+a3*(1^3)=15
1+X3+a2*(2^2)+a3*(2^3)=5


matrix = rbind(c(1,x[1],x[1]^2,x[1]^3,0.000001),c(1,x[2],x[2]^2,x[2]^3,-0.000001),c(1,x[3],x[3]^2,x[3]^3,0.000001),
               c(1,x[4],x[4]^2,x[4]^3,-0.000001), c(1,x[5],x[5]^2,x[5]^3,0.000001))

print (matrix)
b=c()
b=c(b,sin(x))
print(b)

coeficientes=solve(matrix,b)
coeficientes =round(coeficientes,12)
minimax=polynomial(coeficientes)

print(coeficientes)

valorr=polynom(a=coeficientes) # Estos los calculados
print(valorr(x))

original=b
print(valorr)
errorabs = abs(sin(x)-(valorr(x)))  
print(round(errorabs,16))
errorrela = (errorabs/sin(x))*100
print(round(errorrela,16))
