#Función a trabajar ln(x)
f <- function(x){
  return(log(x))
}

fx = 0
x = 0
i = 1

for (j in 1:10){
  evalua <- f(j)
  fx[i]=evalua
  x[i]=j
  i=i+1
}
Tab = data.frame(x,fx)
Tab

#Método de diferencias divididas
diffdiv <- c(0)
i = 1
diffdiv[i] = (fx[i+1]-fx[i])/(x[i+1]-x[i])
i=i+1

while(i<10){
  diffdiv[i] = (fx[i+1]-fx[i])/(x[i+1]-x[i]) 
  i=i+1
}
diffdiv

#Calcular error en el intervalo [1,2]
e <- abs((1.2-x[1])*(1.2-x[2])/factorial(2))

f = expression(log(x))

err<-function(f,n,s){
  x=0;
  while (x < n){
    f = D(f,'x');
    x = x+1;
  }
  return (s*abs(eval(f)));
}

cat("Error estimado", err(f,2,e))