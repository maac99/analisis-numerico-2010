#     *** Analisis Numerico 2020-1 ***
#--Trabajo realizado por:
# Monica A. Alvarez C.
# Santiago Palacios L.
#   Metodo de Horner con numeros complejos y derivadas.


HornerDer <- function(coeficientes, x0) {
 # coeficientes=Cancelacion(coeficientes,0.60) # Esta linea se comenta o no para obtener cancelaciones.
  valor = coeficientes[length(coeficientes)] # Ajuste de las variables
  Der = 0
  multi = 0
  sumas = 0
  j = length(coeficientes) - 1
  while (j > 0) {
    i = coeficientes[j]
    #print(i)
    print(valor)
    Der = valor + (x0 * Der)
    valor <- x0 * valor + i                        
    multi = multi + 1       # En teoria si el coeficiente es 0 no se hace la multiplicacion, tener en cuenta
    sumas = sumas + 1
    j = j - 1
    
  }
  
  return(
    cat(
      "El valor de la derivada es:",Der,"\n El valor del polinomio es: ",valor,"; Multiplicaciones: ",multi,"; Sumas:  ", sumas
    )
  )
  
}

x0 = 10
coeficiente <- c(1,0.5,0.0012)
HornerDer(coeficiente, x0)


HornerComplejo<- function(coeficientes, x0){
  valor = coeficientes[1] # Ajuste de las variables 
  vbool=FALSE
  xbool=FALSE
  ibool=FALSE
  multi =0
  sumas=0
  for(i in coeficientes[2:length(coeficientes)]){ 
    if(Im(valor)>0)
    vbool=TRUE
    if(Im(x0)>0)
      xbool=TRUE
    if(Im(i)>0)
      ibool=TRUE
  
    
    valor <- x0*valor + i
    
    
    if(vbool && xbool && ibool){
      multi = multi + 4
      sumas = sumas + 4  
      print("a")
     
    }else
    {
      if(vbool && xbool){
        multi = multi + 4
        sumas = sumas + 3 
        print("b")
      }else{
        
        if(vbool ||xbool){
          multi = multi + 2
          sumas = sumas + 1 
          print("c")
         
        }else{
          multi = multi + 1
          sumas = sumas + 1 
          print("d")
          
        }
        
      }
    
    }
    
  }
  return(cat("El valor del polinomio es: ", round(valor,4), "\nEl total de operaciones es de: ", 
             multi+sumas,"; Multiplicaciones: ",multi,"; Sumas:  ",sumas))
 
}
x0 = 2+3i
coeficiente <- c(2,2+3i,1+1i,3,0,12)

HornerComplejo(coeficiente, x0)
