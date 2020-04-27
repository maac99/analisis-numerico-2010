#     * Analisis Numerico 2020-1 *
#--Trabajo realizado por:
# Monica A. Alvarez C.
# Santiago Palacios L.
# Paula Piñeros P.

library(pracma)
library(PolynomF)
library(Matrix)
#Datos
x <- c(100, 200, 300, 400, 500, 600)
y <- c(-160,-35, -4.2, 9.0, 16.9, 21.3)

# Polinomio de ajuste (polinomio interpolante en este caso)
datx <- x[1:5]
daty <-y[1:5]
polyAjuste = poly_calc(datx, daty)
polyAjuste
plot(x,y, pch=19, cex=1, col = "blue", asp=1, main = "Polinomio interpolante Lagrange" ) # Representación con puntos
curve(polyAjuste, add=T)

#Coeficiente virial

eval(polyAjuste(450))

#Polinomio de ajuste mediante Lagrange
