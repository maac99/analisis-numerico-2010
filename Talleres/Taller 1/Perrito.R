#     *** Analisis Numerico 2020-1 *** 
#--Trabajo realizado por: 
# Monica A. Alvarez C.
# Santiago Palacios L.
#   Ejercicio 6). Grafica de la silueta de un perro.

  
y=c(3, 3.7, 3.9, 4.5, 5.7, 6.69, 7.12, 6.7, 4.45, 7, 6.1,  5.6, 5.87, 5.15, 4.1,  4.3, 4.1, 3)
 
x=c(1 , 2,   5,   6,  7, 8.1,   10,  13, 17.6, 20, 23.5, 24.5, 25.5,  26.5, 27.5,  28.2, 29, 30)
 x=x*1.2



y1 = y[1:4]   ; x1=x[1:4]
y2 = y[4:6]  ; x2 = x[4:6]
y3 = y[6:9] ; x3 = x[6:9]
y4 = y[9:12] ; x4 = x[9:12]
y5 = y[12:15] ; x5 = x[12:15]
y6 = y[15:18] ; x6 = x[15:18]

plot(x, # Grafica del ejercicio.
     y,
     main = "Silueta del perro",
     xlab = "X",
     ylab = "Y",
     xlim=c(0,36),
     ylim=c(1,8.5))

lines(spline(x1, y1, n = 201), col = "orange4")
lines(spline(x1, y1, n = 201), col = "orange4")
lines(spline(x2, y2, n = 201), col = "orange4")
lines(spline(x3, y3, n = 201), col = "orange4")
lines(spline(x4, y4, n = 201), col = "orange4")
lines(spline(x5,y5, n = 201), col =  "orange4")
lines(spline(x6, y6, n = 201), col = "orange4")
