#     *** Analisis Numerico 2020-1 *** 
#--Trabajo realizado por: 
# Monica A. Alvarez C.
# Santiago Palacios L.
# Paula C. Piñeros P.
#   Ejercicio 5.2). Grafica de la silueta de la mano derecha.

x = c(-7, -6, -7.1, -8, - 11.5, -11.3, -9.5, -7, -4.2, -6.8, -5.3, -3.8, -1.8, -1.9,    0, 1.6,    3,  4,    6,  7,  7.1,    9, 13, 15.4,   16, 15.6, 15, 12,  8, 5.2, 5)
y = c( 0,  3,   11, 15,     22,  24.3,   24, 20, 16.8,   31,   33,   31,   21,   33, 35.2,  33, 19.2, 31, 32.8, 31,   19, 11.1, 17,  18, 17.5,   15, 12,  7,  2, 0.4, 0)
x = x * -1

y1 = y[29:31] ; x1 = x[29:31]
y2 = y[27:29] ; x2 = x[27:29]
y3 = y[25:27] ; x3 = x[25:27]
y4 = y[21:25] ; x4 = x[21:25]
y5 = y[20:21] ; x5 = x[20:21]
y6 = y[18:20] ; x6 = x[18:20]
y7 = y[17:18] ; x7 = x[17:18]
y8 = y[16:17] ; x8 = x[16:17]
y9 = y[14:16] ; x9 = x[14:16]
y10 = y[13:14] ; x10 = x[13:14]
y11 = y[10:13] ; x11 = x[10:13]
y12 = y[17:18] ; x12 = x[17:18]
y13 = y[16:17] ; x13 = x[16:17]
y14 = y[14:16] ; x14 = x[14:16]
y15 = y[13:14] ; x15 = x[13:14]
y16 = y[10:13] ; x16 = x[10:13]
y17 = y[9:10] ; x17 = x[9:10]
y18 = y[7:9] ; x18 = x[7:9]
y19 = y[6:7] ; x19 = x[6:7]
y20 = y[5:6] ; x20 = x[5:6]
y21 = y[4:5] ; x21 = x[4:5]
y22 = y[3:4] ; x22 = x[3:4]
y23 = y[2:3] ; x23 = x[2:3]
y24 = y[1:2] ; x24 = x[1:2]

       
plot(x,
     y,
     main = "Silueta mano derecha",
     xlab = "X",
     ylab = "Y",
     col = "red",
     xlim=c(-17, 12),
     ylim=c(  0, 36)
)

lines(spline(x1, y1, n = 201), col = "black")
lines(spline(x2, y2, n = 201), col = "black")
lines(spline(x3, y3, n = 201), col = "black")
lines(spline(x4, y4, n = 201), col = "black")
lines(spline(x5, y5, n = 201), col = "black")
lines(spline(x6, y6, n = 201), col = "black")
lines(spline(x7, y7, n = 201), col = "black")
lines(spline(x8, y8, n = 201), col = "black")
lines(spline(x9, y9, n = 201), col = "black")
lines(spline(x10, y10, n = 201), col = "black")
lines(spline(x11, y11, n = 201), col = "black")
lines(spline(x12, y12, n = 201), col = "black")
lines(spline(x13, y13, n = 201), col = "black")
lines(spline(x14, y14, n = 201), col = "black")
lines(spline(x15, y15, n = 201), col = "black")
lines(spline(x16, y16, n = 201), col = "black")
lines(spline(x17, y17, n = 201), col = "black")
lines(spline(x18, y18, n = 201), col = "black")
lines(spline(x19, y19, n = 201), col = "black")
lines(spline(x20, y20, n = 201), col = "black")
lines(spline(x21, y21, n = 201), col = "black")
lines(spline(x22, y22, n = 201), col = "black")
lines(spline(x23, y23, n = 201), col = "black")
lines(spline(x24, y24, n = 201), col = "black")
            