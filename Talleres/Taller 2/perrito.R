
#     *** Analisis Numerico 2020-1 *** 
#--Trabajo realizado por: 
# Monica A. Alvarez C.
# Santiago Palacios L.
# Paula C. Piñeros P.
#   Ejercicio 5.1). Grafica de la silueta de un perro.

#Cuerpo
x=c(1,   2,   5,   6,   7,  8.1,   10,  13, 17.6, 20, 23.5, 24.5, 25.5,  26.5, 27.5,  28.2,  29, 30)
y=c(3, 3.7, 3.9, 4.5, 5.7, 6.69, 7.12, 6.7, 4.45,  7,  6.1,  5.6, 5.87,  5.15,  4.1,   4.3, 4.1, 3)
x=x*1.2

#Patica trasera
xx=c(1, 5.3,   6, 5.7, 7.1, 8.8,  11,  12.2,  11.5,   10, 11.5, 8.3,   6)
xx=xx*1.2
yy=c(3, 2.9, 3.4, 2.5,   2, 2.3, 2.1,   2.4,     3,  3.4,  4.5, 5.7, 4.5)

#Oreja
xxx=c(20,   21, 19.8, 18.8, 17.3, 17.6)
xxx=xxx*1.2
yyy=c(  7, 4.45,    3,  2.4,  3.4, 4.45)

#Patica delantera
xxxx=c( 15,  17,  20, 24.5, 25.5, 26.5,  26, 24.4, 19.8)
xxxx=xxxx*1.2
yyyy=c(2.4, 2.9, 2.8,    3,  3.2,  2.9, 2.4,  2.4, 2.4)

#Unión1
xu=c(12.2, 15, 19.8)
xu=xu*1.2
yu=c(2.4, 2.4, 2.4)

#Unión2
xu2=c( 26, 30)
xu2=xu2*1.2
yu2=c(2.4, 3)


y1 = y[1:4] ; x1=x[1:4]
y2 = y[4:6] ; x2 = x[4:6]
y3 = y[6:9] ; x3 = x[6:9]
y4 = y[9:12] ; x4 = x[9:12]
y5 = y[12:15] ; x5 = x[12:15]
y6 = y[15:18] ; x6 = x[15:18]

y7 = yy[1:2] ; x7 = xx[1:2]
y8 = yy[2:3] ; x8 = xx[2:3]
y9 = yy[3:4] ; x9 = xx[3:4]
y10 = yy[4:6] ; x10 = xx[4:6]
y11 = yy[6:8] ; x11 = xx[6:8]
y12 = yy[8:10] ; x12 = xx[8:10]
y13 = yy[10:11] ; x13 = xx[10:11]
y14 = yy[11:13] ; x14 = xx[11:13]

y15 = yyy[1:2] ; x15 = xxx[1:2]
y16 = yyy[2:5] ; x16 = xxx[2:5]
y17 = yyy[5:6] ; x17 = xxx[5:6]

y18 = yyyy[1:2] ; x18 = xxxx[1:2]
y19 = yyyy[2:5] ; x19 = xxxx[2:5]
y20 = yyyy[5:6] ; x20 = xxxx[5:6]
y21 = yyyy[6:8] ; x21 = xxxx[6:8]
y22 = yyyy[8:9] ; x22 = xxxx[8:9]

y23 = yu[1:3]  ; x23 = xu[1:3]

y24 = yu2[1:2] ; x24 = xu2[1:2]


plot(x, # Grafica del ejercicio.
     y,
     main = "Silueta del perro",
     xlab = "X",
     ylab = "Y",
     xlim=c(0,36),
     ylim=c(1,8.5))
par(new = TRUE)
plot(xx,
     yy,
     xlab = "",
     ylab = "",
     xlim=c(0,36),
     ylim=c(1,8.5))
par(new = TRUE)
plot(xxx,
     yyy,
     xlab = "",
     ylab = "",
     xlim=c(0,36),
     ylim=c(1,8.5))
par(new = TRUE)
plot(xxxx,
     yyyy,
     xlab = "",
     ylab = "",
     xlim=c(0,36),
     ylim=c(1,8.5))
par(new = TRUE)
plot(xu,
     yu,
     xlab = "",
     ylab = "",
     xlim=c(0,36),
     ylim=c(1,8.5))
par(new = TRUE)
plot(xu2,
     yu2,
     xlab = "",
     ylab = "",
     xlim=c(0,36),
     ylim=c(1,8.5))

lines(spline(x1, y1, n = 201), col = "black")
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

lines(spline(x15, y15, n = 201), col = "orange4")
lines(spline(x16, y16, n = 201), col = "orange4")
lines(spline(x17, y17, n = 201), col = "orange4")

lines(spline(x18, y18, n = 201), col = "black")
lines(spline(x19, y19, n = 201), col = "black")
lines(spline(x20, y20, n = 201), col = "black")
lines(spline(x21, y21, n = 201), col = "black")
lines(spline(x22, y22, n = 201), col = "black")

lines(spline(x23, y23, n = 201), col = "black")

lines(spline(x24, y24, n = 201), col = "black")

