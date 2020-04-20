#  *** Análisis numérico 2020 - 1 ***
# Trabajo realizado por:
#  Mónica A. Álvarez C.
#  Santiago Palacios L.
#  Paula C. Piñeros P.
# Método de splines para mortero valenciano.

library(rgl)

#Imprime en plot 3D
imprimir = function(x,y,z){
  xa = x*-1
  ya = y*-1
  rgl1<-plot3d(x,y,z,
               type="l",
               size=4,
               add=TRUE,
               xlab="x",ylab="y",zlab="z",
               col = "blue")
  rgl1<-plot3d(xa,ya,z,
               type="l",
               size=4,
               add=TRUE,
               xlab="x",ylab="y",zlab="z",
               col = "blue")
  rgl1<-plot3d(xa,y,z,
               type="l",
               size=4,
               add=TRUE,
               xlab="x",ylab="y",zlab="z",
               col = "blue")
  rgl1<-plot3d(x,ya,z,
               type="l",
               size=4,
               add=TRUE,
               xlab="x",ylab="y",zlab="z",
               col = "blue")
}
#Imprime en 2D puntos los puntos x,y
pl = function(x,y){
  plot(x,
       y,
       xlab = "X",
       ylab = "Y",
       xlim=c(0,1.3),
       ylim=c(-1.3,1.3))
}
#Clcula error por cada spline
error = function(x,y){
  erro = c()
  interpolado = splinefun(x,y)
  for(i in 1:length(x)){
    yExp = interpolado(x)
    err = abs((yExp - y)/y)
    err = err * 100
    erro = c(erro, err)
  }
  return (erro)
}

#1 Baja
y1 = c(-0.551603, -0.540928, -0.457973, -0.388881, -0.104422, 0.003964)
x1 = c(0.00534,    0.108912,  0.309168,  0.393350,  0.545386, 0.556060)

y11 = y1[1:3] ; x11 = x1[1:3]
y21 = y1[3:5] ; x21 = x1[3:5]
y31 = y1[5:6] ; x31 = x1[5:6]

vecx = c(x11, x21, x31)
vecy = c(y11, y21, y31)

z1 = -0.879887
pl(x1, y1)

a11 = spline(x11, y11, n = 201)
a12 = spline(x21, y21, n = 201)
a13 = spline(x31, y31, n = 201)

a111 = a11[1]; a121 = a12[1]; a131 = a13[1]
a111 = unlist(a111, use.names=FALSE); a121 = unlist(a121, use.names=FALSE); a131 = unlist(a131, use.names=FALSE); 
xa1 = c(a111,a121,a131)
a112 = a11[2]; a122 = a12[2]; a132 = a13[2]
a112 = unlist(a112, use.names=FALSE); a122 = unlist(a122, use.names=FALSE); a132 = unlist(a132, use.names=FALSE); 
ya1 = c(a112,a122,a132)

lines(a11, col = "orange4")
lines(a12, col = "orange4")
lines(a13, col = "orange4")

#2Baja
y2 = c(-0.552894, -0.531527,  -0.489773, -0.352577,  -0.162930, -0.001724)
x2 = c( 0.000414,  0.161671,   0.262261,  0.429357,   0.531141,  0.552892)

y12 = y2[1:3] ; x12 = x2[1:3]
y22 = y2[3:5] ; x22 = x2[3:5]
y32 = y2[5:6] ; x32 = x2[5:6]

vecx = c(vecx, x12, x22, x32)
vecy = c(vecy, y12, y22, y32)

z2 = -0.831470
pl(x2,y2)

a21 = spline(x12, y12, n = 201)
a22 = spline(x22, y22, n = 201)
a23 = spline(x32, y32, n = 201)

a211 = a21[1]; a221 = a22[1]; a231 = a23[1]
a211 = unlist(a211, use.names=FALSE); a221 = unlist(a221, use.names=FALSE); a231 = unlist(a231, use.names=FALSE); 
xa2 = c(a211,a221,a231)
a212 = a21[2]; a222 = a22[2]; a232 = a23[2]
a212 = unlist(a212, use.names=FALSE); a222 = unlist(a222, use.names=FALSE); a232 = unlist(a232, use.names=FALSE); 
ya2 = c(a212,a222,a232)

lines(a21, col = "orange4")
lines(a22, col = "orange4")
lines(a23, col = "orange4")

#3Baja
y3 = c(-0.703701, -0.676505, -0.623363, -0.546265, -0.333512, -0.207371, -0.002194)
x3 = c( 0.000527,  0.205769,  0.333794,  0.448993,  0.623513,  0.676015,  0.703698)

y13 = y3[1:3] ; x13 = x3[1:3]
y23 = y3[3:6] ; x23 = x3[3:6]
y33 = y3[6:7] ; x33 = x3[6:7]

vecx = c(vecx, x13, x23, x33)
vecy = c(vecy, y13, y23, y33)

z3 = -0.707107
pl(x3,y3)

a31 = spline(x13, y13, n = 201)
a32 = spline(x23, y23, n = 201)
a33 = spline(x33, y33, n = 201)

xa3 = c(a31[1], a32[1], a33[1])
ya3 = c(a31[2], a32[2], a33[2])

lines(a31, col = "orange4")
lines(a32, col = "orange4")
lines(a33, col = "orange4")

a311 = a31[1]; a321 = a32[1]; a331 = a33[1]
a311 = unlist(a311, use.names=FALSE); a321 = unlist(a321, use.names=FALSE); a331 = unlist(a331, use.names=FALSE); 
xa3 = c(a311,a321,a331)
a312 = a31[2]; a322 = a32[2]; a332 = a33[2]
a312 = unlist(a312, use.names=FALSE); a322 = unlist(a322, use.names=FALSE); a332 = unlist(a332, use.names=FALSE); 
ya3 = c(a312,a322,a332)

#4Baja
y4 = c(-0.827465, -0.827404, -0.795485, -0.732997, -0.642340, -0.527669, -0.392169, -0.243842, -0.002580)
x4 = c( 0.000619,  0.082118,  0.241958,  0.392501,  0.527960,  0.642578,  0.733174,  0.794910,  0.827461)

y14 = y4[1:3] ; x14 = x4[1:3]
y24 = y4[3:6] ; x24 = x4[3:6]
y34 = y4[6:8] ; x34 = x4[6:8]
y44 = y4[8:9] ; x44 = x4[8:9] 

vecx = c(vecx, x14, x24, x34, x44)
vecy = c(vecy, y14, y24, y34, x44)

z4 = -0.555570
pl(x4,y4)

a41 = spline(x14, y14, n = 201)
a42 = spline(x24, y24, n = 201)
a43 = spline(x34, y34, n = 201)
a44 = spline(x44, y44, n = 201)

xa4 = c(a41[1], a42[1], a43[1], a44[1])
ya4 = c(a41[2], a42[2], a43[2], a44[2])

lines(a41, col = "orange4")
lines(a42, col = "orange4")
lines(a43, col = "orange4")
lines(a44, col = "orange4")

a411 = a41[1]; a421 = a42[1]; a431 = a43[1]; a441 = a44[1]
a411 = unlist(a411, use.names=FALSE); a421 = unlist(a421, use.names=FALSE); a431 = unlist(a431, use.names=FALSE); a441 = unlist(a441, use.names=FALSE)
xa4 = c(a411,a421,a431,a441)
a412 = a41[2]; a422 = a42[2]; a432 = a43[2]; a442 = a44[2]
a412 = unlist(a412, use.names=FALSE); a422 = unlist(a422, use.names=FALSE); a432 = unlist(a432, use.names=FALSE); a442 = unlist(a442, use.names=FALSE) 
ya4 = c(a412,a422,a432,a442)

#5Baja
y5 = c(-0.919430, -0.883896, -0.814463, -0.713730, -0.586314, -0.435755, -0.270943, -0.002867)
x5 = c( 0.000688,  0.268849,  0.436123,  0.586637,  0.713995,  0.814660,  0.883257,  0.919426)

y15 = y5[1:3] ; x15 = x5[1:3]
y25 = y5[3:6] ; x25 = x5[3:6]
y35 = y5[6:7] ; x35 = x5[6:7]
y45 = y5[7:8] ; x45 = x5[7:8] 

vecx = c(vecx, x15, x25, x35, x45)
vecy = c(vecy, y15, y25, y35, y45)

z5 = -0.382683
pl(x5,y5)

a51 = spline(x15, y15, n = 201)
a52 = spline(x25, y25, n = 201)
a53 = spline(x35, y35, n = 201)
a54 = spline(x45, y45, n = 201)

xa5 = c(a51[1], a52[1], a53[1], a54[1])
ya5 = c(a51[2], a52[2], a53[2], a54[2])

lines(a51, col = "orange4")
lines(a52, col = "orange4")
lines(a53, col = "orange4")
lines(a54, col = "orange4")

a511 = a51[1]; a521 = a52[1]; a531 = a53[1]; a541 = a54[1]
a511 = unlist(a511, use.names=FALSE); a521 = unlist(a521, use.names=FALSE); a531 = unlist(a531, use.names=FALSE); a541 = unlist(a541, use.names=FALSE)
xa5 = c(a511,a521,a531,a541)
a512 = a51[2]; a522 = a52[2]; a532 = a53[2]; a542 = a54[2]
a512 = unlist(a512, use.names=FALSE); a522 = unlist(a522, use.names=FALSE); a532 = unlist(a532, use.names=FALSE); a542 = unlist(a542, use.names=FALSE) 
ya5 = c(a512,a522,a532,a542)

#6Baja
y6 = c(-0.970937, -0.929325, -0.856322, -0.757691, -0.622428, -0.458151, -0.284868, -0.003030)    
x6 = c( 0.000723,  0.282667,  0.458538,  0.622771,  0.757972,  0.856529,  0.928652,  0.970933)    

y16 = y6[1:3] ; x16 = x6[1:3]
y26 = y6[3:7] ; x26 = x6[3:7]
y36 = y6[7:8] ; x36 = x6[7:8]

vecx = c(vecx, x16, x26, x36)
vecy = c(vecy, y16, y26, y36)

z6 = -0.199724
pl(x6,y6)

a61 = spline(x16, y16, n = 201)
a62 = spline(x26, y26, n = 201)
a63 = spline(x36, y36, n = 201)

xa6 = c(a61[1], a62[1], a62[1])
ya6 = c(a61[2], a62[2], a62[2])

lines(a61, col = "orange4")
lines(a62, col = "orange4")
lines(a63, col = "orange4")

a611 = a61[1]; a621 = a62[1]; a631 = a63[1]
a611 = unlist(a611, use.names=FALSE); a621 = unlist(a621, use.names=FALSE); a631 = unlist(a631, use.names=FALSE); 
xa6 = c(a611,a621,a631)
a612 = a61[2]; a622 = a62[2]; a632 = a63[2]
a612 = unlist(a612, use.names=FALSE); a622 = unlist(a622, use.names=FALSE); a632 = unlist(a632, use.names=FALSE); 
ya6 = c(a612,a622,a632)

#7Baja
y7 = c(-1.044518, -1.044502, -0.9888185, -0.772535, -0.634622, -0.467126, -0.290465, -0.003439)    
x7 = c( 0.000599,  0.089114,  0.3620430,  0.634972,  0.772822,  0.873309,  0.951098,  1.044512)    

y17 = y7[1:3] ; x17 = x7[1:3]
y27 = y7[3:6] ; x27 = x7[3:6]
y37 = y7[6:8] ; x37 = x7[6:8]

vecx = c(vecx, x17, x27, x37)
vecy = c(vecy, y17, y27, y37)

z7 = -0.019903
pl(x7,y7)

a71 = spline(x17, y17, n = 201)
a72 = spline(x27, y27, n = 201)
a73 = spline(x37, y37, n = 201)

lines(a71, col = "orange4")
lines(a72, col = "orange4")
lines(a73, col = "orange4")

a711 = a71[1]; a721 = a72[1]; a731 = a73[1]
a711 = unlist(a711, use.names=FALSE); a721 = unlist(a721, use.names=FALSE); a731 = unlist(a731, use.names=FALSE); 
xa7 = c(a711,a721,a731)
a712 = a71[2]; a722 = a72[2]; a732 = a73[2]
a712 = unlist(a712, use.names=FALSE); a722 = unlist(a722, use.names=FALSE); a732 = unlist(a732, use.names=FALSE); 
ya7 = c(a712,a722,a732)

#8
y8 = c(-1.242903, -1.124318, -0.939459, -0.861914, -0.757691, -0.622428, -0.461142, -0.287489, -0.004606)
x8 = c( 0.000199,  0.046877,  0.285035,  0.461532,  0.622771,  0.757972,  0.862122,  0.938587,  1.242894)

y18 = y8[1:4] ; x18 = x8[1:4]
y28 = y8[4:7] ; x28 = x8[4:7]
y38 = y8[7:9] ; x38 = x8[7:9]

vecx = c(vecx, x18, x28, x38)
vecy = c(vecy, y18, y28, y38)

z8 = 0.08604
pl(x8,y8)

a81 = spline(x18, y18, n = 201)
a82 = spline(x28, y28, n = 201)
a83 = spline(x38, y38, n = 201)

lines(a81, col = "orange4")
lines(a82, col = "orange4")
lines(a83, col = "orange4")

a811 = a81[1]; a821 = a82[1]; a831 = a83[1]
a811 = unlist(a811, use.names=FALSE); a821 = unlist(a821, use.names=FALSE); a831 = unlist(a831, use.names=FALSE); 
xa8 = c(a811,a821,a831)
a812 = a81[2]; a822 = a82[2]; a832 = a83[2]
a812 = unlist(a812, use.names=FALSE); a822 = unlist(a822, use.names=FALSE); a832 = unlist(a832, use.names=FALSE); 
ya8 = c(a812,a822,a832)

#9
y9 = c(-1.124318, -0.939459, -0.861914, -0.757691, -0.622428, -0.461142, -0.287489, -0.050785)
x9 = c( 0.046877,  0.285035,  0.461532,  0.622771,  0.757972,  0.862122,  0.938587,  1.123987)

y19 = y9[1:3] ; x19 = x9[1:3]
y29 = y9[3:7] ; x29 = x9[3:7]
y39 = y9[7:8] ; x39 = x9[7:8]

vecx = c(vecx, x19, x29, x39)
vecy = c(vecy, y19, y29, y39)

z9 = 0.162575
pl(x9,y9)

a91 = spline(x19, y19, n = 201)
a92 = spline(x29, y29, n = 201)
a93 = spline(x39, y39, n = 201)

lines(spline(x19, y19, n = 201), col = "orange4")
lines(spline(x29, y29, n = 201), col = "orange4")
lines(spline(x39, y39, n = 201), col = "orange4")

a911 = a91[1]; a921 = a92[1]; a931 = a93[1]
a911 = unlist(a911, use.names=FALSE); a921 = unlist(a921, use.names=FALSE); a931 = unlist(a931, use.names=FALSE); 
xa9 = c(a911,a921,a931)
a912 = a91[2]; a922 = a92[2]; a932 = a93[2]
a912 = unlist(a912, use.names=FALSE); a922 = unlist(a922, use.names=FALSE); a932 = unlist(a932, use.names=FALSE); 
ya9 = c(a912,a922,a932)


#10
y10 = c(-0.004606, -0.050785)
x10 = c( 1.242894,  1.123987)
z10 = c( 0.086040,  0.162575)

#11
y11 = c(-1.242903, -1.124318)
x11 = c( 0.000199,  0.046877)
z11 = c( 0.086040,  0.162575)


imprimir(xa1,ya1,z1)
imprimir(xa2,ya2,z2)
imprimir(xa3,ya3,z3)
imprimir(xa4,ya4,z4)
imprimir(xa5,ya5,z5)
imprimir(xa6,ya6,z6)
imprimir(xa7,ya7,z7)
imprimir(xa8,ya8,z8)
imprimir(xa9,ya9,z9)
imprimir(x10,y10,z10)
imprimir(x11,y11,z11)

errores = c()
errores = c(error(x11, y11),error(x21, y21),error(x31, y31))
errores = c(errores,error(x12, y12),error(x22, y22),error(x32, y32))
errores = c(errores,error(x13, y13),error(x23, y23),error(x33, y33))
errores = c(errores,error(x14, y14),error(x24, y24),error(x34, y34),error(x44, y44))
errores = c(errores,error(x15, y15),error(x25, y25),error(x35, y35),error(x45, y45))
errores = c(errores,error(x16, y16),error(x26, y26),error(x36, y36))
errores = c(errores,error(x17, y17),error(x27, y27),error(x37, y37))
errores = c(errores,error(x18, y18),error(x28, y28),error(x38, y38))
errores = c(errores,error(x19, y19),error(x29, y29),error(x39, y39))

errorProm = 0
for(i in 1:length(errores)){
  errorProm = errorProm + errores[i]
}
errorProm = errorProm/length(errores)

print(errorProm)

