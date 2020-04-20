#  *** Análisis numérico 2020 - 1 ***
# Trabajo realizado por:
#  Mónica A. Álvarez C.
#  Santiago Palacios L.
#  Paula C. Piñeros P.
# Método de Curvas de Bezier para mortero valenciano.
library(TSdist)
library(bezier)
library(rgl)

funcion<- function(p,t,m,n){
  ps<-list()
  #ps[[1]]<-p[1:16,1:3]
  #ps[[2]]<-p[1:16,1:3]
  #ps[[3]]<-p[1:16,1:3]
  #ps[[4]]<-p[48:64,1:3]
  #print(p)
  print(m)
  print(n)
  i=0
  j=1
  while(j<m+1){
    ps[[j]]<-p[i:(i+n),1:3] 
    j=j+1
    i=i+n
  }
  xS<-vector()
  yS<-vector()
  zS<-vector()
  j=1
  while(j<m+1){
    
    bezier_points <- bezier(t, ps[[j]])   # El deg es el grado de la curvas. 
    # print("Aqui se jodio")
    # print(bezier_points)
    xS<-c(xS,bezier_points[,1])
    yS<-c(yS,bezier_points[,2])
    zS<-c(zS,bezier_points[,3])
    j=j+1
  }
  rgl2<-plot3d(xS,yS,zS,
               type="p",
               col="red",
               size=4,
               add=TRUE,
               xlab="x",ylab="y",zlab="z")
}
################################################
pintarpuntos<-function(p,x,y,t,m,n){
  print(p)
  print(x)
  p[,1]<-x
  p[,2]<-y
  xS<-p[,1]
  yS<-p[,2]
  zS<-p[,3]
  rgl1<-plot3d(xS,yS,zS,
               type="p",
               size=4,
               add=TRUE,
               xlab="x",ylab="y",zlab="z")
  print(p)
  if(m==1 &&n==1){
    print("Solo pintar")
    bezier_points <- bezier(t, p)   # El deg es el grado de la curvas. 
    # print("Aqui se jodio")
    # print(bezier_points)
    xS<-bezier_points[,1]
    yS<-bezier_points[,2]
    zS<-bezier_points[,3]
    rgl2<-plot3d(xS,yS,zS,
                 type="p",
                 col="red",
                 size=4,
                 add=TRUE,
                 xlab="x",ylab="y",zlab="z")
  }else{
     funcion(p,t,m,n)
  }
 
}
###########################################

    y1 = c(-0.551603, -0.546265, -0.540928, -0.457973, -0.4233915,-0.388881, -0.246651, -0.104422,  0.003964)
    x1 = c(0.0053400,  0.057126,  0.108912,  0.309168,  0.351259,  0.393350,  0.469368,  0.545386,  0.556060)
    y2 = c(-0.552894, -0.542210, -0.531527, -0.489773, -0.421175, -0.352577, -0.257753, -0.162930, -0.001724)
    x2 = c( 0.000414,  0.081042,  0.161671,  0.262261,  0.345809,  0.429357,  0.480249,  0.531141,  0.552892)
    x15=(x1+x1)/2
    y15=(y1+y2)/2
    x12=(x1+x15)/2
    y12=(y1+y15)/2
    x17=(x15+x2)/2
    y17=(y15+y2)/2
    y3 = c(-0.703701, -0.676505, -0.623363, -0.546265,  -0.439888,-0.333512, -0.207371,-0.1047825, -0.002194) 
    x3 = c( 0.000527,  0.205769,  0.333794,  0.448993,   0.536253, 0.623513,  0.676015,0.6898565,  0.703698) 
    x25=(x2+x3)/2
    y25=(y2+y3)/2
    x22=(x2+x25)/2
    y22=(y2+y25)/2
    x27=(x25+x3)/2
    y27=(y25+y3)/2
    y4 = c(-0.827465, -0.827404, -0.795485, -0.732997, -0.642340, -0.527669, -0.392169, -0.243842, -0.002580)
    x4 = c( 0.000619,  0.082118,  0.241958,  0.392501,  0.527960,  0.642578,  0.733174,  0.794910,  0.827461)
    x35=(x3+x4)/2
    y35=(y3+y4)/2
    x32=(x3+x35)/2
    y32=(y3+y35)/2
    x37=(x35+x4)/2
    y37=(y35+y4)/2
    y5 = c(-0.919430, -0.883896, -0.814463, -0.713730, -0.650022, -0.586314, -0.435755, -0.270943, -0.002867)
    x5 = c( 0.000688,  0.268849,  0.436123,  0.586637,  0.650316,  0.713995,  0.814660,  0.883257,  0.919426)
    x45=(x4+x5)/2
    y45=(y4+y5)/2
    x42=(x4+x45)/2
    y42=(y4+y45)/2
    x47=(x45+x5)/2
    y47=(y45+y5)/2
    y6 = c(-0.970937, -0.929325, -0.856322, -0.757691, -0.690059,-0.622428, -0.458151, -0.284868, -0.003030)    
    x6 = c( 0.000723,  0.282667,  0.458538,  0.622771,  0.690371,0.757972,  0.856529,  0.928652,  0.970933)
    x55=(x5+x6)/2
    y55=(y5+y6)/2
    x52=(x5+x55)/2
    y52=(y5+y55)/2
    x57=(x55+x6)/2
    y57=(y55+y6)/2
    y7 = c(-1.044518, -1.044502, -0.9888185, -0.880676,-0.772535, -0.634622, -0.467126, -0.290465, -0.003439)    
    x7 = c( 0.000599,  0.089114,  0.3620430,  0.498507, 0.634972,  0.772822,  0.873309,  0.951098,  1.044512)
    x65=(x6+x7)/2
    y65=(y6+y7)/2
    x62=(x6+x65)/2
    y62=(y6+y65)/2
    x67=(x65+x7)/2
    y67=(y65+y7)/2
    y8 = c(-1.242903, -1.242854, -1.124318, -1.008098, -0.939459, -0.861914, -0.757691, -0.622428, -0.461142, -0.287489, -0.029849, -0.004606)
    x8 = c( 0.000199,  0.025441,  0.046877,  0.120758,  0.285035,  0.461532,  0.622771,  0.757972,  0.862122,  0.938587,  1.242845,  1.242894)
    y8medio = c(-1.242903, -1.242854, -1.124318,  -0.939459, -0.861914, -0.757691,  -0.461142, -0.287489, -0.004606)
    x8medio = c( 0.000199,  0.025441,  0.046877,  0.285035,  0.461532,  0.622771,    0.862122,  0.938587,  1.242894)
    x75=(x7+x8medio)/2
    y75=(y7+y8medio)/2
    x72=(x7+x75)/2
    y72=(y7+y75)/2
    x77=(x65+x8medio)/2
    y77=(y65+y8medio)/2
    y9 = c(-1.124318, -1.008098, -0.939459, -0.861914, -0.757691, -0.622428, -0.461142, -0.287489, -0.050785)
    x9 = c( 0.046877,  0.120758,  0.285035,  0.461532,  0.622771,  0.757972,  0.862122,  0.938587,  1.123987)
    x85=(x8medio+x9)/2
    y85=(y8medio+y9)/2
    x82=(x8medio+x85)/2
    y82=(y8medio+y85)/2
    x87=(x85+x9)/2
    y87=(y85+y9)/2
    
    
    y10 = c(-0.004606, -0.050785)
    x10 = c( 1.242894,  1.123987)
    y11 = c(-1.242903, -1.124318)
    x11 = c( 0.000199,  0.046877)
    z10 = c( 2.4,  2.6) #primera de 8 y segunda de 9 

t <- seq(0,1, length=100)
m=3 #m es la cantidad de particiones 
puntos=6
n=puntos/m #n es la cantidad de puntos de cada particion 


p <- matrix(c(1,1,1), ncol=3, nrow=9, byrow=TRUE)
p12<- matrix(c(1.05,1.05,1.05), ncol=3, nrow=9, byrow=TRUE)
p15<- matrix(c(1.1,1.1,1.1), ncol=3, nrow=9, byrow=TRUE)
p17<- matrix(c(1.15,1.15,1.15), ncol=3, nrow=9, byrow=TRUE)

p2 <- matrix(c(1.2,1.2,1.2), ncol=3, nrow=9, byrow=TRUE)
p22<- matrix(c(1.25,1.25,1.25), ncol=3, nrow=9, byrow=TRUE)
p25<- matrix(c(1.3,1.3,1.3), ncol=3, nrow=9, byrow=TRUE)
p27<- matrix(c(1.35,1.35,1.35), ncol=3, nrow=9, byrow=TRUE)

p3 <- matrix(c(1.4,1.4,1.4), ncol=3, nrow=9, byrow=TRUE)
p32 <- matrix(c(1.45,1.45,1.45), ncol=3, nrow=9, byrow=TRUE)
p35<- matrix(c(1.5,1.5,1.5), ncol=3, nrow=9, byrow=TRUE)
p37<- matrix(c(1.55,1.55,1.55), ncol=3, nrow=9, byrow=TRUE)

p4 <- matrix(c(1.6,1.6,1.6), ncol=3, nrow=9, byrow=TRUE)
p42 <- matrix(c(1.65,1.65,1.65), ncol=3, nrow=9, byrow=TRUE)
p45<- matrix(c(1.7,1.7,1.7), ncol=3, nrow=9, byrow=TRUE)
p47<- matrix(c(1.75,1.75,1.75), ncol=3, nrow=9, byrow=TRUE)

p5 <- matrix(c(1.8,1.8,1.8), ncol=3, nrow=9, byrow=TRUE)
p52 <- matrix(c(1.85,1.85,1.85), ncol=3, nrow=9, byrow=TRUE)
p55<- matrix(c(1.9,1.9,1.9), ncol=3, nrow=9, byrow=TRUE)
p57<- matrix(c(1.95,1.95,1.95), ncol=3, nrow=9, byrow=TRUE)

p6 <- matrix(c(2,2,2), ncol=3, nrow=9, byrow=TRUE)
p62 <- matrix(c(2.05,2.05,2.05), ncol=3, nrow=9, byrow=TRUE)
p65<- matrix(c(2.1,2.1,2.1), ncol=3, nrow=9, byrow=TRUE)
p67<- matrix(c(2.15,2.15,2.15), ncol=3, nrow=9, byrow=TRUE)

p7 <- matrix(c(2.2,2.2,2.2), ncol=3, nrow=9, byrow=TRUE)
p72 <- matrix(c(2.25,2.25,2.25), ncol=3, nrow=9, byrow=TRUE)
p75<- matrix(c(2.3,2.3,2.3), ncol=3, nrow=9, byrow=TRUE)
p77<- matrix(c(2.35,2.35,2.35), ncol=3, nrow=9, byrow=TRUE)

p8 <- matrix(c(2.4,2.4,2.4), ncol=3, nrow=12, byrow=TRUE)#
p82 <- matrix(c(2.45,2.45,2.45), ncol=3, nrow=12, byrow=TRUE)#
p85<- matrix(c(2.5,2.5,2.5), ncol=3, nrow=9, byrow=TRUE)
p87<- matrix(c(2.55,2.55,2.55), ncol=3, nrow=9, byrow=TRUE)
p9 <- matrix(c(2.6,2.6,2.6), ncol=3, nrow=9, byrow=TRUE)

p10 <- matrix(c(0,0,0), ncol=3, nrow=2, byrow=TRUE)
p10[,3]<-z10
p11 <- matrix(c(0,0,0), ncol=3, nrow=2, byrow=TRUE)
p11[,3]<-z10



p <- matrix(c(1,1,1), ncol=3, nrow=9, byrow=TRUE)
p[,1]<-x1
p[,2]<-y1
xS<-p[,1]
yS<-p[,2]
zS<-p[,3]
rgl1<-plot3d(xS,yS,zS,
             type="p",
             size=4,
             xlab="x",ylab="y",zlab="z")#,xlim=c(0,1),ylim=c(-1,1) )

funcion(p,t,m,3)
pintartodo<-function(){
  pintarpuntos(p,x1,y1,t,m,3)
  pintarpuntos(p12,x12,y12,t,m,3)
  pintarpuntos(p15,x15,y15,t,m,3)
  pintarpuntos(p17,x17,y17,t,m,3)
  pintarpuntos(p2,x2,y2,t,m,3)
  pintarpuntos(p22,x22,y22,t,m,3)
  pintarpuntos(p25,x25,y25,t,m,3)
  pintarpuntos(p27,x27,y27,t,m,3)
  pintarpuntos(p3,x3,y3,t,m,3)
  pintarpuntos(p32,x32,y32,t,m,3)
  pintarpuntos(p35,x35,y35,t,m,3)
  pintarpuntos(p37,x37,y37,t,m,3)
  pintarpuntos(p4,x4,y4,t,m,3)
  pintarpuntos(p42,x42,y42,t,m,3)
  pintarpuntos(p45,x45,y45,t,m,3)
  pintarpuntos(p47,x47,y47,t,m,3)
  pintarpuntos(p5,x5,y5,t,m,3)
  pintarpuntos(p52,x52,y52,t,m,3)
  pintarpuntos(p55,x55,y55,t,m,3)
  pintarpuntos(p57,x57,y57,t,m,3)
  pintarpuntos(p6,x6,y6,t,m,3)
  pintarpuntos(p62,x62,y62,t,m,3)
  pintarpuntos(p65,x65,y65,t,m,3)
  pintarpuntos(p67,x67,y67,t,m,3)
  pintarpuntos(p7,x7,y7,t,m,3)
  pintarpuntos(p72,x72,y72,t,m,3)
  pintarpuntos(p75,x75,y75,t,m,3)
  pintarpuntos(p77,x77,y77,t,m,3)
  pintarpuntos(p8,x8,y8,t,4,3)
  #pintarpuntos(p82,x82,y82,t,m,3)
  pintarpuntos(p85,x85,y85,t,m,3)
  pintarpuntos(p87,x87,y87,t,m,3)
  pintarpuntos(p9,x9,y9,t,m,3)
  pintarpuntos(p10,x10,y10,t,1,1)
  pintarpuntos(p11,x11,y11,t,1,1)
  
  pintarpuntos(p,-x1,y1,t,m,3)
  pintarpuntos(p12,-x12,y12,t,m,3)
  pintarpuntos(p15,-x15,y15,t,m,3)
  pintarpuntos(p17,-x17,y17,t,m,3)
  pintarpuntos(p2,-x2,y2,t,m,3)
  pintarpuntos(p22,-x22,y22,t,m,3)
  pintarpuntos(p25,-x25,y25,t,m,3)
  pintarpuntos(p27,-x27,y27,t,m,3)
  pintarpuntos(p3,-x3,y3,t,m,3)
  pintarpuntos(p32,-x32,y32,t,m,3)
  pintarpuntos(p35,-x35,y35,t,m,3)
  pintarpuntos(p37,-x37,y37,t,m,3)
  pintarpuntos(p4,-x4,y4,t,m,3)
  pintarpuntos(p42,-x42,y42,t,m,3)
  pintarpuntos(p45,-x45,y45,t,m,3)
  pintarpuntos(p47,-x47,y47,t,m,3)
  pintarpuntos(p5,-x5,y5,t,m,3)
  pintarpuntos(p52,-x52,y52,t,m,3)
  pintarpuntos(p55,-x55,y55,t,m,3)
  pintarpuntos(p57,-x57,y57,t,m,3)
  pintarpuntos(p6,-x6,y6,t,m,3)
  pintarpuntos(p62,-x62,y62,t,m,3)
  pintarpuntos(p65,-x65,y65,t,m,3)
  pintarpuntos(p67,-x67,y67,t,m,3)
  pintarpuntos(p7, -x7,y7,t,m,3)
  pintarpuntos(p72,-x72,y72,t,m,3)
  pintarpuntos(p75,-x75,y75,t,m,3)
  pintarpuntos(p77,-x77,y77,t,m,3)
  pintarpuntos(p8, -x8,y8,t,4,3)
  #pintarpuntos(p82,-x82,y82,t,m,3)
  pintarpuntos(p85,-x85,y85,t,m,3)
  pintarpuntos(p87,-x87,y87,t,m,3)
  pintarpuntos(p9, -x9,y9,t,m,3)
  pintarpuntos(p10,-x10,y10,t,1,1)
  pintarpuntos(p11,-x11,y11,t,1,1)
  
  pintarpuntos(p,-x1,-y1,t,m,3)
  pintarpuntos(p12,-x12,-y12,t,m,3)
  pintarpuntos(p15,-x15,-y15,t,m,3)
  pintarpuntos(p17,-x17,-y17,t,m,3)
  pintarpuntos(p2,-x2,-y2,t,m,3)
  pintarpuntos(p22,-x22,-y22,t,m,3)
  pintarpuntos(p25,-x25,-y25,t,m,3)
  pintarpuntos(p27,-x27,-y27,t,m,3)
  pintarpuntos(p3,-x3,-y3,t,m,3)
  pintarpuntos(p32,-x32,-y32,t,m,3)
  pintarpuntos(p35,-x35,-y35,t,m,3)
  pintarpuntos(p37,-x37,-y37,t,m,3)
  pintarpuntos(p4,-x4,-y4,t,m,3)
  pintarpuntos(p42,-x42,-y42,t,m,3)
  pintarpuntos(p45,-x45,-y45,t,m,3)
  pintarpuntos(p47,-x47,-y47,t,m,3)
  pintarpuntos(p5,-x5,-y5,t,m,3)
  pintarpuntos(p52,-x52,-y52,t,m,3)
  pintarpuntos(p55,-x55,-y55,t,m,3)
  pintarpuntos(p57,-x57,-y57,t,m,3)
  pintarpuntos(p6,-x6,-y6,t,m,3)
  pintarpuntos(p62,-x62,-y62,t,m,3)
  pintarpuntos(p65,-x65,-y65,t,m,3)
  pintarpuntos(p67,-x67,-y67,t,m,3)
  pintarpuntos(p7,-x7,-y7,t,m,3)
  pintarpuntos(p72,-x72,-y72,t,m,3)
  pintarpuntos(p75,-x75,-y75,t,m,3)
  pintarpuntos(p77,-x77,-y77,t,m,3)
  pintarpuntos(p8,-x8,-y8,t,4,3)
  #pintarpuntos(p82,-x82,-y82,t,m,3)
  pintarpuntos(p85,-x85,-y85,t,m,3)
  pintarpuntos(p87,-x87,-y87,t,m,3)
  pintarpuntos(p9,-x9,-y9,t,m,3)
  pintarpuntos(p10,-x10,-y10,t,1,1)
  pintarpuntos(p11,-x11,-y11,t,1,1)
  
  pintarpuntos(p,x1,-y1,t,m,3)
  pintarpuntos(p12,x12,-y12,t,m,3)
  pintarpuntos(p15,x15,-y15,t,m,3)
  pintarpuntos(p17,x17,-y17,t,m,3)
  pintarpuntos(p2,x2,-y2,t,m,3)
  pintarpuntos(p22,x22,-y22,t,m,3)
  pintarpuntos(p25,x25,-y25,t,m,3)
  pintarpuntos(p27,x27,-y27,t,m,3)
  pintarpuntos(p3,x3,-y3,t,m,3)
  pintarpuntos(p32,x32,-y32,t,m,3)
  pintarpuntos(p35,x35,-y35,t,m,3)
  pintarpuntos(p37,x37,-y37,t,m,3)
  pintarpuntos(p4,x4,-y4,t,m,3)
  pintarpuntos(p42,x42,-y42,t,m,3)
  pintarpuntos(p45,x45,-y45,t,m,3)
  pintarpuntos(p47,x47,-y47,t,m,3)
  pintarpuntos(p5,x5,-y5,t,m,3)
  pintarpuntos(p52,x52,-y52,t,m,3)
  pintarpuntos(p55,x55,-y55,t,m,3)
  pintarpuntos(p57,x57,-y57,t,m,3)
  pintarpuntos(p6,x6,-y6,t,m,3)
  pintarpuntos(p62,x62,-y62,t,m,3)
  pintarpuntos(p65,x65,-y65,t,m,3)
  pintarpuntos(p67,x67,-y67,t,m,3)
  pintarpuntos(p7,x7,-y7,t,m,3)
  pintarpuntos(p72,x72,-y72,t,m,3)
  pintarpuntos(p75,x75,-y75,t,m,3)
  pintarpuntos(p77,x77,-y77,t,m,3)
  pintarpuntos(p8,x8,-y8,t,4,3)
 # pintarpuntos(p82,x82,-y82,t,m,3)
  pintarpuntos(p85,x85,-y85,t,m,3)
  pintarpuntos(p87,x87,-y87,t,m,3)
  pintarpuntos(p9,x9,-y9,t,m,3)
  pintarpuntos(p10,x10,-y10,t,1,1)
  pintarpuntos(p11,x11,-y11,t,1,1)
  

}

pintartodo()

########################x
