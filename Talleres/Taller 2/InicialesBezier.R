#  *** Análisis numérico 2020 - 1 ***
# Trabajo realizado por:
#  Mónica A. Álvarez C.
#  Santiago Palacios L.
#  Paula C. Piñeros P.
# Aproximacion de iniciales mediante curvas de Bezier.
########################################
library(bezier)
p <- matrix(c(0,0.5,
              2.6,0.7,
              2.6,1.3,
              -0.6,2,
              0,2.8,
              2,2.8
), ncol=2, byrow=TRUE)

t <- seq(0,1, length=500)

bezier_points <- bezier(t, p)

bezier_points[,1]
bezier_points[,2]
plot(bezier_points[,1],
     bezier_points[,2],
     xlab = "X",
     ylab = "Y")
########################################
p <- matrix(c(0,0,
              0.4,1.2,
              1,1.2,
              1.2,-1.7,
              1.4,1.2,
              2,1.2,
              2.4,0

), ncol=2, byrow=TRUE)
t <- seq(0,1, length=500)

bezier_points <- bezier(t, p)

bezier_points[,1]
bezier_points[,2]
plot(bezier_points[,1],
     bezier_points[,2],
     xlab = "X",
     ylab = "Y")
#######################################
p <- matrix(c(2,   0.25,
              1.5, 0.25,
              1,  -0.5,
              0.5, 0.75,
              0,   1,
              0.5, 1.25,
              1,   2.5,
              1.5, 1.75,
              2,   1.75
              
), ncol=2, byrow=TRUE)
t <- seq(0,1, length=500)

bezier_points <- bezier(t, p)

plot(bezier_points[,1],
     bezier_points[,2],
     xlab = "X",
     ylab = "Y")
