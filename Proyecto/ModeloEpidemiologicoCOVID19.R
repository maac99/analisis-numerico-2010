
#     *** Analisis Numerico 2020-1 *** 
#--Trabajo realizado por: 
# Monica A. Alvarez C.
# Santiago Palacios L.
# Paula C. Piñeros P.
#   Reto 3. Transmisión de COVID 19.

library(shiny)
library(shinyjs)
library(deSolve)
ruta = file.choose()
database = read.csv(ruta)
dias=nrow(database)-1
print("help")
SeMuestran=TRUE
n = database[88,"dimessi_guariti"] + database[88,"deceduti"]+database[88,"totale_positivi"]

ui <- fluidPage(
    useShinyjs(),
    titlePanel("COVID-19"),

    sidebarLayout(
        sidebarPanel(
            sliderInput("poblacion",
                        "Población:",
                        min = 1,
                        max = n,
                        value = n ),
            sliderInput("infectados",
                        "Infectados:",
                        min = 1,
                        max = n ,
                        value = 1),
            sliderInput("recuperados",
                        "Recuperados:",
                        min = 0,
                        max = n,
                        value = 0),
            sliderInput("transmision",
                        "Tasa transmisión:",
                        min = 0,
                        max = 1,
                        value = 0.5),
            sliderInput("recuperacion",
                        "Tasa recuperación:",
                        min = 0,
                        max = 1,
                        value = 0.1),
            sliderInput("tiempodias",
                        "Tiempo en días:",
                        min = 1,
                        max = dias,
                        value = dias),
        ),
        
        mainPanel(
            actionButton("botonCalcularSI", "Calcular SI"),
            actionButton("botonCalcularSIR", "Calcular SIR"),
            actionButton("botonReales", "Mostrar Reales"),
            plotOutput("distPlot"),
            plotOutput("distPlot2"),
            plotOutput("distPlotERROR"),
            plotOutput("distPlot2ERROR"),
            plotOutput("distPlot3ERROR"),
        )
    )
)

server <- function(input, output, session) {
    
    
    observeEvent(input$botonCalcularSI, {
        print(class(input$poblacion))
        print(input$poblacion)
        CalcularSI()
    })
    
    observeEvent(input$botonCalcularSIR, {
        print(class(input$poblacion))
        print(input$poblacion)
        calcularSIR()
    })
    observeEvent(input$botonReales, {
        print(class(input$poblacion))
        print(input$poblacion)
        mostrarReales()
    })
    si <- function(times, iniciales, tasa)          #Función que calcula el modelo SI
    {
      with(as.list(c(iniciales, tasa)), 
           {
             dS <- - beta*S*I/(S+I)             #Ecuación diferencial para los susceptibles
             dI <- beta*S*I/(S+I)               #Ecuación diferencial para los infectados
             return(list(c(dS, dI)))
           })
    }
    
    sir <- function(times, iniciales, tasas)         #Función que calcula el modelo SIR
    {
      with(as.list(c(iniciales, tasas)), 
           {
             dS <- - beta*S*I/(S+I+R)           #Ecuación diferencial para los susceptibles
             dI <- beta*S*I/(S+I+R) - gamma*I   #Ecuación diferencial para los infectados
             dR <- gamma*I                      #Ecuación diferencial para los recuperados
             return(list(c(dS, dI, dR)))
           })
    }
    
    CalcularSI <- function(){
  
        iniciales <- c(S = input$poblacion,          #Iniciales contiene en S la población acorde al deslizador  
                       I = input$infectados)         #Iniciales contiene en I los infectados acorde al deslizador
      
        tasa<-c(beta = input$transmision)      #Tasa contiene en Beta la tasa de de transmisión acorde al deslizador

        output$distPlot3ERROR <- NULL
        
        output$distPlot <- renderPlot({         #Imprime la gráfica dada por la solución del método Runge Kutta 4
            tf = input$tiempodias
              
            
            times <- seq(0, tf, by = 0.1)
            simulacionM1SI.si <- as.data.frame(ode(y=iniciales, times=times, func=si,parms=tasa,method = "rk4"))
            attach(simulacionM1SI.si)

            plot(times, S, type="l", col="blue", ylim=c(0,sum(iniciales)), xlab="Tiempo (en días)", ylab="Población",main = "Metodo 1: Runge Kutta 4")
            lines(times, I, type="l", col="red")
            legend(x = "topright", legend=c("Susceptibles", "Infectados"), col=c("blue", "red"), lty=rep(1, 2)) 
        })
        
        output$distPlot2 <- renderPlot({        #Imprime la gráfica dada por la solución del método de Euler
            tf = input$tiempodias
            times <- seq(0, tf, by = 0.1)
            simulacionM2SI.si <- as.data.frame(ode(y=iniciales, times=times, func=si,parms=tasa,method = "euler"))
            attach(simulacionM2SI.si)
            plot(times, S, type="l", col="blue", ylim=c(0,sum(iniciales)), xlab="Tiempo (en días)", ylab="Población", main = "Metodo 2: Euler")
            lines(times, I, type="l", col="red")
            legend(x = "topright", legend=c("Susceptibles", "Infectados"), col=c("blue", "red"), lty=rep(1, 2)) 
            
        })
        
        output$distPlotERROR <- renderPlot({    #Imprime la gráfica del error del método de Runge Kutta 4 contra Euler
            tf = input$tiempodias
            times <- seq(0, tf, by = 0.1)
            simulacionM1SI.si <- as.data.frame(ode(y=iniciales, times=times, func=si,parms=tasa,method = "rk4"))
            simulacionM2SI.si <- as.data.frame(ode(y=iniciales, times=times, func=si,parms=tasa,method = "euler"))
            
            i =1
            erroresM1S <- c()
            for (x in simulacionM1SI.si$S) {
                erroresM1S <- c(erroresM1S, ( (abs(simulacionM2SI.si$S[i]-x) )/simulacionM2SI.si$S[i] ) *100)
                i = i +1;
            }
            x <- seq(1,tf)
            plot(simulacionM1SI.si$time , erroresM1S, col="blue", type="l", xlab="Tiempo (días)", ylab="Error relativo", main = "Error E vs RK4", ylim = c(0,100))
            
        })
        
        output$distPlot2ERROR <- renderPlot({    #Imprime la gráfica del error del método de Euler  contra Runge Kutta 4
            tf = input$tiempodias
            times <- seq(0, tf, by = 0.1)
            simulacionM1SI.si <- as.data.frame(ode(y=iniciales, times=times, func=si,parms=tasa,method = "rk4"))
            simulacionM2SI.si <- as.data.frame(ode(y=iniciales, times=times, func=si,parms=tasa,method = "euler"))
            i =1
            erroresM1I <- c()
            for (x in simulacionM1SI.si$I) {
                erroresM1I <- c(erroresM1I, ( (abs(simulacionM2SI.si$I[i]-x) )/simulacionM2SI.si$I[i] ) *100 )
                i = i +1;
            }
            
            x <- seq(1,tf)
            plot(simulacionM1SI.si$time , erroresM1I, col="red", type="l", xlab="Tiempo (en días)", ylab="Error relativo", main = "Error RK4 vs E", ylim = c(0,100))
            
        })
        
    }
    
    calcularSIR <- function(){
        if(SeMuestran){
            SeMuestran=TRUE
            
        }else{
            toggle("distPlotERROR")
            toggle("distPlot2")
            toggle("distPlot2ERROR")
            SeMuestran=TRUE
        }
      
        iniciales <- c(S = input$poblacion,         #Iniciales contiene en S la población acorde al deslizador
                       I = input$infectados,        #Iniciales contiene en I los infectados acorde al deslizador
                       R = input$recuperados)       #Iniciales contiene en R los recuperados acorde al deslizador
      
        tasas<-c(beta = input$transmision,      #Tasa contiene en beta la tasa de transmisión acorde al deslizador
                 gamma = input$recuperacion)    #Tasa contiene en gamma la tasa de recuperacón acorde al deslizador
      
        output$distPlot <- renderPlot({         #Imprime la gráfica dada por la solución del método Runge Kutta 4
            tf = input$tiempodias
            times <- seq(0, tf, by = 0.1)
            simulacionM1SI.sir <- as.data.frame(ode(y=iniciales, times=times, func=sir,parms=tasas,method = "rk4"))
            attach(simulacionM1SI.sir)
            
            plot(times, S, type="l", col="blue", ylim=c(0,sum(iniciales)), xlab="Tiempo (en días)", ylab="Población",main = "Metodo 1: Runge Kutta 4")
            lines(times, I, type="l", col="red")
            lines(times, R, type="l", col="green")
            legend(x = "topright", legend=c("Susceptibles", "Infectados", "Recuperados"), col=c("blue", "red", "green"), lty=rep(1, 2, 3)) 
        })
        
        output$distPlot2 <- renderPlot({        #Imprime la gráfica dada por la solución del método de Euler
            tf = input$tiempodias
            times <- seq(0, tf, by = 0.1)
             
            simulacionM2SI.sir <- as.data.frame(ode(y=iniciales, times=times, func=sir, parms=tasas, method = "euler"))
            attach(simulacionM2SI.sir)
            
            plot(times, S, type="l", col="blue", ylim=c(0,sum(iniciales)), xlab="Tiempo (en días)", ylab="Población", main = "Metodo 2: Euler")
            lines(times, I, type="l", col="red")
            lines(times, R, type="l", col="green")
            legend(x = "topright", legend=c("Susceptibles", "Infectados", "Recuperados"), col=c("blue", "red"), lty=rep(1, 2)) 
            
        })
        
        output$distPlotERROR <- renderPlot({    #Imprime la gráfica del error de los infectados del método de Runge kutta 4 contra la base de datos
          
            tf = input$tiempodias
            times <- seq(0, tf, by = 0.1)
            timesaux <- seq(0, (tf-2), by = 1)
            Steorico = database[["totale_positivi"]]       
            Iteorico = database[["totale_positivi"]]
            Rteorico = database[["dimessi_guariti"]]
            Mteorico = database[["deceduti"]] 
            
            n = database[88,"dimessi_guariti"] + database[88,"deceduti"]+database[88,"totale_positivi"]
            
            pobla <- -88:0; 
            pobla[pobla<0] <- n;
            Steorico=pobla-(Rteorico+Iteorico+Mteorico)
            Iteorico=Iteorico[-(tf:length(Iteorico))]
            TeoricoI <- Iteorico

            simulacionM2SI.sir <- as.data.frame(ode(y=iniciales, times=timesaux, func=sir,parms=tasas,method = "rk4"))
            
            i =1
            erroresM1I <- c()
            for (x in TeoricoI) { 
                erroresM1I <- c(erroresM1I, ( (abs(x-simulacionM2SI.sir$I[i]) )/x ) *100)
                i = i +1;
            }
            
            x <- seq(1,tf)
            plot(simulacionM2SI.sir$time , erroresM1I, col="blue", type="l", xlab="Dias", ylab="Error relativo", main = "Error población Infectada")
            
        })
        
        output$distPlot2ERROR <- renderPlot({   #Imprime la gráfica del error de la población susceptible del método de Runge kutta 4 contra la base de datos
            tf = input$tiempodias
            times <- seq(0, tf, by = 0.1)
            timesaux <- seq(0, (tf-2), by = 1)
            Steorico = database[["totale_positivi"]]         
            Iteorico = database[["totale_positivi"]]  
            Rteorico = database[["dimessi_guariti"]] 
            Mteorico = database[["deceduti"]] 
            
            n = database[88,"dimessi_guariti"] + database[88,"deceduti"]+database[88,"totale_positivi"]
            
            pobla <- -88:0; 
            pobla[pobla<0] <- n;
            Steorico=pobla-(Rteorico+Iteorico+Mteorico)
            Steorico=Steorico[-(tf:length(Steorico))]
            TeoricoS <- Steorico

            simulacionM2SI.sir <- as.data.frame(ode(y=iniciales, times=timesaux, func=sir,parms=tasas,method = "rk4"))
            i =1
            erroresM1S <- c()
            for (x in TeoricoS) {
                erroresM1S <- c(erroresM1S, ( (abs(x-simulacionM2SI.sir$S[i]) )/x ) *100 )
                i = i +1;
            }
            x <- seq(1,tf)
            plot(simulacionM2SI.sir$time , erroresM1S, col="red", type="l", xlab="Dias", ylab="Error relativo", main = "Error población Susceptible")
            
        })
        
        output$distPlot3ERROR <- renderPlot({   #Imprime la gráfica del error de los recuperados del método de Runge kutta 4 contra la base de datos
            tf = input$tiempodias
            times <- seq(0, tf, by = 0.1)
            timesaux <- seq(0, (tf-2), by = 1)
            Steorico = database[["totale_positivi"]]       
            Iteorico = database[["totale_positivi"]]
            Rteorico = database[["dimessi_guariti"]]
            Mteorico = database[["deceduti"]] 
            
            n = database[88,"dimessi_guariti"] + database[88,"deceduti"]+database[88,"totale_positivi"]
            
            pobla <- -88:0; 
            pobla[pobla<0] <- n;
            Steorico=pobla-(Rteorico+Iteorico+Mteorico)
            Rteorico=Rteorico[-(tf:length(Rteorico))]
            TeoricoR <- Rteorico

            simulacionM2SI.sir <- as.data.frame(ode(y=iniciales, times=timesaux, func=sir,parms=tasas,method = "rk4"))

            i =1
            erroresM1R <- c() 
            for (x in TeoricoR) {
                erroresM1R <- c(erroresM1R, ( (abs(x-simulacionM2SI.sir$R[i]) )/x ) *100 )
                i = i +1;
            }
            x <- seq(1,tf)
            plot(simulacionM2SI.sir$time , erroresM1R, col="red", type="l", xlab="Días", ylab="Error relativo", main = "Error población Recuperada")
            
        })
        
    }
    mostrarReales <- function(){
        output$distPlot2 <- NULL
        output$distPlotERROR <- NULL
        output$distPlot2ERROR <- NULL
        output$distPlot3ERROR <- NULL
        
        output$distPlot <- renderPlot({
            tf = input$tiempodias
            times <- seq(0+2, tf, by = 1)
            S = database[["totale_positivi"]] 
            I = database[["totale_positivi"]] 
            R = database[["dimessi_guariti"]] 
            M = database[["deceduti"]] 
            
            n = database[88,"dimessi_guariti"] + database[88,"deceduti"]+database[88,"totale_positivi"]
            
            pobla <- -88:0; 
            pobla[pobla<0] <- n;
            S=pobla-(R+I+M)
        
            S=S[-(tf:length(S))]
            R=R[-(tf:length(R))]
            I=I[-(tf:length(I))]
            M=M[-(tf:length(M))]

            plot(times, S, type="l", col="blue", xlab="Tiempo (en días)", ylab="Población",main = "Datos Reales", ylim=c(0,n))
            lines(times, I, type="l", col="red")
            lines(times, R, type="l", col="green")
            lines(times, M, type="l", col="black")
            legend(x = "topright", legend=c("Susceptibles", "Infectados", "Recuperados", "Muertos"), col=c("blue", "red", "green", "black"), lty=rep(1, 2, 3, 4)) 
        })
    
    }
}

shinyApp(ui = ui, server = server)
