#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(deSolve)
ruta = file.choose()
database = read.csv(ruta)
dias=nrow(database)-1

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Modelo SI"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("poblacion",
                        "Poblacion:",
                        min = 1,
                        max = 61000000,
                        value = 61000000),
            sliderInput("infectados",
                        "Infectados:",
                        min = 1,
                        max = 61000000,
                        value = 1),
            sliderInput("recuperados",
                        "Recuperados:",
                        min = 0,
                        max = 61000000,
                        value = 0),
            sliderInput("transmision",
                        "Tasa transmision:",
                        min = 0,
                        max = 1,
                        value = 0.5),
            sliderInput("recuperacion",
                        "Tasa recuperacion:",
                        min = 0,
                        max = 1,
                        value = 0.1),
            sliderInput("tiempodias",
                        "Tiempo en dias:",
                        min = 1,
                        max = dias,
                        value = dias),
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            actionButton("botonCalcularSI", "Calcular SI"),
            actionButton("botonCalcularSIR", "Calcular SIR"),
            plotOutput("distPlot"),
            plotOutput("distPlotERROR"),
            plotOutput("distPlot2"),
            plotOutput("distPlot2ERROR")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    
    observeEvent(input$botonCalcularSI, {
        print(class(input$poblacion))
        print(input$poblacion)
        CalcularEDOSI()
    })
    
    observeEvent(input$botonCalcularSIR, {
        print(class(input$poblacion))
        print(input$poblacion)
        calcularEDOSIR()
    })
    
    CalcularEDOSI <- function(){
        output$distPlot <- renderPlot({
            tf = input$tiempodias
            #ESTADO INICIAL, mover deslizadores solo cambia estos valores, lo iniciales. 
            init <- c(S = input$poblacion , # Susceptibles
                      I = input$infectados)  # Infectados

            param<-c(beta = input$transmision) #Hay que cambiar esto para que se ingrese con valores de los deslizadores. TOCA INVESTIGAR QUE VALORES SON MASO REALES. 
            #crear la funcion con las ODE.TOCA MIRAR SI ESTAN BIEN, Ellos interpretan gamma como probabilidad de infeccion y Wikipedia lo define como tiempo promedio de infeccion.
            si <- function(times, init, param) #Por eso toca mirar bien como se estan construyendo esas ED
            {
                with(as.list(c(init, param)), 
                     {
                         #ecuaciones diferenciales   
                         dS <- - beta*S*I/(S+I)
                         dI <- beta*S*I/(S+I)
                         #resultados de las tasas de cambio    
                         return(list(c(dS, dI)))
                     })
            }
            
            #intervalo de tiempo y resolucion
            times <- seq(0, tf, by = 0.1)
            #Mediante la funcion "ode" resolvemos el sistema de ecuaciones diferenciales y generamos un data frame 
            simulacionM1SI.si <- as.data.frame(ode(y=init, times=times, func=si,parms=param,method = "rk4"))      #Aqui se solucionan, mejor DEJAR ESTAS TRES LINEAS QUIETAS
            #Esta orden nos permite hacer referencia de manera directa a #las columnas de los resultados obtenidos
            attach(simulacionM1SI.si)

            #Representamos graficamente los resultados obtenidos CON ESTO SI PODEMOS JUGAR UN POCO MAS, HAY QUE CAMBIAR LOS TITULOS Y ESAS COSAS. 
            plot(times, S, type="l", col="blue", ylim=c(0,sum(init)), xlab="Tiempo (en días)", ylab="Numero de host",main = "Metodo 1: Runge Kutta 4")
            lines(times, I, type="l", col="red")
            #title("Modelo SI")
            legend(x = "topright", legend=c("Susceptibles", "Infectados"), col=c("blue", "red"), lty=rep(1, 2)) 
        })
        
        output$distPlot2 <- renderPlot({
            tf = input$tiempodias
            
            init <- c(S = input$poblacion , # Suceptibles 
                      I = input$infectados)  # Infectados

            param<-c(beta = input$transmision)
            #crear la funcion con las ODE
            si <- function(times, init, param) 
            {
                with(as.list(c(init, param)), 
                     {
                         #ecuaciones diferenciales   
                         dS <- -beta*S*I/(S+I)
                         dI <- beta*S*I/(S+I)
                         #resultados de las tasas de cambio    
                         return(list(c(dS, dI)))
                     })
            }
            
            #intervalo de tiempo y resolucion
            times <- seq(0, tf, by = 0.1)
            #Mediante la funcion "ode" resolvemos el sistema de ecuaciones diferenciales y generamos un data frame 
            simulacionM2SI.si <- as.data.frame(ode(y=init, times=times, func=si,parms=param,method = "euler"))
            #Esta orden nos permite hacer referencia de manera directa a #las columnas de los resultados obtenidos
            attach(simulacionM2SI.si)
            #Representamos graficamente los resultados obtenidos
            plot(times, S, type="l", col="blue", ylim=c(0,sum(init)), xlab="Tiempo (en días)", ylab="Numero de host", main = "Metodo 2: Euler")
            lines(times, I, type="l", col="red")
            #title("Modelo SI")
            legend(x = "topright", legend=c("Susceptibles", "Infectados"), col=c("blue", "red"), lty=rep(1, 2)) 
            
        })
        
        output$distPlotERROR <- renderPlot({
            tf = input$tiempodias
            init <- c(S = input$poblacion , # Suceptibles 
                      I = input$infectados)  # Infectados 
            
            param<-c(beta = input$transmision)
            #crear la funcion con las ODE
            si <- function(times, init, param) 
            {
                with(as.list(c(init, param)), 
                     {
                         #ecuaciones diferenciales   
                         dS <- - beta*S*I/(S+I)
                         dI <- beta*S*I/(S+I)
                         #resultados de las tasas de cambio    
                         return(list(c(dS, dI)))
                     })
            }
            
            #intervalo de tiempo y resolucion
            times <- seq(0, tf, by = 0.1)
            #Mediante la funcion "ode" resolvemos el sistema de ecuaciones diferenciales y generamos un data frame 
            simulacionM1SI.si <- as.data.frame(ode(y=init, times=times, func=si,parms=param,method = "rk4"))
            #simulacionM1SI.si
            
            #Mediante la funcion "ode" resolvemos el sistema de ecuaciones diferenciales y generamos un data frame 
            simulacionM2SI.si <- as.data.frame(ode(y=init, times=times, func=si,parms=param,method = "euler"))
            #simulacionM2SI.si
            
            # SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
            i =1
            erroresM1S <- c()
            for (x in simulacionM1SI.si$S) {
                #En M1, el teorico es M2
                erroresM1S <- c(erroresM1S, ( (abs(simulacionM2SI.si$S[i]-x) )/simulacionM2SI.si$S[i] ) *100)
                i = i +1;
            }
            x <- seq(1,70)
            plot(simulacionM1SI.si$time , erroresM1S, col="blue", type="l", xlab="Tiempo (en horas)", ylab="Error relativo", main = "Error", ylim = c(0,100))
            
        })
        
        output$distPlot2ERROR <- renderPlot({
            tf = input$tiempodias
            init <- c(S = input$poblacion , # Suceptibles
                      I = input$infectados)  # Infectados 
            
            #parametros del modelo (coeficientes de las variables)
            param <- c( beta = input$transmision) # Tasa de infeccion ( contagiados por unidad de tiempo) input$bethaSI
            #crear la funcion con las ODE
            si <- function(times, init, param) 
            {
                with(as.list(c(init, param)), 
                     {
                         #ecuaciones diferenciales   
                         dS <- -beta*S*I/(S+I)
                         dI <- beta*S*I/(S+I)
                         #resultados de las tasas de cambio    
                         return(list(c(dS, dI)))
                     })
            }
            
            #intervalo de tiempo y resolucion
            times <- seq(0, tf, by = 0.1)
            #Mediante la funcion "ode" resolvemos el sistema de ecuaciones diferenciales y generamos un data frame 
            simulacionM1SI.si <- as.data.frame(ode(y=init, times=times, func=si,parms=param,method = "rk4"))
            #simulacionM1SI.si
            
            #Mediante la funcion "ode" resolvemos el sistema de ecuaciones diferenciales y generamos un data frame 
            simulacionM2SI.si <- as.data.frame(ode(y=init, times=times, func=si,parms=param,method = "euler"))
            #simulacionM2SI.si
            
            i =1
            erroresM1I <- c()
            for (x in simulacionM1SI.si$I) {
                #En M1, el teorico es M2
                erroresM1I <- c(erroresM1I, ( (abs(simulacionM2SI.si$I[i]-x) )/simulacionM2SI.si$I[i] ) *100 )
                i = i +1;
            }
            
            x <- seq(1,701)
            plot(simulacionM1SI.si$time , erroresM1I, col="red", type="l", xlab="Tiempo (en días)", ylab="Error relativo", main = "Error", ylim = c(0,100))
            
        })
        
    }
    
    calcularEDOSIR <- function(){
        output$distPlot <- renderPlot({
            tf = input$tiempodias
            #ESTADO INICIAL, mover deslizadores solo cambia estos valores, lo iniciales. 
            init <- c(S = input$poblacion , # Susceptibles
                      I = input$infectados,  # Infectados
                      R = input$recuperados) # Recuperados
            
            param<-c(beta = input$transmision, #Hay que cambiar esto para que se ingrese con valores de los deslizadores. TOCA INVESTIGAR QUE VALORES SON MASO REALES. 
                     gamma = input$recuperacion)
            #crear la funcion con las ODE.TOCA MIRAR SI ESTAN BIEN, Ellos interpretan gamma como probabilidad de infeccion y Wikipedia lo define como tiempo promedio de infeccion.
            sir <- function(times, init, param) #Por eso toca mirar bien como se estan construyendo esas ED
            {
                with(as.list(c(init, param)), 
                     {
                         #ecuaciones diferenciales   
                         dS <- - beta*S*I/(S+I+R)
                         dI <- beta*S*I/(S+I+R) - gamma*I
                         dR <- gamma*I
                         #resultados de las tasas de cambio    
                         return(list(c(dS, dI, dR)))
                     })
            }
            
            #intervalo de tiempo y resolucion
            times <- seq(0, tf, by = 0.1)
            #Mediante la funcion "ode" resolvemos el sistema de ecuaciones diferenciales y generamos un data frame 
            simulacionM1SI.sir <- as.data.frame(ode(y=init, times=times, func=sir,parms=param,method = "rk4"))      #Aqui se solucionan, mejor DEJAR ESTAS TRES LINEAS QUIETAS
            #Esta orden nos permite hacer referencia de manera directa a #las columnas de los resultados obtenidos
            attach(simulacionM1SI.sir)
            
            #Representamos graficamente los resultados obtenidos CON ESTO SI PODEMOS JUGAR UN POCO MAS, HAY QUE CAMBIAR LOS TITULOS Y ESAS COSAS. 
            plot(times, S, type="l", col="blue", ylim=c(0,sum(init)), xlab="Tiempo (en días)", ylab="Numero de host",main = "Metodo 1: Runge Kutta 4")
            lines(times, I, type="l", col="red")
            lines(times, R, type="l", col="green")
            #title("Modelo SIR")
            legend(x = "topright", legend=c("Susceptibles", "Infectados", "Recuperados"), col=c("blue", "red", "green"), lty=rep(1, 2, 3)) 
        })
        
        output$distPlot2 <- renderPlot({
            tf = input$tiempodias
            
            init <- c(S = input$poblacion , # Suceptibles 
                      I = input$infectados,  # Infectados
                      R = input$recuperados) # Recuperados
            
            param<-c(beta = input$transmision,
                     gamma = input$recuperacion)
            #crear la funcion con las ODE
            sir <- function(times, init, param) 
            {
                with(as.list(c(init, param)), 
                     {
                         #ecuaciones diferenciales   
                         dS <- - beta*S*I/(S+I+R)
                         dI <- beta*S*I/(S+I+R) - gamma*I
                         dR <- gamma*I
                         #resultados de las tasas de cambio    
                         return(list(c(dS, dI, dR)))
                     })
            }
            
            #intervalo de tiempo y resolucion
            times <- seq(0, tf, by = 0.1)
            #Mediante la funcion "ode" resolvemos el sistema de ecuaciones diferenciales y generamos un data frame 
            simulacionM2SI.sir <- as.data.frame(ode(y=init, times=times, func=sir, parms=param, method = "euler"))
            #Esta orden nos permite hacer referencia de manera directa a #las columnas de los resultados obtenidos
            attach(simulacionM2SI.sir)
            #Representamos graficamente los resultados obtenidos
            plot(times, S, type="l", col="blue", ylim=c(0,sum(init)), xlab="Tiempo (en días)", ylab="Numero de host", main = "Metodo 2: Euler")
            lines(times, I, type="l", col="red")
            lines(times, R, type="l", col="green")
            #title("Modelo SIR")
            legend(x = "topright", legend=c("Susceptibles", "Infectados", "Recuperados"), col=c("blue", "red"), lty=rep(1, 2)) 
            
        })
        
        output$distPlotERROR <- renderPlot({
            tf = input$tiempodias
            init <- c(S = input$poblacion , # Suceptibles 
                      I = input$infectados,  # Infectados 
                      R = input$recuperados) # Recuperados
            
            param<-c(beta = input$transmision,
                     gamma = input$recuperacion)
            #crear la funcion con las ODE
            sir <- function(times, init, param) 
            {
                with(as.list(c(init, param)), 
                     {
                         #ecuaciones diferenciales   
                         dS <- - beta*S*I/(S+I+R)
                         dI <- beta*S*I/(S+I+R) - gamma*I
                         dR <- gamma*I
                         #resultados de las tasas de cambio    
                         return(list(c(dS, dI, dR)))
                     })
            }
            
            #intervalo de tiempo y resolucion
            times <- seq(0, tf, by = 0.1)
            #Mediante la funcion "ode" resolvemos el sistema de ecuaciones diferenciales y generamos un data frame 
            simulacionM1SI.sir <- as.data.frame(ode(y=init, times=times, func=sir,parms=param,method = "rk4"))
            #simulacionM1SI.sir
            
            #Mediante la funcion "ode" resolvemos el sistema de ecuaciones diferenciales y generamos un data frame 
            simulacionM2SI.sir <- as.data.frame(ode(y=init, times=times, func=sir,parms=param,method = "euler"))
            #simulacionM2SI.sir
            
            # SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
            i =1
            erroresM1S <- c()
            for (x in simulacionM1SI.sir$S) {
                #En M1, el teorico es M2
                erroresM1S <- c(erroresM1S, ( (abs(simulacionM2SI.sir$S[i]-x) )/simulacionM2SI.sir$S[i] ) *100)
                i = i +1;
            }
            x <- seq(1,70)
            plot(simulacionM1SI.sir$time , erroresM1S, col="blue", type="l", xlab="Tiempo (en días)", ylab="Error relativo", main = "Error", ylim = c(0,100))
            
        })
        
        output$distPlot2ERROR <- renderPlot({
            tf = input$tiempodias
            init <- c(S = input$poblacion , # Suceptibles
                      I = input$infectados,  # Infectados
                      R = input$recuperados) # Recuperados
            
            #parametros del modelo (coeficientes de las variables)
            param <- c( beta = input$transmision, # Tasa de infeccion ( contagiados por unidad de tiempo) input$bethaSI
                        gamma = input$recuperacion)
            #crear la funcion con las ODE
            sir <- function(times, init, param) 
            {
                with(as.list(c(init, param)), 
                     {
                         #ecuaciones diferenciales   
                         dS <- - beta*S*I/(S+I+R)
                         dI <- beta*S*I/(S+I+R) - gamma*I
                         dR <- gamma*I
                         #resultados de las tasas de cambio    
                         return(list(c(dS, dI, dR)))
                     })
            }
            
            #intervalo de tiempo y resolucion
            times <- seq(0, tf, by = 0.1)
            #Mediante la funcion "ode" resolvemos el sistema de ecuaciones diferenciales y generamos un data frame 
            simulacionM1SI.sir <- as.data.frame(ode(y=init, times=times, func=sir,parms=param,method = "rk4"))
            #simulacionM1SI.si
            
            #Mediante la funcion "ode" resolvemos el sistema de ecuaciones diferenciales y generamos un data frame 
            simulacionM2SI.sir <- as.data.frame(ode(y=init, times=times, func=sir,parms=param,method = "euler"))
            #simulacionM2SI.si
            
            i =1
            erroresM1I <- c()
            for (x in simulacionM1SI.sir$I) {
                #En M1, el teorico es M2
                erroresM1I <- c(erroresM1I, ( (abs(simulacionM2SI.sir$I[i]-x) )/simulacionM2SI.sir$I[i] ) *100 )
                i = i +1;
            }
            
            x <- seq(1,701)
            plot(simulacionM1SI.sir$time , erroresM1I, col="red", type="l", xlab="Tiempo (en días)", ylab="Error relativo", main = "Error", ylim = c(0,100))
            
        })
        
    }
}


# Run the application 
shinyApp(ui = ui, server = server)
