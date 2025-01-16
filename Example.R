source("Utils_.R")

M=15
TT=20
theta=.01
beta=.9
K=3
mu=3
rho=0
phi=.8
P=20
Pcat=5
pGap=.1
pNAs=0.02
seed=1234

Y=generate_spatio_temporal_data(M, TT, theta, beta, K,
                                          mu,rho,phi,
                                          P,Pcat,seed,
                                          pGap,
                                          pNAs)

D=Y$dist_matrix
YY=Y$Y.NA

# Fit ST-JM
fit=STjumpDist(YY,3,
               D,
               jump_penalty=.1,
               spatial_penalty=.05,
               timeflag=T)

# Retrive true state sequence (S_true) and fitted state sequence (S_fit)
S_true=Y$S[unique(Y$Y.NA$t),]
S_fit=fit$best_s

# Order states for comparison
S_true_vec=order_states_condMean(YY$V6,S_true)
S_fit_vec=order_states_condMean(YY$V6,S_fit)

# Compute BAC
BAC_fit=BAC(S_true_vec,S_fit_vec);BAC_fit

S_fit_df=data.frame(S_fit)
colnames(S_fit_df)=1:M
S_fit_df$t=unique(Y$Y.NA$t)
S_fit_df$tt=1:dim(S_fit_df)[1]

coords=Y$spatial_points

# Dynamic plot of spatio-temporal clustering

P <- ncol(S_fit_df) - 2  

# Shiny UI
ui <- fluidPage(
  titlePanel("Spatio-Temporal Point Visualization"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("time", "Select Time:", 
                  min = min(S_fit_df$tt), 
                  max = max(S_fit_df$tt),
                  value = min(S_fit_df$tt),
                  step = 1,
                  animate = TRUE)
    ),
    
    mainPanel(
      plotOutput("spatioPlot")
    )
  )
)

# Shiny Server
server <- function(input, output) {
  output$spatioPlot <- renderPlot({
    # Estrai i dati per il tempo selezionato
    selected_data <- S_fit_df[S_fit_df$tt == input$time, 1:P]
    
    # Combina con le coordinate
    coords_selected <- data.frame(
      coords, 
      Category = as.factor(as.vector(t(selected_data)))  # Converte i dati selezionati in un vettore
    )
    
    # Crea il grafico
    ggplot(coords_selected, aes(x = x, y = y, color = Category)) +
      geom_point(size = 8) +
      scale_color_manual(values = c("1" = "lightgreen", "2" = "purple", "3" = "orange")) +
      labs(title = paste("Spatio-Temporal Visualization at Time", input$time),
           x = "X Coordinate", y = "Y Coordinate", color = "State") +
      theme_minimal()
  })
}

# Avvia l'app Shiny
shinyApp(ui = ui, server = server)

# Long format
S_fit_long=reshape2::melt(S_fit_df,id.vars="t")
colnames(S_fit_long)=c("t","m","State")

# Merge with original data
Y_res=merge(YY,S_fit_long,by=c("t","m"))

# Order by t and m
Y_res=Y_res[order(Y_res$t,Y_res$m),]

# State conditional prototypes
tapply(Y_res$V1,Y_res$State,Mode,na.rm=T)
tapply(Y_res$V6,Y_res$State,mean,na.rm=T)
