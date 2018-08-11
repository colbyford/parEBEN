library(EBEN)
library(parEBEN)

## Register your local cluster.
library(doParallel)  
no_cores <- detectCores() 
cl <- makeCluster(no_cores)
#clusterExport(cl, c("parEBEN.cv.doParallel"))
registerDoParallel(cl)

## Generate your test matrix
#tests <- data.frame( n = c(25,25,50,50),
#                     k = c(100,200,100,200),
#                     nFolds = c(3,3,3,3),
#                     ser_user_time = c(0,0,0,0),
#                     ser_system_time = c(0,0,0,0),
#                     ser_elapsed_time = c(0,0,0,0),
#                     par_user_time = c(0,0,0,0),
#                     par_system_time = c(0,0,0,0),
#                     par_elapsed_time = c(0,0,0,0)
#                     )

tests <- read.csv("testmatrix_2-7-2018_1_gaussian.csv")

## Loop throught the tests in the matrix and report the time
for (i in 1:nrow(tests)){
  data(BASIS)
  data(y)
  n <- tests[i,1]
  k <- tests[i,2]
  nFolds <- tests[i,3]
  #N <- length(y)
  set.seed(1)
  #set <- sample(N,n)
  BASISset <- BASIS[1:n,1:k]
  yset  <- y[1:n]
  cat("Iteration:",i ," (n:", n,"k:", k,"nFolds:", nFolds,")\n")
  sertime <- system.time(EBelasticNet.GaussianCV(BASIS,y, nFolds, Epis = "no"))
  tests[i,4] <- sertime[1]
  tests[i,5] <- sertime[2]
  tests[i,6] <- sertime[3]
  partime <- system.time(parEBEN.cv(BASIS, y, nFolds, Epis = "no", parMethod = "doParallel", prior = "gaussian"))
  tests[i,7] <- partime[1]
  tests[i,8] <- partime[2]
  tests[i,9] <- partime[3]
}

write.csv(tests,"testoutput_2-7-2018_1_gaussian.csv")

## Visualize the Results
## Line Graphs
# n vs. elapsed time
n_plot_3_400 <- tests %>%
  filter(nFolds == 3 & k == 400) %>%
  plot_ly(x = ~n,
          y = ~ser_elapsed_time,
          name = 'Serial Elapsed Time',
          type = 'scatter',
          line = list(shape = "spline")) %>%
  add_trace(y = ~par_elapsed_time,
            name = 'Parallel Elapsed Time')
# k vs. elapsed time
k_plot_3_50 <- tests %>%
  filter(nFolds == 3 & n == 50) %>%
  plot_ly(x = ~k,
          y = ~ser_elapsed_time,
          name = 'Serial Elapsed Time',
          type = 'scatter',
          line = list(shape = "spline")) %>%
  add_trace(y = ~par_elapsed_time,
            name = 'Parallel Elapsed Time')
subplot(n_plot_3_400,
        k_plot_3_50,
        titleX = TRUE,
        titleY = TRUE) %>%
  layout(showlegend = TRUE)

##Reshape the Data
library(reshape2)
tests_elapsed_graph <- tests %>%
  melt(id.vars = c("n", "k", "nFolds")) %>%
  filter(variable == "ser_elapsed_time" | variable == "par_elapsed_time")
## 3D Scatterplot
library(plotly)
tests_elapsed_graph %>%
  filter(nFolds == 3) %>%
  plot_ly(x = ~n,
          y = ~k,
          z = ~value,
          color = ~variable,
          colors = c('#BF382A', '#00CCFF')) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = "n"),
                      yaxis = list(title = "k"),
                      zaxis = list(title = "Elapsed Time")))


