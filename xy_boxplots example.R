devtools::install_github("hammocar/hammond", force = T)
library(hammond)
a <- rnorm(100, 70, 40)
b <- rpois(100, 40)
c <- rbinom(100,2, .2)
data<-data.frame(a = a,
             b = b,
             c = c)


xy_boxplots(data, x = "a", y = "b", color = factor(c)) 
