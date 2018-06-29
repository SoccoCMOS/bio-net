###### Test distributions of data ####

setwd("../data")
library(vcd)
library(reshape2)
library(ggplot2)

data<-read.csv2("adjacence_esp_count.csv",row.names = 1, header=TRUE)

data$ids=rownames(data)

longdata=melt(data,value.name="count",id.vars = "ids")
v=longdata$count
hist(v)
  
gf=goodfit(v$count,type="poisson",method="MinChisq")
summary(gf)
