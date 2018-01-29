source('C:/Liang/Googlebox/Research/Project1/R_pro1/Final/Nindex.R', echo=TRUE)
source('C:/Liang/Googlebox/Research/Project1/R_pro1/Final/sddsim.R', echo=TRUE)
source('C:/Liang/Googlebox/Research/Project1/R_pro1/Final/event_matrix.R', echo=TRUE)
library(DDD)
library(MASS)
library(rgl)
library(stringr)
library(matrixcalc)
library("reshape2")
library('Matrix')
library(plyr) 
library(twitteR)
tree = sddsim(n=2,parsN=c(2,0),age=15,pars=c(0.8,0.3,10) , seed_fun = 29, lambda_allo0 = 0.2, M0=0,K_fix = 1)
print(dim(tree$L))

L = tree$L
time.list = c(sort(c(L[,1],L[which(L[,4]!= -1),4]),decreasing = TRUE),0)
#total number of species
num.species = nrow(L)
trait.table = matrix(0,nrow = length(time.list)-1,ncol = nrow(L)+1)
time.branching = match(L[,1],time.list)
time.end = match(L[,4],time.list)
time.end[is.na(time.end)] = length(time.list)
survival.time = cbind(time.branching,time.end)
timelist = as.data.frame(time.list)
timebranch = as.data.frame(time.branching)
timeend = as.data.frame(time.end)


for(i in 1:num.species){
  
  trait.table[,i+1][time.branching[i]:(time.end[i]-1) ] = 1
}
trait.table = rbind(trait.table,trait.table[dim(trait.table)[1],])
trait.table[,1] = time.list
existing.species.table = trait.table[-1,-1]


write.csv(timelist, file = "C:\\Liang\\Googlebox\\Python\\Project2\\R-tree_sim\\timelist.csv")
write.csv(timebranch, file = "C:\\Liang\\Googlebox\\Python\\Project2\\R-tree_sim\\timebranch.csv")
write.csv(timeend, file = "C:\\Liang\\Googlebox\\Python\\Project2\\R-tree_sim\\timeend.csv")
write.csv(existing.species.table, file = "C:\\Liang\\Googlebox\\Python\\Project2\\R-tree_sim\\traittable.csv")
write.csv(L, file = "C:\\Liang\\Googlebox\\Python\\Project2\\R-tree_sim\\Ltable.csv")


output = list( timelist= time.list, timebranch = time.branching, timeend = time.end, survivaltime = survival.time)
save(output, file = "C:\\Liang\\Googlebox\\Python\\Project2\\R-tree_sim\\treedata.Rdata")

