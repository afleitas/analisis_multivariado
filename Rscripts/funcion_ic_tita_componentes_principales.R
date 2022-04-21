A = matrix(c(1,2,4,1,2,5,8.5,4,4,8.5,17.25,5,1,4,5,6),4,4)

data("USArrests")


df = read.table("C:\\Users\\Dell7400\\Documents\\Ale\\Facu\\Tecnicas_reduccion\\datos\\ratones.txt")


df = na.omit(df[,-1])

df$V3 = as.numeric(df$V3)
df$V4 = as.numeric(df$V4)
df$V5 = as.numeric(df$V5)
df$V6 = as.numeric(df$V6)
df$V7 = as.numeric(df$V7)
df$V8 = as.numeric(df$V8)
df$V9 = as.numeric(df$V9)


A = var(na.omit(df[,-1]))

autovalores = eigen(A)$values


sumfun<-function(x,start,end){
  return(sum(x[start:end]))
  }


na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}




ic_tita = function(po,autoval,n,q,alfa){

    
  
tita = (1-po)*(sumfun(autoval,1,q)-na.zero(po*(sumfun(autoval,q+1,length(autovalores)))))


sigma_2 = (2*(1-po)**2)*(sumfun(autoval**2,1,q)+2*(po**2)*na.zero((sumfun(autoval**2,q+1,length(autovalores)))))


li = tita-(sqrt(sigma_2)/sqrt(n))*qnorm(alfa, mean = 0, sd = 1, lower.tail = FALSE)
ls = tita+(sqrt(sigma_2)/sqrt(n))*qnorm(alfa, mean = 0, sd = 1, lower.tail = FALSE)

ic = c(li,ls)

return(ic) 
  
}



a = ic_tita(0.8,autovalores,10,1,0.05)


