df_diesel = read.table("C:\\Users\\Dell7400\\Documents\\Ale\\Facu\\Multivariado\\datos\\diesel.txt")
df_nafta = read.table("C:\\Users\\Dell7400\\Documents\\Ale\\Facu\\Multivariado\\datos\\nafta.txt")


# me fijo si existe distribucion normal multivariada en cada poblacion
#test para normalidad multivariada

library(MVN)
library(dplyr)



mvn(df_diesel, mvnTest = "hz", multivariateOutlierMethod = "quan")




# calculo vectores de medias


x_barra_diesel = apply(df_diesel, 2, mean)
x_barra_nafta = apply(df_nafta, 2, mean)
resta_barraa = x_barra_diesel-x_barra_nafta


df_diesel-df_nafta

cov(df_diesel,df_nafta)








  
##### para una muestra el ejemplo q tengo es

x = cbind(peso,altura)
rownames(x)=1:length(peso)
plot(peso,altura, col="lightblue", pch=19)
text(peso, altura, labels=rownames(x), cex=0.9, font=2)
points(peso[39],altura[39],col="red",pch=20)
y = x[-39,]
oy = apply(y,2,mean)


wsigma = cov(y)
mu0 = cbind(63.64,1615.38)
n = dim(y)[1]
p = dim(y)[2]
T0 = n*(oy-mu0)%*% solve(wsigma)%*%t(oy-mu0)
F0 = T0*(n-p)/(p*(n-1))
falpha = qf(0.001, p, n-p, ncp=0, lower.tail = FALSE, log.p = FALSE)
pvalor = pf(F0, p, n-p, ncp=0, lower.tail = FALSE, log.p = FALSE)

