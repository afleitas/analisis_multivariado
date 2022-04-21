df= read.table("C:\\Users\\Dell7400\\Documents\\Ale\\Facu\\Multivariado\\datos\\P2-ejCORCHO.txt",header=TRUE)
alfa = 0.005
n = nrow(df)
p = 3 

# Cargo la matriz A

A_t =  matrix(c(1,-1,1,-1,0,1,0,-1,1,0,-1,0),4)
A = t(A_t)

s = var(df)

xbarra = sapply(Filter(is.numeric, df), mean)


x_barra_resta_1 = as.numeric(xbarra[1])-as.numeric(xbarra[2])+as.numeric(xbarra[3])-as.numeric(xbarra[4])
x_barra_resta_2 = as.numeric(xbarra[2])-as.numeric(xbarra[4])
x_barra_resta_3 = as.numeric(xbarra[1])-as.numeric(xbarra[3])



xbarravector = cbind(x_barra_resta_1,x_barra_resta_2,x_barra_resta_3)

S_inv_total = solve(A%*%s%*%A_t) # aca tengo la cuenta de (A*S*A_t)^-1



To_2 = n*xbarravector%*%S_inv_total%*%t(xbarravector)

gl1 = p
gl2 = n-p

# Calculo el p valor transformando el estadistico hotelling a un F snedecor


To_2_c = To_2*((n-p)/(p*(n-1)))

p_valor = pf(To_2_c, gl1, gl2,lower.tail = F, log.p = F)




a_sombrero = S_inv_total%*%t(xbarravector)

s_total = A%*%s%*%A_t


#### Intervalo hotelling

a_t = as.matrix(c(1,0,0))

error_muestreo = sqrt(t(a_t)%*%s_total%*%a_t*((p*(n-1))/(n-p))*qf(1-alfa, gl1, gl2)*(1/n))

at_xbarra_1 = t(a_t)%*%t(xbarravector)


Li = at_xbarra_1-error_muestreo
Ls = at_xbarra_1+error_muestreo

Li
Ls
amplitud = Ls-Li
amplitud



a_t = as.matrix(c(0,1,0))

error_muestreo = sqrt(t(a_t)%*%s_total%*%a_t*((p*(n-1))/(n-p))*qf(1-alfa, gl1, gl2)*(1/n))

at_xbarra_1 = t(a_t)%*%t(xbarravector)


Li = at_xbarra_1-error_muestreo
Ls = at_xbarra_1+error_muestreo

Li
Ls
amplitud = Ls-Li
amplitud



a_t = as.matrix(c(0,0,1))

error_muestreo = sqrt(t(a_t)%*%s_total%*%a_t*((p*(n-1))/(n-p))*qf(1-alfa, gl1, gl2)*(1/n))

at_xbarra_1 = t(a_t)%*%t(xbarravector)


Li = at_xbarra_1-error_muestreo
Ls = at_xbarra_1+error_muestreo

Li
Ls
amplitud = Ls-Li
amplitud





#### Bonferroni

a_t = as.matrix(c(1,0,0))

error_muestreo = sqrt(t(a_t)%*%s_total%*%a_t/n)*qt(1-alfa/6, n-1)


at_xbarra_1 = t(a_t)%*%t(xbarravector)

Li = at_xbarra_1-error_muestreo
Ls = at_xbarra_1+error_muestreo

Li
Ls
amplitud = Ls-Li
amplitud





a_t = as.matrix(c(0,1,0))

error_muestreo = sqrt(t(a_t)%*%s_total%*%a_t/n)*qt(1-alfa/6, n-1)


at_xbarra_1 = t(a_t)%*%t(xbarravector)

Li = at_xbarra_1-error_muestreo
Ls = at_xbarra_1+error_muestreo

Li
Ls
amplitud = Ls-Li
amplitud




a_t = as.matrix(c(0,0,1))

error_muestreo = sqrt(t(a_t)%*%s_total%*%a_t/n)*qt(1-alfa/6, n-1)


at_xbarra_1 = t(a_t)%*%t(xbarravector)

Li = at_xbarra_1-error_muestreo
Ls = at_xbarra_1+error_muestreo

Li
Ls
amplitud = Ls-Li
amplitud


### test normalidad

library(MVN)
test = mvn(df, mvnTest = "hz", multivariateOutlierMethod = "quan")
pvalue = test$multivariateNormality[3]
pvalue
ifelse(pvalue<alfa,"Rechazo, no hay normalidad","No rechazo, hay normalidad multivariada")









