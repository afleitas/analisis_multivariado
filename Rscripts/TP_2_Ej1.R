

#####
##### EJERCICIO 1 
#####


df = read.table("C:\\Users\\Dell7400\\Documents\\Ale\\Facu\\Multivariado\\datos\\hombres.txt",header=TRUE)

head(df)



s = var(df)


eigen = eigen(s)

eigen$vectors

eigen$vectors[,1][1]
eigen$vectors[,2][2]
eigen$vectors[,3][3]
eigen$vectors[,4][4]
eigen$vectors[,5][5]
eigen$vectors[,6][6]

# segun la regla pedida, elgiendo el vector i esimo de modo que su coordenada i sea positiva

autovectores= cbind(eigen$vectors[,1]*(-1),eigen$vectors[,2]*(-1),eigen$vectors[,3]*(-1),eigen$vectors[,4],eigen$vectors[,5],eigen$vectors[,6])



# donde las columnas de la matriz autovectores son las direcciones principales



# las coordenadas de la primer compoenten principal son todas positivas (o negativas dependiendo la direccion) porque la matriz 
# S tiene todas sus coordenadas positivas


## para mirar cuantas variables hacen falta para explicar un determinado % de varianzas construyo un vector con las varianzas
## usando la traza de la matriz como denominador


traza = sum(eigen$values)
varianzas = (1/traza)*c(eigen$values[1],eigen$values[2],eigen$values[3],eigen$values[4],eigen$values[5],eigen$values[6])


# con la funcion cumsum puedo ver la varianza acumulada explicada, vemos por ejemplo que con las 3 primeras componentes principales se alcanza
# casi el 80% de la variabilidad total

cumsum(varianzas)


#0.4217673 0.6278385 0.7824267 0.8817436 0.9587639 1.00




### PREGUNTA CUANTAS PARECEN SUFICIENTES? HAY UNA REGLA O ES A OJO?

# no puedo discriminar entre los subespacios





#signo = c()

#for (i in 1:length(eigen$vectors[,1]))
#{

#  print(ifelse(eigen$vectors[,i][i]>0,eigen$vectors[,i],eigen$vectors[,i]*(-1)))
#}


ifelse(eigen$vectors[,1][1]>0,eigen$vectors[,1],eigen$vectors[,1]*(-1))



#b )  test de normalidad

library(mvShapiroTest)

test_shapiro = mvShapiro.Test(as.matrix(df))

test_shapiro$p.value

# El p valor es 0.7592, grande, no rechazo Ho, puedo asumir normalidad multivariada




# c) test de esfericidad



# Ho: matriz_sigma = sigma*1 (1 es la matriz identidad de igual dimension que la matriz sigma)
# H1: matriz_sigma != sigma*1

# es un test asintotico, me baso del test del apunte inferencia_normal_II pagina 49.

alfa = 0.05
p = ncol(df)
n = nrow(df)
q = (n-1)*s
gl = p*(p+1)*0.5-1
k_alfa = exp((-1/(n*p))*qchisq(alfa,gl,lower.tail = FALSE))



traza = sum(eigen(q)$values)


estadistico = ((det(q)^{1/p})/((1/p)*traza))


ifelse(estadistico<=k_alfa,"rechazo Ho, no podemos asegurar esfericidad","No rechazo Ho, no podemos rechazar esfericidad")


# como rechazo Ho, entonces sigo con el item 2, me baso en el test presentado en la slide de PCA pagina 40

# tenemos un test cociente de maximaverosimilitud

# Ho: lambda2=lambda3=lambda4=lambda5

r = 1
h = 4
n = nrow(df)
gl = (h*(h+1))/2-1

lambda1_est = eigen$values[1]
lambda2_est = eigen$values[2]
lambda3_est = eigen$values[3]
lambda4_est = eigen$values[4]
lambda5_est = eigen$values[5]


estadistico = (lambda2_est*lambda3_est*lambda4_est*lambda5_est)/(((1/h)*(lambda2_est+lambda3_est+lambda4_est+lambda5_est))^{h})
  
nlog_estadistico = -n*log(estadistico)


pvalue = pchisq(nlog_estadistico,gl,lower.tail = FALSE)


# como el pvalue = 0.001245 < 0.05 , rechazo Ho




# agregar lambda1=lambda2, lambda5=lambda6

# parece q me puedo quedar con lambda 1 o lambda 1 a 5 



# como rechazo la Ho, la idea ahora es hacer test de dos en dos para entender que pares de lambda son iguales


# Ho: lambda1=lambda2

r = 0
h = 2
n = nrow(df)
gl = (h*(h+1))/2-1

estadistico = (lambda1_est*lambda2_est)/(((1/h)*(lambda1_est+lambda2_est))^{h})

nlog_estadistico = -n*log(estadistico)


pvalue = pchisq(nlog_estadistico,gl,lower.tail = FALSE)



# pvalue= 0.002564963 no rechazo Ho



# Ho: lambda2=lambda3

r = 1
h = 2
n = nrow(df)
gl = (h*(h+1))/2-1

estadistico = (lambda2_est*lambda3_est)/(((1/h)*(lambda2_est+lambda3_est))^{h})

nlog_estadistico = -n*log(estadistico)


pvalue = pchisq(nlog_estadistico,gl,lower.tail = FALSE)

# pvalue= 0.3761 no rechazo Ho



# Ho: lambda3=lambda4

r = 2
h = 2
n = nrow(df)
gl = (h*(h+1))/2-1

estadistico = (lambda3_est*lambda4_est)/(((1/h)*(lambda3_est+lambda4_est))^{h})

nlog_estadistico = -n*log(estadistico)


pvalue = pchisq(nlog_estadistico,gl,lower.tail = FALSE)

# pvalue= 0.09966478 no rechazo Ho



# Ho: lambda4=lambda5

r = 3
h = 2
n = nrow(df)
gl = (h*(h+1))/2-1

estadistico = (lambda4_est*lambda5_est)/(((1/h)*(lambda4_est+lambda5_est))^{h})

nlog_estadistico = -n*log(estadistico)


pvalue = pchisq(nlog_estadistico,gl,lower.tail = FALSE)

# pvalue =  0.4650727 no rechazo Ho,




####
#### Conclusion ejercicio c
####


# rechazo igualdad conjunta de los lambda
# no rechazo igualdad de a pares

# PREGUNTAR ESTO



## d) depende lo que entiendo arriba 




## e) Realzio el test basandonos en la regla de decision explicada en la pagina 46 del apunte de PCA 



# Hipotesis

# Ho: tita_po <= 0 vs H1: tita_po > 0 

# me creo esta dos funciones primero que las voy a usar en tita_test, sigma_2_test


# bajo Ho tenemos que 

po = 0.9
q = 5
autoval = eigen$values



sumfun<-function(x,start,end){
  return(sum(x[start:end]))
}


na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}






tita_est = (1-po)*(sumfun(autoval,1,q)-na.zero(po*(sumfun(autoval,q+1,length(autovalores)))))


sigma_2_est = (2*(1-po)**2)*(sumfun(autoval**2,1,q)+2*(po**2)*na.zero((sumfun(autoval**2,q+1,length(autovalores)))))



estadistico = sqrt(n)*(abs(tita_est)/sqrt(sigma_2_est))



# luego como el test es a dos colas, comparo con el percentil de una normal usando todo el nivel alfa



z_alfa = qnorm(1-alfa)



# como estadistica > z_alfa, entonces rechazo Ho, tenemos evidencia estadistica para afirmar que las primeras cinco componentes principales explican mas del
# %90 de la variabilidad




#####
##### Algunas tecnicas de visualizacion
#####

head(df)
ggcorr(df)

