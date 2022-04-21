
# Ejercicio 2 

#a

df= read.table("C:\\Users\\Dell7400\\Documents\\Ale\\Facu\\Multivariado\\datos\\P2-ejPESO-ALT.txt",header=TRUE)
df$Altura = df$Altura/1000


### Hipotesis
# Ho: (mu_peso,mu_altura) = (63,1.60) vs H1: (mu_peso,mu_altura) != (63,1.60)


n = nrow(df)
p = ncol(df)

xbarra = c(mean(df$Peso),mean(df$Altura))
muo = c(63,1.60)
s = var(df) #cov cuidado porq creo q divide por n
alfa = 0.05

# To_2 = n (xbarra - muo)t (S-1) (xbarra-mu)

to_2 = n * t(xbarra-muo) %*% solve(s) %*% (xbarra-muo)

gl1 = p
gl2 = n-p

fcritico = qf(1-alfa, gl1, gl2, lower.tail = T, log.p = F)

t2_critico = ((p*(n-1))/(n-p))*fcritico


# Regla de decision

ifelse(to_2>t2_critico,"rechazo","norechazo")


#b

# Elipsoide de confianza

confianza = 0.95

n*(n-p)/(p*(n-1))
solve(s)

# n * (n-p)/(p(n-1))*(xbarra-mu)t S-1 (xbarra-mu) <= fp,n-p(alfa)


#c

# test shapiro-wilks



library(MVN)

test = mvn(df, mvnTest = "hz", multivariateOutlierMethod = "quan")
pvalue = test$multivariateNormality[3]

ifelse(pvalue<alfa,"rechazo, no hay normalidad","no rechazo, hay normalidad multivariada")




# Ejercicio 7


df= read.table("C:\\Users\\Dell7400\\Documents\\Ale\\Facu\\Multivariado\\datos\\P2-ejCORCHO.txt",header=TRUE)
alfa = 0.01

head(df)

#a)  test de normalidad


library(MVN)

test = mvn(df, mvnTest = "hz", multivariateOutlierMethod = "quan")
pvalue = test$multivariateNormality[3]
pvalue
ifelse(pvalue<alfa,"rechazo, no hay normalidad","no rechazo, hay normalidad multivariada")


#b) 


#i
muo = as.vector(c(45,42,45,42))


n = nrow(df)
p = ncol(df)


xbarra = sapply(Filter(is.numeric, df), mean)


s = var(df) #cov cuidado porq creo q divide por n
alfa = 0.05

# To_2 = n (xbarra - muo)t (S-1) (xbarra-mu)

to_2 = n * t(xbarra-muo) %*% solve(s) %*% (xbarra-muo)

gl1 = p
gl2 = n-p

fcritico = qf(1-alfa, gl1, gl2, lower.tail = T, log.p = F)
t2_critico = ((p*(n-1))/(n-p))*fcritico
ifelse(to_2>t2_critico,"rechazo","norechazo")


# o bien usando p valor

p_valor = 1-pf(to_2*((n-p)/(p*(n-1))), gl1, gl2,lower.tail = T, log.p = F)

# Regla de decision

ifelse(to_2>t2_critico,"rechazo","norechazo")

### No se rechaza Ho, suponemos igualdad conjunta en las medias 


#ii)

# me olvide de hacerlo pero la idea es armar una matriz A mas grande para comprar





#iii)  aca la idea es usar lo visto en el punto 3, 


head(df)

# cargo la matriz A

A_t =  matrix(c(1,0,-1,0,0,1,0,-1),4,2)
A = t(A)
s

xbarra_n_s = as.numeric(xbarra[1])- as.numeric(xbarra[3])
xbarra_e_o = as.numeric(xbarra[2])- as.numeric(xbarra[4])

xbarravector = cbind(xbarra_n_s,xbarra_e_o)

S_inv_total = solve(A%*%s%*%A_t) # aca tengo la cuenta de (A*S*A_t)^-1

To_2 = n*xbarravector%*%S_inv_total%*%t(xbarravector)

gl1 = p
gl2 = n-p

# lo hago distinto ahora para cambiar, dejo el F critico y paso el T0_2 a una F haciendo la cuenta inversa

fcritico = qf(alfa, p, n-p, ncp=0, lower.tail = FALSE, log.p = FALSE)

t2_critico = ((p*(n-1))/(n-p))*fcritico #---> asi habia hecho antes




To_2_c = To_2*((n-p)/(p*(n-1)))
p_valor = 1-pf(To_2_c, gl1, gl2,lower.tail = F, log.p = F)



# el p valor es casi 0, entonces rechazo Ho, no existe evidencia para afirmar las igualdades simultaneas


#iv Intervalos de confianza simultaneos


# uso la funcion que armo Dani que se la ve simple


### Para Hotelling = Sheffe se tiene

hotelling.IC1m <- function(X, at, alpha){
  # Dada una muestra, un vector at y un nivel alpha, 
  # devuelve el IC para at mu de nivel 1-alpha; adem?s de Xraya y Sx
  
  # Datos y dimensiones
  X <- as.matrix(X)
  n <- dim(X)[1]
  p <- dim(X)[2]
  at <- t(as.matrix(at))
  
  # Valores muestrales
  xraya <- colMeans(X)
  Sx <- cov(X)
  
  # L?mites del IC
  T2 <- ((p*(n-1))/(n-p)) * qf(alpha, p, n-p, ncp=0, lower.tail = FALSE, log.p = FALSE)
  linf <- at%%xraya - sqrt(at%%Sx%*%t(at)) * sqrt(T2) / sqrt(n)
  lsup <- at%%xraya + sqrt(at%%Sx%*%t(at)) * sqrt(T2) / sqrt(n)
  
  IC <- cbind(linf, lsup)
  # Salidas
  return(list(media=xraya, var=Sx, IC_atmu = IC))
}









