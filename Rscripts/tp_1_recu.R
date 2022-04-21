library(mvShapiroTest)


###### Ejercicio 1 Recu TP_1

df = read.table("C:\\Users\\Dell7400\\Documents\\Ale\\Facu\\Multivariado\\datos\\operarios1.txt")

# a) i)

# Ante todo voy a realizar un shapiro test para ver si los datos son normales

test_shapiro = mvShapiro.Test(as.matrix(df))

test_shapiro$p.value

# El p valor es 0.357, grande, no rechazo Ho, puedo asumir normalidad multivariada


n = nrow(df)
head(df)
var(df)



### Construyo un test aproximado

## voy a usar la distribucion asintotica de -2log(alfa*) que converge en distribucion
## a una chi cuadrado con v grados de libertad

q = (n-1)*var(df)

q11 = (n-1)*var(df[1:3])
q22 = (n-1)*var(df[4:6])

estadistico = det(q)/(det(q11)*det(q22))

# estadistico = 0.7570775

p = ncol(q)
p1 = ncol(q11)
p2 = ncol(q22)


### regla de decision y calculo de p valor

alfa = 0.05
v = (1/2)*(p^2-(p1^2+p2^2))


k_alfa = exp(((-1/n)*qchisq(alfa,v,lower.tail = FALSE)))

ifelse(estadistico <= k_alfa,"Rechazo Ho, la matriz de varianzas no es triangular, es decir, no garantizamos independencia",
       "No Rechazo Ho, podemos garantizar independencia entre las variables del vector aleatorio")



# para el calculo de p valor, la idea es calcular P(estadistico < k_alfa)
# despejando nos queda

# P(chi < - log(estadistico)/1.75), donde chi es una variable aleatoria chi cuadrado con 9 grados de libertad


p_valor = 1-pchisq(log(estadistico)*(-76), df=9)
p_valor

# p_valor  = 0.01200074






### b) i)

##### Test para el vector de medias, usamos el estadistico de Hotelling

muo = as.vector(c(10,25,13,31,28,9))


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


Fo =  (n-p)/(p*(n-1)) * to_2
ifelse(to_2>fcritico,"rechazo","norechazo")


#o bien

t2_critico = ((p*(n-1))/(n-p))*fcritico
ifelse(to_2>t2_critico,"rechazo","norechazo")


# o bien usando p valor

p_valor = 1-pf(to_2*((n-p)/(p*(n-1))), gl1, gl2,lower.tail = T, log.p = F)

# Regla de decision

ifelse(to_2>t2_critico,"rechazo","norechazo")

### Rechaza Ho, suponemos igualdad en el vector de medias vs el valor de mu0 dado


## la direccion "a" que mejor divide a la poblacion es 

a = solve(s)%*%(xbarra-muo)


## ii)

# test para la suma de las medias o 1*medias



unos = c(1,1,1,1,1,1)
fcritico = qf(1-alfa, gl1, gl2, lower.tail = T, log.p = F)

to_2 = n * t(t(unos)%*%xbarra-t(unos)%*%muo) %*% solve(t(unos)%*%s%*%unos) %*% (t(unos)%*%xbarra-t(unos)%*%muo)

p = 1
gl1 = p
gl2 = n-p


t2_critico = ((p*(n-1))/(n-p))*fcritico
ifelse(to_2>t2_critico,"rechazo","norechazo")



Fo =  (n-p)/(p*(n-1)) * to_2
ifelse(to_2>fcritico,"rechazo","norechazo")



# o bien usando p valor

p_valor = 1-pf(to_2*((n-p)/(p*(n-1))), gl1, gl2,lower.tail = T, log.p = F)





#### Intervalos de confianza simultaneos


## Ic hotelling

# creo una funcion 

ic_hotelling = function(n,p,s,xbarravector,a,alfa) {
  
  a_t = as.matrix(a)
  
  error_muestreo = sqrt(t(a_t)%*%s%*%a_t*((p*(n-1))/(n-p))*qf(1-alfa, p, n-p)*(1/n))
  
  at_xbarra_1 = t(a_t)%*%xbarravector
  
  Li = at_xbarra_1-error_muestreo
  Ls = at_xbarra_1+error_muestreo
  
  return(c(Li,Ls))
  }


# cargo los inputs del intervalo
n = nrow(df)
p = ncol(df)
xbarra = sapply(Filter(is.numeric, df), mean)
s = var(df) 
alfa = 0.05

ich1=ic_hotelling(n,p,s,xbarra,c(1,0,0,0,0,0),alfa)
ich2=ic_hotelling(n,p,s,xbarra,c(0,1,0,0,0,0),alfa)
ich3=ic_hotelling(n,p,s,xbarra,c(0,0,1,0,0,0),alfa)
ich4=ic_hotelling(n,p,s,xbarra,c(0,0,0,1,0,0),alfa)
ich5=ic_hotelling(n,p,s,xbarra,c(0,0,0,0,1,0),alfa)
ich6=ic_hotelling(n,p,s,xbarra,c(0,0,0,0,0,1),alfa)



## Ic de bonferroni




ic_bonferroni = function(n,s,xbarravector,a,k,alfa) {
  
  
  a_t = as.matrix(a)
  
  error_muestreo = sqrt(t(a_t)%*%s%*%a_t/n)*qt(1-alfa/(2*k), n-1)
  
  
  at_xbarra_1 = t(a_t)%*%xbarravector
  
  Li = at_xbarra_1-error_muestreo
  Ls = at_xbarra_1+error_muestreo
  
  return(c(Li,Ls))

}

icb1=ic_bonferroni(n,s,xbarra,c(1,0,0,0,0,0),6,alfa)
icb2=ic_bonferroni(n,s,xbarra,c(0,1,0,0,0,0),6,alfa)
icb3=ic_bonferroni(n,s,xbarra,c(0,0,1,0,0,0),6,alfa)
icb4=ic_bonferroni(n,s,xbarra,c(0,0,0,1,0,0),6,alfa)
icb5=ic_bonferroni(n,s,xbarra,c(0,0,0,0,1,0),6,alfa)
icb6=ic_bonferroni(n,s,xbarra,c(0,0,0,0,0,1),6,alfa)


amplitud_bonferroni = c(icb1[2]-icb1[1],icb2[2]-icb2[1],icb3[2]-icb3[1],icb4[2]-icb4[1],icb5[2]-icb5[1],icb6[2]-icb6[1])
amplitud_hotelling = c(ich1[2]-ich1[1],ich2[2]-ich2[1],ich3[2]-ich3[1],ich4[2]-ich4[1],ich5[2]-ich5[1],ich6[2]-ich6[1])


plot(amplitud_hotelling, type="l", col="red")
lines(amplitud_bonferroni , type="l", col="green",title("amplitud de intervalos"))
legend("topright",c("H","B"),fill=c("red","green"))




# se ve en el grafico que la amplitud del intervalo de hotelling es siempre mayor a la de bonferroni
