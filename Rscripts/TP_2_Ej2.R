df = read.table("C:\\Users\\Dell7400\\Documents\\Ale\\Facu\\Multivariado\\datos\\centollas.txt")
head(df)


## para el test de igualdad de matrices de var cov voy a usar el test de maximaverosimilitud pag 53 del apunte de inferencia normal II, test asintotico

## como paso 1 separo las poblaciones


#### Distribuciones asociadas
# qi son wishart (sigma,p.n_i-1)
# u es wihart (sigma,p,n-k)
# k es la cantidad de q_i que hay, osea la cantidad de poblaciones


df_1 = df[df$V1==1,][,-1]
df_2 = df[df$V1==2,][,-1]
df_3 = df[df$V1==3,][,-1]

n1=nrow(df_1)
n2=nrow(df_2)
n3=nrow(df_3)
n=n1+n2+n3
k = 3
p = ncol(df[,-1])

v = p*k+(p*(p+1))*(k/2)
v1 = p*k+(p*(p+1))*(1/2)



q1 = (n1-1)*var(df_1)
q2 = (n2-1)*var(df_2)
q3 = (n3-1)*var(df_3)
u = q1+q2+q3


alfa_estrella = ((det(q1/n1)^{n1/2})*(det(q2/n2)^{n2/2})*(det(q3/n3)^{n3/2}))/((det(u/n))^{n/2})

# -2*log(alfa_estrella) converge en distribucion a una chi(v-v1)


gl = v-v1
  
pvalue = pchisq(-2*log(alfa_estrella),gl,lower.tail = FALSE)

# pvalue 0.7088, no rechazo Ho, podemos asumir igualdad de matrices de varianzas y covarianzas



#b) Estadistico para igualdad de medias 


# como nos pide nivel exacto, voy a tener que usar el criterio de wilks, nos basamos en el test presentado en pagina 55 del apunte de inferencia normal II



# como k=3, entonces r=k-1

r = k-1


# la distribucion bajo Ho es wilks(n-1,p,k-1)



media_1 = sapply(df_1,mean)
media_2 = sapply(df_2,mean)
media_3 = sapply(df_3,mean)


gran_media = (1/n)*(n1*media_1+n2*media_2+n3*media_3)

h = n1*(media_1 - gran_media)%*%t(media_1 - gran_media)+n2*(media_2 - gran_media)%*%t(media_2 - gran_media)+n3*(media_3 - gran_media)%*%t(media_3 - gran_media)



wilks = det(u)/det(u+h)


estadistico = (1-wilks^{1/2})/(wilks^{1/2}) * (n-1-p-1) * (1/p)


v1 = 2*r
v2 = n-1-r-1


# calculo el p valor

pf(estadistico,v1,v2,lower.tail = FALSE)


# el p valor da casi 0, entonces rechazo Ho, no puedo suponer igualdad de medias

# esto me da a sospechar, siendo que el p valor va a ser muy chico, que deberia encontrar un buen discriminante lineal



## c)


# recordar que u/n estima a sigma_w, h/m estima a sigma_b, entonces la coordenada discriminante alfa_j*x donde alfa_j es el autovector j de la matriz sigma_w^{-1}*sigma_b



# PREGUNTA, ES RAZONABLE HACER UN PLOT DE LAS PRMERAS 3?




#### 
#### de esta manera me estan dando valores con numeros imaginarios
####

# encuentro las direcciones canonicas

eigen_ = eigen(solve(u/n)%*%(h/n))


alfa_1 = eigen_$vector[,1]*(-1)
alfa_2 = eigen_$vector[,2]*(-1)
alfa_3 = eigen_$vector[,3]*(-1)
alfa_4 = eigen_$vector[,4]*(-1)
alfa_5 = eigen_$vector[,5]*(-1)
####################################################################






#### FORMA 2 uso decomposition de cholesky

c = chol(u)
c_t = t(c)

b = t(solve(c))%*%(h/n)%*%solve(c)

eigen_b = eigen(b)

# alfa_j = c^{-1} * b_j, con b_j autovector de b


alfa_1 = solve(c)%*%eigen_b$vector[,1]
alfa_2 = solve(c)%*%eigen_b$vector[,2]*(-1)
alfa_3 = solve(c)%*%eigen_b$vector[,3] 
alfa_4 = solve(c)%*%eigen_b$vector[,4]
alfa_5 = solve(c)%*%eigen_b$vector[,5]




## d ) 

# sea z_j variable canonica o discriminante


A = as.matrix(cbind(alfa_1,alfa_2))
  
xproy =  t(A)%*%t(df[,-1])
x1proy = t(A)%*%t(df_1)
x2proy = t(A)%*%t(df_2)
x3proy = t(A)%*%t(df_3)



zraya1 = t(A)%*%media_1
zraya2 = t(A)%*%media_2
zraya3 = t(A)%*%media_3


plot(xproy[1,],xproy[2,],pch=20, ylab="")
points(x1proy[1,],x1proy[2,], col="blue", pch=20, ylab="")
points(x2proy[1,],x2proy[2,], col="yellow", pch=20, ylab="")
points(x3proy[1,],x3proy[2,], col="green", pch=20, ylab="")
points(zraya1[1,],zraya1[2,], pch=17, lwd=4, col="blue")
points(zraya2[1,],zraya2[2,], pch=17, lwd=4, col="yellow")
points(zraya3[1,],zraya2[2,], pch=17, lwd=4, col="green")


# que observa?  tiene sentido proyectar en ambas coordenadas, las proyecciones discriminan en ambos ejes

# PREGUNTA
#Que interpretacion le daria a la primer
#variable canonica, en relacion a las medidas del caparazon?

# ver los pesos del alfa
  

##e)

J<-zraya1+zraya3
C<-zraya1-zraya3
M<--C[1]/C[2]
oo<-t(C)%*%J/(2*C[2])
#abline(oo,M,col="violet")



segments(x0 = 0.27,
         x1 = 0.1565,
         y0 = -.3,
         y1 = 0.078,
         lwd = 2,
         col = "violet") 


J<-zraya1+zraya2
C<-zraya1-zraya2
M<--C[1]/C[2]
oo<-t(C)%*%J/(2*C[2])
#abline(oo,M,col="black")


segments(x0 = 0.1565,
         x1 = -0.4,
         y0 = 0.078,
         y1 = -0.20547,
         lwd = 2,
         col = "black") 


J<-zraya2+zraya3
C<-zraya2-zraya3
M<--C[1]/C[2]
oo<-t(C)%*%J/(2*C[2])
#abline(oo,M,col="red")

segments(x0 = 0.1565,
         x1 = 1,
         y0 = 0.078,
         y1 = 7,
         lwd = 2,
         col = "red") 





## f)

# armo una funcion distancia euclidea


euc_dist = function(x1, x2) sqrt(sum((x1 - x2) ^ 2))


# funcion de clasificacion


clasificador = function(x,A)
  
# x vector de variables
# A matriz de proyeccion


  {
  
  proy = t(A)%*%x
  
  
  zraya1 = t(A)%*%media_1
  zraya2 = t(A)%*%media_2
  zraya3 = t(A)%*%media_3
  
  distancia = c(euc_dist(proy,zraya1),euc_dist(proy,zraya2),euc_dist(proy,zraya3))
  
  regla = ifelse(euc_dist(proy,zraya1) == min(distancia), "grupo1",
                  ifelse(euc_dist(proy,zraya2) == min(distancia),"grupo2",
                         ifelse(euc_dist(proy,zraya3) == min(distancia),"grupo3")
                         )
  )
  
  return(c(regla,distancia))
}





## g)

clasificador(c(11.2, 10.7, 27.6, 31.8, 11),A)


# grupo 2 con distancia 0,367




## h, calcular las componentes principales, primero centro los valores de Xs

library(dplyr)
library(xml2)
library(rvest)





# para df_1

s1 = var(df_1)
v1 = eigen(s1)$vectors[,1]*(-1)
v2 = eigen(s1)$vectors[,2]*(-1)
v3 = eigen(s1)$vectors[,3]*(-1)
v4 = eigen(s1)$vectors[,4]*(-1)
v5 = eigen(s1)$vectors[,5]*(-1)






pca1 = as.vector(v1%*%t(scale(df_1, center = TRUE, scale = FALSE)))
pca2 = as.vector(v2%*%t(scale(df_1, center = TRUE, scale = FALSE)))
pca3 = as.vector(v3%*%t(scale(df_1, center = TRUE, scale = FALSE)))
pca4 = as.vector(v4%*%t(scale(df_1, center = TRUE, scale = FALSE)))
pca5 = as.vector(v5%*%t(scale(df_1, center = TRUE, scale = FALSE)))






# para df_2

s2 = var(df_2)
v1 = eigen(s2)$vectors[,1]*(-1)
v2 = eigen(s2)$vectors[,2]*(-1)
v3 = eigen(s2)$vectors[,3]*(-1)
v4 = eigen(s2)$vectors[,4]*(-1)
v5 = eigen(s2)$vectors[,5]*(-1)






pca1 = as.vector(v1%*%t(scale(df_2, center = TRUE, scale = FALSE)))
pca2 = as.vector(v2%*%t(scale(df_2, center = TRUE, scale = FALSE)))
pca3 = as.vector(v3%*%t(scale(df_2, center = TRUE, scale = FALSE)))
pca4 = as.vector(v4%*%t(scale(df_2, center = TRUE, scale = FALSE)))
pca5 = as.vector(v5%*%t(scale(df_2, center = TRUE, scale = FALSE)))





# para df_3

s3 = var(df_3)
v1 = eigen(s3)$vectors[,1]*(-1)
v2 = eigen(s3)$vectors[,2]*(-1)
v3 = eigen(s3)$vectors[,3]*(-1)
v4 = eigen(s3)$vectors[,4]*(-1)
v5 = eigen(s3)$vectors[,5]*(-1)






pca1 = as.vector(v1%*%t(scale(df_3, center = TRUE, scale = FALSE)))
pca2 = as.vector(v2%*%t(scale(df_3, center = TRUE, scale = FALSE)))
pca3 = as.vector(v3%*%t(scale(df_3, center = TRUE, scale = FALSE)))
pca4 = as.vector(v4%*%t(scale(df_3, center = TRUE, scale = FALSE)))
pca5 = as.vector(v5%*%t(scale(df_3, center = TRUE, scale = FALSE)))

