df_1 = read.table("C:\\Users\\Dell7400\\Documents\\Ale\\Facu\\Multivariado\\datos\\ejemplo_tipo_1.txt",header = TRUE)
df_2 = read.table("C:\\Users\\Dell7400\\Documents\\Ale\\Facu\\Multivariado\\datos\\ejemplo_tipo_2.txt",header = TRUE)

n=nrow(df_2)


df_1_medias = apply(df_1,2,mean)
df_2_medias = apply(df_2,2,mean)

resta = df_1_medias-df_2_medias

s1 = cov(df_1)
s2 = cov(df_2)

s = ((n-1)*s1+(n-1)*s2)/(n+n-2)

# construyo el valor del estadistico To_2

To_2 = ((n*n)/(n+n))*t(resta)%*%solve(s)%*%resta 

p=5
Fo= ((n+n-p-1)/((n+n-2)*p)) * To_2

#busco el fractil de la t student cn p,n1+n2-p-1 grados de libertad

gl1 = p
gl2 = n+n-p-1

Fcritico = qf(0.99, gl1, gl2, lower.tail = T, log.p = F)




