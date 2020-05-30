
########## Ordenamiento No Métrico Multi-dimensional NMDS - básico vegan ###########

# Datos

## Para este ejemplo se usan los datos de :

# dune del paquete vegan
library(vegan)
df <- data("dune") # abundancias
df_env <- data("dune.env") # ambiental
df <- dune
df_env <- dune.env

# NMDS Básico

# La función metaMDS de vegan ejecuta el NMDS por defecto de vegan, es la función recomendada por el autor de vegan
# otras funciones para ejecutar NMDS son monoMDS (paquete: vegan)  e isoMDS (paquete: MASS)
# Vegan emplea por defecto:

# Estandarización: Wisconsin
# Transformación: raíz cuadrada sqrt()
# Distancia: Bray-Curtis
# NMDS: monoMDS

ord <- metaMDS(df, k = 2, trymax = 100)

# Figura NMDS por defecto de vegan

plot(ord)

# Diagrama de shpard para explorar el ajuste de las distancias observadas vs las distancias calculadas
stressplot(ord)

ordiplot(ord, type = "n")
orditorp(ord, display = "species", col = "red", air = 0.01)
orditorp(ord, display = "sites", cex = 1.25, air = 0.01)

## Diagramas complementados en base a Factores

# Factor:
groups <- df_env$Management

# Polígono:
ordiplot(ord, type = "n")
ordihull(ord, groups = groups,  col = 1:4, lwd = 3)
orditorp(ord, display = "species", col = "red", air = 0.01)
orditorp(ord, display = "sites", cex = 1.25, air = 0.01)

# Elipses ehull
ordiplot(ord, type = "n")
ordiellipse(ord, groups, col = 1:4, kind = "ehull", lwd = 3)

# Elipses 
ordiplot(ord, type = "n")
ordiellipse(ord, groups, col=1:4, draw="polygon")

# Spiders
ordiplot(ord, type = "n")
ordispider(ord, groups, col = 1:4, label = TRUE)

# Combinación de diagramas
ordiplot(ord, type = "n")
ordihull(ord, groups = groups,  col = 1:4, lwd = 3)
ordiellipse(ord, groups, col = 1:4, kind = "ehull", lwd = 3)
ordiellipse(ord, groups, col=1:4, draw="polygon")
ordispider(ord, groups, col = 1:4, label = TRUE)
points(ord, disp="sites", pch=21, col="red", bg="yellow", cex=1.3)

## Diagramas complementados en base a variables numéricas

# Diagrama de burbujas
ordiplot(ord, type = "n")
ordisurf(ord ~ A1, df_env, bubble = 5)

# Diagrama vs variable ambiental
ordiplot(ord, type = "n")
ordisurf(ord ~ A1, df_env,  col = "blue", add = FALSE,
         select = FALSE, method = "GCV.Cp")

# Diagrama de abundancia relativa de una especie

str(df)
ordiplot(ord, type = "n")

fit <- ordisurf(ord ~  Achimill, df, family=quasipoisson)

# Diagrama complementado con un cluster
cl <- hclust(vegdist(df))
plot(ord, type = "p", display="sites")
ordicluster(ord, cl, prune=3, col = cutree(cl, 4))
