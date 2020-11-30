
########### Análisis de componentes Principales ############


### 1. Paquetes ###########

library(vegan)
library(tidyverse)

## 2. Datos #######

# Para este ejemplo se usan los datos de:
data(dune)
data(dune.env)


df <- dune
df_env <- dune.env

### 1. Calculamos el PCA ########

ord <- rda(df, scale = TRUE)

summary(ord)

### Interpretación:

## (1) Inertia: Es el término general para la "variación" en los datos. En PCA es ala suma de varaianzas de las variablesl
## (2) Unconstrained: Indica que el análisis no limitado por otras variables, contiene un set de variables autoexplicativas.
## (3) Eigenvalues (Autovalores): Valores de la importancia de los ejes del PCA (PC1, PC2, ....)
## (4) Proportion Explained: o proporciones de variación contadas para cada eje, dividiendo cada autovalor por la "incia total".
## (5) Scaling: No confundir con el argumento "scale", empleado para la estandarización de variables.
      # Hace referencia a la forma en que los resultados son mostrados:
      # Escalamiento 1: Biplot de distancia: Los autovectores son escalados a la unidad. Distancias entre objetos en el biplot son aproximaciones de sus distancias euclidianas. Los ángulos no reflejan las correlaciones
      # Escalamiento 2: Biplot de correlación: Cada autovector es escalado a la raíz cuadrada de su autovalor. Las distancias no son aproximaciones de su distancia euclidiana. Los ángulos reflejan sus correlaciones.
      # Escalamiento 3: Escalamiento simétrico: consiste en ambos escalamientos, permite la presentación simultánea de los sitios y puntajes.

## (6) Species scores: Coordenadas de las flechas de las variables. Por raones históricas, las variables en vegan reciben el nombre de "species", sin importar lo que representen.
## (7) Site scores: Coordenadas de los sitios en el diagrama de ordenación. Los objetos son siempre llamados "Sitios" en los resultados de vegan.

### 2. Ploteando ########

# Funciones previas

# Introducimos las funciones ord_labesls y scale_arrow , extraídas del paquete ggordiplots

################# Función ord_labels ############################

ord_labels <-
  function(ord){
    ev <- vegan::eigenvals(ord)
    if (!is.na(ev)[1]) {
      tol <- -(1e-07)*ev[1]
      ord.labels <- rep("", length(ev))
      if ((any(is.na(ev))) | (any(ev < tol))) {
        for ( i in 1:length(ev)) {
          ord.labels[i] <- paste("DIM", i, sep = "")
        }
      }
      else {
        ev.pc <- round(100*(ev/sum(ev)), 2)
        axis.names <- names(ev)
        if (is.null(axis.names)) {
          for ( i in 1:length(ev.pc)) {
            ord.labels[i] <- paste("DIM", i, " ", sprintf(ev.pc[i], fmt = '%#.1f'), "%", sep="")
          }
        } else {
          for (i in 1:length(ev.pc)){
            ord.labels[i] <- paste(axis.names[i], " ", ev.pc[i],"%", sep="")
          }
        }
      }
    } else {
      ord.labels <- colnames(vegan::scores(ord))
    }
    
    return(ord.labels)
  }


################## Función scale_arrow ###########################

scale_arrow <- function(arrows, data, at = c(0, 0), fill = 0.75) {
  u <- c(range(data[,1], range(data[,2])))
  u <- u - rep(at, each = 2)
  r <- c(range(arrows[, 1], na.rm = TRUE), range(arrows[, 2], na.rm = TRUE))
  rev <- sign(diff(u))[-2]
  
  if (rev[1] < 0) {
    u[1:2] <- u[2:1]
  }
  if (rev[2] < 0) {
    u[3:4] <- u[4:3]
  }
  u <- u/r
  u <- u[is.finite(u) & u > 0]
  invisible(fill * min(u))
}

choices = c(1,2)

# Extraemos los scores

df_pca <- as.data.frame(scores(ord, display = "sites", scaling = 1, choices = c(1,2)))
df_pca

# Elegimos etiquetas de sitios y ejes
axis.labels <- ord_labels(ord)[choices]
xlab <- axis.labels[1]
ylab <- axis.labels[2]
colnames(df_pca) <- c("x", "y")
df_pca$site <- rownames(df_pca)

# Diagrama Básico
(
  ggplot(data = df_pca, aes(x = x, y = y)) +
    geom_point(size = 3)+
    xlab(xlab) +
    ylab(ylab) +
    geom_text(mapping = aes(label = site, x = x, y = y),
              size = 2, vjust=-1) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()
)
