##################### CURVAS DE RAREFACIÓN Y EXTRAPOLACIÓN BASADAS EN NÚMEROS DE HILL #################

# FUENTES: 
# 1. T. C. Hsieh  K. H. Ma  Anne Chao (2016). iNEXT: an R package for rarefaction and extrapolation of species diversity (Hill numbers) 
#  https://doi.org/10.1111/2041-210X.12613

# 2. Appendix iNEXT: An R package for rarefaction and extrapolation of species diversity (Hill numbers)
# https://besjournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2F2041-210X.12613&file=mee312613-sup-0001-AppendixS1.pdf

## 1. Paquetes ###############

# Instalar sólo si no se ha instalado previamente:

## install.packages("iNEXT")

# Importar paquetes

library(iNEXT)
library(tidyverse)

## 2. Datos ##############

### 2.1 Datos de abundancia de individuos ############

# Datos de cada ensamblaje/sitio con abundancias en una muestra empírica de n individuos.
# Cuando hay N ensamblajes, los datos consisten de una matriz S x N o N listas de abundancias.

## Datos de ejemplo: 

# Tipo lista
data(spider)

# Tipo matriz (data frame) # Nota: Los nombres de las especies deben ir como nombres de filas o no ir.
data(bird)


### 2.2 Datos de incidencia ##############



## 3. Ejemplo de Rarefacción y Extrapolación (Abundancia) #########

data(spider)
str(spider)

# Cálculo

iNEXT(spider, q = 0, datatype = "abundance")

# Interpretación 

## 0. Sitio de referencia

# Si se analizan varios sitios a la vez, el sitio de referencia es aquel con la menor abundancia total de todos.
# Y es aquel sobre el que se realizan las comparaciones rarificadas.

## 1. $class: iNEXT , es el tipo de objeto.

## 2. $DataInfo: información básica del modelo para cada sitio (site)
# (n): tamaño de la muestra de referencia
# (S.obs): Riqueza observada de especies
# (SC): Estimado de cobertura de muestra
# (f1-f10): 10 primeros conteos de frecuencias. 
# f1: Especies representadadas por 1 sólo individuo (singletons)
# f2: Especies representadas por 2 individuos (doubletons)
# fk: número de especies representadas por k individuos

# Nota: el argumento knots está en 40 por defecto, y separa las abundancias en 40 puntos o knots (igualmente espaciados)
# en una longitud de 1 a 2(abundancia), si una muestra tiene 168 individuos el sistema toma 40 puntos (1, 10, 19,.., 336)
# Los estimados de diversidad se realizan para cada knot (nudo).

## 3. $ iNextEst: Estimados de diversidad con muestras rarificadas y extrapoladas

# (m): Tamaño de muestra de cada uno los 40 knots (por defecto)
# (method): Puede ser interpolado (< obs), observado (= obs), extrapolado (> obs)
# (order): Orden de diversidad según Números de Hill q = 0,1,2
# (qD): Estimado de la diversidad de orden q
# (qD.LCL): Límite inferior del intervalo de confianza del estimado de diversidad de orden q
# (qD.UCL): Límite superior del intervalo de confianza del estimado de diversidad de orden q
# (SC): Estimado de la Cobertura de muestra (Sample coverage) 
# (SC.LCL): Límite inferior del intervalo de confianza del estimado de Cobertura de muestra 
# (SC.UCL): Límite inferior del intervalo de confianza del estimado de Cobertura de muestra 

## 4. $AsyEst Estimados de diversidad Asintótica y estadísticos relacionados
# Lista la diversidad observada, los estimados asintóticos, y los intervalos de confianza al 95%

# (Site): Sitio 
# (Diversity): Los estimados son calculados mediante las funciones ChaoRichness(), ChaoShannon() y ChaoSimpson()
# (Observed): Valor observado de diversidad
# (Estimator): Estimador calculado
# (s.e.): Error estándar
# (LCL): Límite inferior del intervalo de confianza
# (UCL): Límite superior del intervalo de confianza


## Nota:

# El argumento endpoint especifica un tamaño máximo para el cálculo de R/E.
# Para la riqueza este método es confiable hasta el doble del tamaño del sitio de referencia.

# El argumento knots puede ser grande, sin embargo el tiempo de cálculo puede ser demasiado largo.

# El argumento size recibe valores de abundancia en los que se quiere analizar

m <- c(1, 5, 20 , 50, 100, 200, 400)
iNEXT(spider, q = 0, datatype = "abundance", size = m)

# Se pueden calcular todos los ordenes de q al mismo tiempo

(out <- iNEXT(spider, q= c(0,1,2), datatype = "abundance", size = m))

### 4. Funciones gráficas ggiNEXT() #############

# La función ggiNEXT emplea ggplot2, usa como base un objeto tipo iNEXT y el resultado se puede manipular con opciones
# de ggplot2.

# Genera 3 tipos de curva:

# (1) R/E basado en tamaño de muestra (type = 1): Plotea los estimados de diversidad con intervalos de confianza (se = TRUE)
# al doble del tamaño del sitio de referencia.

# (2) Curva de completitud de muestra (type = 2): Con los intervalos de confianza (se = TRUE). Plotea la cobertura de muestra 
# respecto del tamaño de muestra para el mismo rango descrito en (1).

# (3) Curva basada en cobertura (type = 3): Esta curva plotea los estimados de diversidad con intervalos de confianza (se = TRUE)
# como una función de cobertura de muestra llevada al tamaño máximo de muestra descrito en (1).


## Argumentos facet.var("none", "order", "site" o "both")
# usado para crear un plot separado para cada valor especificado
# si, facet.var = "both" , podemos usar el argumento color.var("none", "order", "site" o "both"), para mostrar las curvas de color diferente para cada valor
# también puede usarse grey = TRUE para plotear figuras en blanco y negro


## Curva basada en tamaño de muestra (type = 1)

# separando sitios
ggiNEXT(out, type = 1, facet.var = "site")

# separando por orden de q
ggiNEXT(out, type = 1, facet.var = "order") 

## Curva basada en completitud de la muestra (type = 2)
ggiNEXT(out, type = 2)

## Curva basada en cobertura (type = 3)
ggiNEXT(out, type = 3, facet.var = "site") 


### 5. Otras opciones gráficas ###########

library(gridExtra)
library(grid)

# Blanco y negro
ggiNEXT(out, type = 1, facet.var = "order", grey = TRUE)

# Remover leyenda
ggiNEXT(out, type = 3, facet.var = "site") +
  theme(legend.position = "none")

# Usar temas, theme_bw()
ggiNEXT(out, type = 1, facet.var = "site") +
  theme_bw(base_size = 18)

# Cambiar la forma de los puntos

ggiNEXT(out, type = 1, facet.var = "site") +
  scale_shape_manual(values = c(19, 19, 19))

# Liberar la escala del eje Y

ggiNEXT(out, type = 1, facet.var = "order") +
  facet_wrap(~order, scales = "free")


## Nuevo gráfico
out <- iNEXT(spider, q = 0, datatype = "abundance")

g <- ggiNEXT(out, type = 1)
g

## Cambiando formas de puntos,  tipo de línea y color

g1 <- g + scale_shape_manual(values = c(11,12)) +
  scale_linetype_manual(values = c(1,2))

g2 <- g + scale_colour_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue"))


grid.arrange(g1, g2, ncol = 2)


# Para cambiar el tamaño del punto de referencia o la curva de R/E  se puede modificar el objeto ggplot.

# Panel izquierdo: cambiar el tamaño del punto de referencia a 10 (tamaño por defecto: 5)
gb3 <- ggplot_build(g)
gb3$data[[1]]$size <- 10
gt3 <- ggplot_gtable(gb3)
#grid.draw(gt3)

# Panel derecho: cambiar el tamaño de línea a 3 (tamaño por defecto es 1.5)

gb4 <- ggplot_build(g)
gb4$data[[2]]$size <- 3
gt4 <- ggplot_gtable(gb4)
#grid.draw(gt4)

grid.arrange(gt3, gt4, ncol = 2)

## Personalizar theme

# Panel izquierdo: cambiar a tema black and white
g5 <- g +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(legend.position = "bottom")

# Panel derecho: cambiar al thema clásico blanco y negro

g6 <- g +
  theme_classic() +
  theme(legend.position = "bottom")

grid.arrange(g5, g6, ncol = 2)

# Más temas

library(ggthemes)

g7 <- g + theme_hc(style = "darkunica") +
  scale_color_hc("darkunica")

g8 <- g + theme_economist() +
  scale_colour_economist() 

grid.arrange(g7, g8, ncol = 2)


# Tema Black-White

g9 <- g +
  theme_bw(base_size = 18) +
  scale_fill_grey(start = 0, end = .4) +
  scale_colour_grey(start = .2, end = .2) +
  theme(legend.position = "bottom",
        legend.title = element_blank())

g10 <- g +
  theme_tufte(base_size = 12) +
  scale_fill_grey(start = 0, end = .4) +
  scale_colour_grey(start = .2, end = .2) +
  theme(legend.position = "bottom",
        legend.title = element_blank ())

grid.arrange(g9, g10, ncol = 2)

#### 6. Dibujar la curva R/E ##########

df <- fortify(out, type = 1)
head(df)

df.point <- df[which(df$method =="observed"),]
df.line <- df[which(df$method != "observed"), ]

df.line$method <- factor(df.line$method,
                         c("interpolated", "extrapolated"),
                         c("interpolation", "extrapolation"))

ggplot(df, aes(x = x, y = y, colour = site)) +
  geom_point(aes(shape = site), size = 5, data = df.point) +
  geom_line(aes(linetype = method), lwd = 1.5, data = df.line) +
  geom_ribbon(aes(ymin = y.lwr, ymax = y.upr,
                  fill = site, color = NULL), alpha = 0.2) +
  labs(x = "Número de individuos", y = "Diversidad de especies") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        text = element_text(size = 18))

