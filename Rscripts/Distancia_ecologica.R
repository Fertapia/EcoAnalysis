
########## Distancia ecológica (di-) similaridades ###########

# Fuente: McCune & Grace (2002). Analysis of ecological communities. Chapter 6: Distance Measures
#         Borcard  et al. (2018). Numerical Ecology with R

# La semejanza (resemblance) puede ser medida ya sea como distancia (disimilaridad) o como similaridad.
# La mayoría de las medidas de distancia pueden ser convertidas en similaridades y viceversa.
# Las medidas de distancia pueden ser aplicadas a datos binarios (presencia-ausencia) o cuantitativos.
# Puede calcularse las distancias entre filas o columnas. Con datos comunitarios, los datos analizadados son entre unidades de muestreo en un espacio de especies o entre especies en un espacio de muestreo.

### Paquetes ##############

library(vegan)
library(tidyverse)
library(ade4)
library(adespatial)

install.packages("adespatial")

### Datos #################

df <- read.csv("Data/A-Amoco/Amoco-Cadiz.csv", row.names = 1, header = TRUE)
df <- as.data.frame(t(df))

####### IMPORTANTE ############

# Es importante saber el dominio de datos aceptables para cada medida de distancia. Muchas medidas no son compatibles con números negativos, y otras asumen que los datos son proporciones que van de 0 a 1.


# Es importante tener en cuenta las principales categorías de medidas para elegir la más apropiada:

# Nota: Se emplean las palabras "meadida" ("measure"), "índice" ("index"), y "coeficiente" ("coefficient"); como sinónimos al referirse a las cantidades usadas para compara pares de objetos o variables.


# Antes de elegir una medida responder:

# 1. Estás comparando sitios u objetos (Q-mode) o variables especies (R-mode)?

# 2. Estás tratando con datos de especies (coeficientes asimétricos) u otros tipos de variables (coeficientes simétricos)?

# 3. Tus datos son binarios (coeficientes binarios) o cuantitativos (coeficientes cuantitativos) o de otros tipos (coeficientes ordinales y especiales)?

#### Modo Q: Comparar objetos (sitios). Disimilaridad o similaridad (distancia Euclidiana, Jaccard). ###############################

## Problema del doble cero: En datos de especies (abundancias o presencia-ausencia), el significado del cero no implica ausencia absoluta (la especie puede estar presente pero no haber sido observada). La clave es:

# 1. La ausencia de especies de ambos sitios (doble cero) no puede ser contada como un indicador de similitud entre sitios.

# 2. El número de no-interpretables doble ceros depende del número de especies y aumenta considerablemente con las especies raras.

# Las medidas que consideran el doble cero como similitud entre sitio son SIMÉTRICAS y las otras son ASIMÉTRICAS.
# Es preferible emplear medidas asimétricas debido a la naturaleza de los datos de abundancia de especies.

########## 1. Q-mode (Semi-) Cuantitativos Asimétricos ##################

### Disimilariada Bray-Curtis (alias Porcentaje de similaridad) D-14 #############

# Puede ser calculado directamente de datos crudos, también a partir de transformaciones.
# Da la misma importancia a diferncias absolutas en la abundancia indiferente del orden de magnitud:
# la diferencia de 5 individuos de 3 a 8 tiene el mismo valor que de 6203 a 6208.

# En datos crudos
spe.db <- vegdist(df) # method = "bray" (por defecto)
head(spe.db)

# En datos log-transformados
spe.dbln <- vegdist(log1p(df))
head(spe.dbln)

#### Distancia de Chord #######
# Distancia euclidiana calculada en sitios normalizados de longitud 1 (transformación de chord).
# La función decostand() de vegan es hecha mediante el argumento normalize.
# El paquete adespatial hace el cálculo directamente con la funciín dist.ldc() y el argumento chord.

spe.dc <- dist.ldc(df, "chord")
head(spe.dc)

# Alternativa de dos pasos en vegan
spe.norm <- decostand(df, "nor")
spe.dc <- dist(spe.norm)
head (spe.dc)

### Distancia de Hellinger #######
# Es la distancia Euclidiana entre sitios donde los valores de abundancia son primero dividios por la abundancia total, y el resultado es transformado por raíz cuadrada. 
# es obtenido por la función decostand(), argumento hellinger.

spe.dh <- dist.ldc(df) # Distancia de Hellinger por defecto.
head(spe.dh)

### log-chord distance ###########
# Es la distancia de chord a la que se le aplicó una transformación ln(y+1).
# dist.ldc() del paquete adespatial argumento log.chord.

spe.logchord <- dist.ldc(df, "log.chord")
head(spe.logchord)


########## 2. Q-mode (Semi-) Cuantitativos Simétricos##################

### Distancia Euclidiana ##########

# Calculada mediante la fórmula de pitágoras, a partir de puntos de sitios posicionados en un espacio p-dimensional llamado espacio Euclidiano.
# No tiene límite superior, y está fuertemente influenciado por la escala del descriptor. Su uso en datos curdos está restringido a datos homogéneos. 
# Se debe estandarizar en z-scores.
# Se recomienda su uso en variable ambientales. Nota: Sólo para el ejemplo se emplearán los datos del data frame df.

# En datos crudos
Euc.df <- dist(df)

#transformados
Euc.df2 <- dist(scale(df))


########## 3. Q-mode Cualitativos Asimétricos##################

# Relevante sólo cuando los datos disponibles son binarios, o cuando las abundancias son irrelevantes, o cuando los datos contienen valores cuantitativos de calidad incierta. Los análisis son realizados en datos de presencia-ausencia (1-0).

# Nota: No es necesario transformar previamente los datos a binario, las funciones realizan esta transformación automáticamente, la función vegdist requiere el arugmento binary= TRUE.


#### Disimilaridad de Jaccard ######
# Es el radio de similaridad entre el número de dobles 1s y el número de especies, excluyendo los doble ceros en el par considerado.
# un jaccard de 0.25 significa que el 25% del total de especies observadas en dos sitios estuvieron presentes en ambos sitios y el 75% sólo en uno.

spe.dj <- vegdist(df, "jac", binary = TRUE)
head(spe.dj)

#### Disimilaridad de Sorensen #####
# Brinda doble peso al número de dobles 1s, su recíproco (complemento a 1) es equivalente al porcentaje de diferencia (Bray Curtis).

spe.ds2 <- vegdist(df, method = "bray", binary = TRUE)

#### Disimilaridad de Ochiai ############

# Relacionado a chord y Hellinger, calcula la distancia en presencia ausencia seguida por una división por raíz de 2. 

spe.och <- dist.ldc(df, "ochiai")
head(spe.och)

########## 3. Q-mode Tipos Mixtos categóricos y cuantitativos ##################


###### Similiaridad de Gower ##########

# Puede manejar variabes de varios tipos matemáticos. Cada una recibe un tratamiento acorde a su categoría.

# La (di)similaridad entre dos objetos es obtenida promediando las  (di) similaridades de todas las variable separadamente. 

# Gower es una medida simétrica.

# Cuando una variable es declarada como un factor en el data frame, aplica una regla simple de coincidencias (0 si es diferente, 1 si son iguales).

# La función daisy() del paquete cluster calcula la medida. Evitar usar vegdist(method = "gower"), no está diseñada para variables multiclase.

### Creando un conjunto de datos ficticio

# Variable aleatoria normal con media cero y desviación estándar 1 

var.g1 <- rnorm(30, 0, 1)

# Aleatorios de 0 a 5 
var.g2 <- runif(30, 0, 5)

# Factor de 3 niveles (10 objetos cada uno) 

var.g3 <- gl(3, 10, labels = c("A", "B", "C"))

# Factor con 2 niveles, ortogonal a var.g3
var.g4 <- gl(2, 5, 30, labels = c("D", "E"))

dat2 <- data.frame(var.g1, var.g2, var.g3, var.g4)
summary(data2)

### Calculando la disimilaridad de Gower

dat2.S15 <- daisy(dat2, "gower")

#### Modo R: Comparar variables o descriptores (especies). Depedencia (covarianza  o coeficiente de correlación). #####################

##### R-mode para datos de abundancias de especies #################


## Distancia Chi-cuadrado ##########

## trasponer

df.t <- t(df)

# Pretransformación Chi-cuadrado seguida por Distancia Euclidiana

spe.t.chi <- decostand(spe.t, "chi.square")
spe.t.D16 <- dist(spe.t.chi)

##### R mode para datos de especies  presencia-ausencia #########

# Jaccard, Sorensen y Ochiai pueden ser usados en mod R.

# Índice Jaccard presencia ausencia

spe.t.S7 <- vegdist(df.t, "jaccard", binary = TRUE)

##### R mode para datods cuantitativos y ordinales #########

### Coeficiente de correlación de Pearson #########

# Recomendado cuando las variables no son dimensionalmente homogéneas, ya que la correlación r es de hecho la covarianza calculada en variables estandarizadas.

# Correlación de pearson

Pearson <- cor(df) # método por defecto = "pearson"

# La comparación entre variables ordinales o variables cuantitativas que pueden ser monotónicas pero no linearmente relacionadas, puede ser logrado usando un correlación de clasificación como Spearman o Kendall.

Kendall <- cor(df, method = "kendall")

