
########## Distancia ecol�gica (di-) similaridades ###########

# Fuente: McCune & Grace (2002). Analysis of ecological communities. Chapter 6: Distance Measures
#         Borcard  et al. (2018). Numerical Ecology with R

# La semejanza (resemblance) puede ser medida ya sea como distancia (disimilaridad) o como similaridad.
# La mayor�a de las medidas de distancia pueden ser convertidas en similaridades y viceversa.
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

# Es importante saber el dominio de datos aceptables para cada medida de distancia. Muchas medidas no son compatibles con n�meros negativos, y otras asumen que los datos son proporciones que van de 0 a 1.


# Es importante tener en cuenta las principales categor�as de medidas para elegir la m�s apropiada:

# Nota: Se emplean las palabras "meadida" ("measure"), "�ndice" ("index"), y "coeficiente" ("coefficient"); como sin�nimos al referirse a las cantidades usadas para compara pares de objetos o variables.


# Antes de elegir una medida responder:

# 1. Est�s comparando sitios u objetos (Q-mode) o variables especies (R-mode)?

# 2. Est�s tratando con datos de especies (coeficientes asim�tricos) u otros tipos de variables (coeficientes sim�tricos)?

# 3. Tus datos son binarios (coeficientes binarios) o cuantitativos (coeficientes cuantitativos) o de otros tipos (coeficientes ordinales y especiales)?

#### Modo Q: Comparar objetos (sitios). Disimilaridad o similaridad (distancia Euclidiana, Jaccard). ###############################

## Problema del doble cero: En datos de especies (abundancias o presencia-ausencia), el significado del cero no implica ausencia absoluta (la especie puede estar presente pero no haber sido observada). La clave es:

# 1. La ausencia de especies de ambos sitios (doble cero) no puede ser contada como un indicador de similitud entre sitios.

# 2. El n�mero de no-interpretables doble ceros depende del n�mero de especies y aumenta considerablemente con las especies raras.

# Las medidas que consideran el doble cero como similitud entre sitio son SIM�TRICAS y las otras son ASIM�TRICAS.
# Es preferible emplear medidas asim�tricas debido a la naturaleza de los datos de abundancia de especies.

########## 1. Q-mode (Semi-) Cuantitativos Asim�tricos ##################

### Disimilariada Bray-Curtis (alias Porcentaje de similaridad) D-14 #############

# Puede ser calculado directamente de datos crudos, tambi�n a partir de transformaciones.
# Da la misma importancia a diferncias absolutas en la abundancia indiferente del orden de magnitud:
# la diferencia de 5 individuos de 3 a 8 tiene el mismo valor que de 6203 a 6208.

# En datos crudos
spe.db <- vegdist(df) # method = "bray" (por defecto)
head(spe.db)

# En datos log-transformados
spe.dbln <- vegdist(log1p(df))
head(spe.dbln)

#### Distancia de Chord #######
# Distancia euclidiana calculada en sitios normalizados de longitud 1 (transformaci�n de chord).
# La funci�n decostand() de vegan es hecha mediante el argumento normalize.
# El paquete adespatial hace el c�lculo directamente con la funci�n dist.ldc() y el argumento chord.

spe.dc <- dist.ldc(df, "chord")
head(spe.dc)

# Alternativa de dos pasos en vegan
spe.norm <- decostand(df, "nor")
spe.dc <- dist(spe.norm)
head (spe.dc)

### Distancia de Hellinger #######
# Es la distancia Euclidiana entre sitios donde los valores de abundancia son primero dividios por la abundancia total, y el resultado es transformado por ra�z cuadrada. 
# es obtenido por la funci�n decostand(), argumento hellinger.

spe.dh <- dist.ldc(df) # Distancia de Hellinger por defecto.
head(spe.dh)

### log-chord distance ###########
# Es la distancia de chord a la que se le aplic� una transformaci�n ln(y+1).
# dist.ldc() del paquete adespatial argumento log.chord.

spe.logchord <- dist.ldc(df, "log.chord")
head(spe.logchord)


########## 2. Q-mode (Semi-) Cuantitativos Sim�tricos##################

### Distancia Euclidiana ##########

# Calculada mediante la f�rmula de pit�goras, a partir de puntos de sitios posicionados en un espacio p-dimensional llamado espacio Euclidiano.
# No tiene l�mite superior, y est� fuertemente influenciado por la escala del descriptor. Su uso en datos curdos est� restringido a datos homog�neos. 
# Se debe estandarizar en z-scores.
# Se recomienda su uso en variable ambientales. Nota: S�lo para el ejemplo se emplear�n los datos del data frame df.

# En datos crudos
Euc.df <- dist(df)

#transformados
Euc.df2 <- dist(scale(df))


########## 3. Q-mode Cualitativos Asim�tricos##################

# Relevante s�lo cuando los datos disponibles son binarios, o cuando las abundancias son irrelevantes, o cuando los datos contienen valores cuantitativos de calidad incierta. Los an�lisis son realizados en datos de presencia-ausencia (1-0).

# Nota: No es necesario transformar previamente los datos a binario, las funciones realizan esta transformaci�n autom�ticamente, la funci�n vegdist requiere el arugmento binary= TRUE.


#### Disimilaridad de Jaccard ######
# Es el radio de similaridad entre el n�mero de dobles 1s y el n�mero de especies, excluyendo los doble ceros en el par considerado.
# un jaccard de 0.25 significa que el 25% del total de especies observadas en dos sitios estuvieron presentes en ambos sitios y el 75% s�lo en uno.

spe.dj <- vegdist(df, "jac", binary = TRUE)
head(spe.dj)

#### Disimilaridad de Sorensen #####
# Brinda doble peso al n�mero de dobles 1s, su rec�proco (complemento a 1) es equivalente al porcentaje de diferencia (Bray Curtis).

spe.ds2 <- vegdist(df, method = "bray", binary = TRUE)

#### Disimilaridad de Ochiai ############

# Relacionado a chord y Hellinger, calcula la distancia en presencia ausencia seguida por una divisi�n por ra�z de 2. 

spe.och <- dist.ldc(df, "ochiai")
head(spe.och)

########## 3. Q-mode Tipos Mixtos categ�ricos y cuantitativos ##################


###### Similiaridad de Gower ##########

# Puede manejar variabes de varios tipos matem�ticos. Cada una recibe un tratamiento acorde a su categor�a.

# La (di)similaridad entre dos objetos es obtenida promediando las  (di) similaridades de todas las variable separadamente. 

# Gower es una medida sim�trica.

# Cuando una variable es declarada como un factor en el data frame, aplica una regla simple de coincidencias (0 si es diferente, 1 si son iguales).

# La funci�n daisy() del paquete cluster calcula la medida. Evitar usar vegdist(method = "gower"), no est� dise�ada para variables multiclase.

### Creando un conjunto de datos ficticio

# Variable aleatoria normal con media cero y desviaci�n est�ndar 1 

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

#### Modo R: Comparar variables o descriptores (especies). Depedencia (covarianza  o coeficiente de correlaci�n). #####################

##### R-mode para datos de abundancias de especies #################


## Distancia Chi-cuadrado ##########

## trasponer

df.t <- t(df)

# Pretransformaci�n Chi-cuadrado seguida por Distancia Euclidiana

spe.t.chi <- decostand(spe.t, "chi.square")
spe.t.D16 <- dist(spe.t.chi)

##### R mode para datos de especies  presencia-ausencia #########

# Jaccard, Sorensen y Ochiai pueden ser usados en mod R.

# �ndice Jaccard presencia ausencia

spe.t.S7 <- vegdist(df.t, "jaccard", binary = TRUE)

##### R mode para datods cuantitativos y ordinales #########

### Coeficiente de correlaci�n de Pearson #########

# Recomendado cuando las variables no son dimensionalmente homog�neas, ya que la correlaci�n r es de hecho la covarianza calculada en variables estandarizadas.

# Correlaci�n de pearson

Pearson <- cor(df) # m�todo por defecto = "pearson"

# La comparaci�n entre variables ordinales o variables cuantitativas que pueden ser monot�nicas pero no linearmente relacionadas, puede ser logrado usando un correlaci�n de clasificaci�n como Spearman o Kendall.

Kendall <- cor(df, method = "kendall")

