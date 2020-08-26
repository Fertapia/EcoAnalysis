
########## Tramsformaciones y Estandarizaciones ###########

# Fuente: https://chrischizinski.github.io/SNR_R_Group/2016-08-10-Data-Transformations 



### Paquetes ##############

library(vegan)
library(tidyverse)
library(cowplot)

### Datos #################

df <- read.csv("Data/A-Amoco/Amoco-Cadiz.csv", row.names = 1, header = TRUE)
df <- as.data.frame(t(df))


### 1. Transformaciones (Monotónicas)----------
# Propósito:

## Estadístico
# Mejoran losupuestos de normalidad y linearidad, homogeneidad de varianza, etc.
# Hacen las unidades de las variables comparables cuando son medidas a escalas diferentes.

## Ecológica
# Hacen que las distancias ecológicas trabajen mejor
# Reducen el efecto de la cantidad total en las unidades de muestreo, para poner el foco en las cantidades relativas.
# Igualan (o alteran) la importancia relativa de las variables (p.ej. comunes y raras)
# Enfatizan en variables informativas (especies) a expensa de variables no informativas (especies).

## Datos sin transformar
(
  ST <- ggplot(df, aes(x = sp_1)) +
  geom_histogram()
 )

# Aplican a todo elemento de la matriz de datos
## Funciones monotónicas (conservan el orden de las estaciones y no modifican las relaciones entre estaciones)
# Según Clarke & Warwick (2016), en orden de intensidad de la transformación

## Raíz cuadrada
# Comprime valores altos y amplía valores bajos expresando valores en orden de magnitud. Efecto menos dramático que la transformación logarítmica, a menudo usado en datos de conteo (distribución de Poisson)

Sqrt2 <- sqrt(df)

(
  T_2 <- ggplot(Sqrt2, aes(x = sp_1)) +
  geom_histogram()
)

## Doble raíz cuadrada (raíz cuarta)
Sqrt4 <- sqrt(sqrt(df))

# Similar a la transformación de raíz cuadrada. Más severa que esta pero menos que la transformación logarítmica.

(
  T_4 <- ggplot(Sqrt4, aes(x = sp_1)) +
  geom_histogram()
)

## Logaritmo (Log + 1)
# Comprime valores altos y amplía valores bajos expresando valores en orden de magnitud. Útil cuando hay un alto grado de variación; proporción del más grande y del más pequeño > 10; datos muy positivamente sesgados.
Log <- log1p(df)


(
  T_L<- ggplot(Log, aes(x = sp_1)) +
  geom_histogram()
)

## Arcoseno de la raíz cuadrada
# Comprime los extremos y amplía los valores intermedios, útil cuando existe sesgo positivo (Y también negativo)

prop.data <- df / apply(df,1,sum) # Primero estandariza dividiendo por el total del sitio (coloca los datos en un rango de 0 a 1)

# Función de transformación: 1ro. Raíz cuadrada, 2do. función arcoseno
acsin_trans<-function(x){
  x<- 2/pi *asin(sqrt(x))
  return(x)
}

Arcsin <- apply(prop.data, c(1,2), function(x) acsin_trans(x))
Arcsin


## Presencia - Ausencia
# Aplicable para datos de especies, más útil cuando hay poca información cuantitativa presente. Transformación severa.
PA <- decostand(df, "pa")

(
  T_PA <- ggplot(PA, aes(x = sp_1)) +
    geom_histogram()
)
  
plot_grid(ST, T_2, T_4, T_L, T_PA)


### 1.1 REGLAS DE TRANSFORMACIONES MONOTONICAS ##############

# Según Dr. Kevin McGarigals en su curso Applied Multivariate Analysis.

### 1. Usar transfomraciones log o raíz cuadrada cuando hay sesgo positivo (2 o más ordenes de magnitud)

### 2. Usar transformación arcoseno raíz cuadrada cuando los datos son de proporciones.

### 3. Si se aplica a un conjunto de variables relacionadas emplear la misma transformación, tal que las variables tengan el mismo escalamiento.

### 4. Considerar la transformación Presencia - Ausencia, cuando: 
  # El % de ceros es mayor al 50%
  # El valor de variables significativas es menor de 10%
  # La diversidad beta es alta > 5. B = (Especies/promedio de especies) - 1

### 2. Estandarizaciones ----------

## Cuándo estandarizar?
# Estandarizar para igualar las condiciones de unidades de muestreo o variables altamente desiguales.
# Para representar mejor los patrones de interés.

## Qué estandarización?
# Depende del objetivo (ajuste de muestra o variable) y técnica estadística (Ordenación, cluster, etc?).
### Cálculos rápidos (algunas estandarizaciones están basadas en cálculos diversos de las columnas y filas del data.frame)





## Estandarización por filas
# Cuando la principal preocupación es ajustar diferencias (abundancia total, diversidad) entre unidades de muestreo en orden de ponerlas en igualdad de condiciones.
# Cuando la atención está en el perfil dentro de la unidad de muestreo.
# Filas (Especies)
colSums(df) # Suma de filas
apply(df,2,max) # Valores máximos

## Estandarización por columnas
# Cuando la principal preocupación es ajustar diferencias (varianzas, abundancia total, ubicuidad) entre variables (especies) en orden de igualar las condiciones.
# Cuando la atención está en el perfil entre unidades de muestreo.
# Columnas (Sitios)
rowSums(df) # Suma de filas
apply(df,1, max) # Valores máximos

# Ajustan los datos en base a un estadístico de las filas o columnas (p.ej. max, sum, mean)

# MARGIN = 1 ## Por filas o sitios
# MARGIN = 2 ## Por columnas o especies

## Estadarización por z-scores 

# Convierte los datos z-scores (media = 0), (varianza = 1)
# Comúmente usado para poner variables en igualdad de condiciones
# Esencial cuando las variables tienen diferentes escalas o unidades de medida

std_z_score <- decostand(df, "standardize")

## Estandarización dividiendo por el total de columnas (especies)
# Comúnmente usado con datos de especies para ajustar abundancias desiguales entre especies.
# Iguala áreas bajo los perfiles de respuesta de curva de especies.

Total_2 <- decostand(df, "total", MARGIN = 2)
colSums(Total)

## Estandarización dividiendo por el máximo de columnas (especies)
# Similar al total de columnas, excepto que: iguala los altos de los picos de las curvas de especies. 
# Basados en valores extremos que pueden introducir ruido.
# Pueden exhacerbar importancia de especies raras.

Max_2 <- decostand(df, "max", MARGIN = 2)

## Estandarización dividiendo por el total de filas (sitios)
# Comúnmente usado con datos de especies para ajustar abundancias desiguales entre unidades de muestreo.
# Iguala areas bajo las ruvas de perfil de muestras.
# Enfatiza la abundancia relativa dentro de la unidad de muestreo.
# Perfiles de abundancia relativa de las muestras son independientes.

Total_1 <- decostand(df, "total", MARGIN = 1)
rowSums(Total)

## Estandarización dividiendo por el máximo de filas (sitios)
# Similar al total de filas; excepto:
# Iguala la altura de los picos de los perfiles de muestras.
# Basados en valore extremos que pueden reducir el ruido

Max_1 <- decostand(df, "max", MARGIN = 1)



## Estandarización de Wisconsin

# Palacio 2018 indica que la doble estandarización de Wisconsin, cada elemento se divide por el máximo de su columna y luego por el total de sus filas correspondientes de la nueva matriz (Cottam et al. 1978). Los resultados varían entre 0 y 1.
# La función metaMDS del paquete vegan, efectúa la estadarización de Wisconsin cuando la abundancia es mayor a nueve y además la transformación raíz cuadrada cuando la abundancia es mayor a 50.
# Atractivo, pero tiene el costo de disminuir el significado intuitivo de los valores individuales.

Wis_S <- wisconsin(df)


### REGLAS ##########
# El efecto de la estandarización depende de la variabilidad en las columnas/filas
# Estandarización de filas para, datos de especies comunmente:
  # Distancia euclidiana (Normalización por filas)
  # Distancia chi-cuadrado (Estandarización Chi cuadrado) para CA y CCA
  # Perfil de distancia de especies (Estandarización por total de filas)
  # Distancia de Hellinger (Distancia de Hellinger)

# Estandarización de columnas para:
  # Para igualar variables medidas en diferentes unidades y escalas
  # Estandarización por columnas z-score
  # Normalización por columnas
  # Estandarización total de columnas
  # Estandarización por rango

# Algunas estandarizaciones pueden ser irrelevantes en estos casos:
  # PCA de datos estandarizados por columna
  # CA está construido sobre la base de la estandarización chi-cuadrado

# No hay base teórica para la selección de un método de estandarización u otro, deben de ser justificados sobre la base biológica o mediante análisis de sensibilidad.

