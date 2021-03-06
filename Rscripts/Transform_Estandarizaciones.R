
########## Tramsformaciones y Estandarizaciones ###########

# Fuente: https://chrischizinski.github.io/SNR_R_Group/2016-08-10-Data-Transformations 



### Paquetes ##############

library(vegan)
library(tidyverse)
library(cowplot)

### Datos #################

df <- read.csv("Data/A-Amoco/Amoco-Cadiz.csv", row.names = 1, header = TRUE)
df <- as.data.frame(t(df))


### 1. Transformaciones (Monot�nicas)----------
# Prop�sito:

## Estad�stico
# Mejoran losupuestos de normalidad y linearidad, homogeneidad de varianza, etc.
# Hacen las unidades de las variables comparables cuando son medidas a escalas diferentes.

## Ecol�gica
# Hacen que las distancias ecol�gicas trabajen mejor
# Reducen el efecto de la cantidad total en las unidades de muestreo, para poner el foco en las cantidades relativas.
# Igualan (o alteran) la importancia relativa de las variables (p.ej. comunes y raras)
# Enfatizan en variables informativas (especies) a expensa de variables no informativas (especies).

## Datos sin transformar
(
  ST <- ggplot(df, aes(x = sp_1)) +
  geom_histogram()
 )

# Aplican a todo elemento de la matriz de datos
## Funciones monot�nicas (conservan el orden de las estaciones y no modifican las relaciones entre estaciones)
# Seg�n Clarke & Warwick (2016), en orden de intensidad de la transformaci�n

## Ra�z cuadrada
# Comprime valores altos y ampl�a valores bajos expresando valores en orden de magnitud. Efecto menos dram�tico que la transformaci�n logar�tmica, a menudo usado en datos de conteo (distribuci�n de Poisson)

Sqrt2 <- sqrt(df)

(
  T_2 <- ggplot(Sqrt2, aes(x = sp_1)) +
  geom_histogram()
)

## Doble ra�z cuadrada (ra�z cuarta)
Sqrt4 <- sqrt(sqrt(df))

# Similar a la transformaci�n de ra�z cuadrada. M�s severa que esta pero menos que la transformaci�n logar�tmica.

(
  T_4 <- ggplot(Sqrt4, aes(x = sp_1)) +
  geom_histogram()
)

## Logaritmo (Log + 1)
# Comprime valores altos y ampl�a valores bajos expresando valores en orden de magnitud. �til cuando hay un alto grado de variaci�n; proporci�n del m�s grande y del m�s peque�o > 10; datos muy positivamente sesgados.
Log <- log1p(df)


(
  T_L<- ggplot(Log, aes(x = sp_1)) +
  geom_histogram()
)

## Arcoseno de la ra�z cuadrada
# Comprime los extremos y ampl�a los valores intermedios, �til cuando existe sesgo positivo (Y tambi�n negativo)

prop.data <- df / apply(df,1,sum) # Primero estandariza dividiendo por el total del sitio (coloca los datos en un rango de 0 a 1)

# Funci�n de transformaci�n: 1ro. Ra�z cuadrada, 2do. funci�n arcoseno
acsin_trans<-function(x){
  x<- 2/pi *asin(sqrt(x))
  return(x)
}

Arcsin <- apply(prop.data, c(1,2), function(x) acsin_trans(x))
Arcsin


## Presencia - Ausencia
# Aplicable para datos de especies, m�s �til cuando hay poca informaci�n cuantitativa presente. Transformaci�n severa.
PA <- decostand(df, "pa")

(
  T_PA <- ggplot(PA, aes(x = sp_1)) +
    geom_histogram()
)
  
plot_grid(ST, T_2, T_4, T_L, T_PA)


### 1.1 REGLAS DE TRANSFORMACIONES MONOTONICAS ##############

# Seg�n Dr. Kevin McGarigals en su curso Applied Multivariate Analysis.

### 1. Usar transfomraciones log o ra�z cuadrada cuando hay sesgo positivo (2 o m�s ordenes de magnitud)

### 2. Usar transformaci�n arcoseno ra�z cuadrada cuando los datos son de proporciones.

### 3. Si se aplica a un conjunto de variables relacionadas emplear la misma transformaci�n, tal que las variables tengan el mismo escalamiento.

### 4. Considerar la transformaci�n Presencia - Ausencia, cuando: 
  # El % de ceros es mayor al 50%
  # El valor de variables significativas es menor de 10%
  # La diversidad beta es alta > 5. B = (Especies/promedio de especies) - 1

### 2. Estandarizaciones ----------

## Cu�ndo estandarizar?
# Estandarizar para igualar las condiciones de unidades de muestreo o variables altamente desiguales.
# Para representar mejor los patrones de inter�s.

## Qu� estandarizaci�n?
# Depende del objetivo (ajuste de muestra o variable) y t�cnica estad�stica (Ordenaci�n, cluster, etc?).
### C�lculos r�pidos (algunas estandarizaciones est�n basadas en c�lculos diversos de las columnas y filas del data.frame)





## Estandarizaci�n por filas
# Cuando la principal preocupaci�n es ajustar diferencias (abundancia total, diversidad) entre unidades de muestreo en orden de ponerlas en igualdad de condiciones.
# Cuando la atenci�n est� en el perfil dentro de la unidad de muestreo.
# Filas (Especies)
colSums(df) # Suma de filas
apply(df,2,max) # Valores m�ximos

## Estandarizaci�n por columnas
# Cuando la principal preocupaci�n es ajustar diferencias (varianzas, abundancia total, ubicuidad) entre variables (especies) en orden de igualar las condiciones.
# Cuando la atenci�n est� en el perfil entre unidades de muestreo.
# Columnas (Sitios)
rowSums(df) # Suma de filas
apply(df,1, max) # Valores m�ximos

# Ajustan los datos en base a un estad�stico de las filas o columnas (p.ej. max, sum, mean)

# MARGIN = 1 ## Por filas o sitios
# MARGIN = 2 ## Por columnas o especies

## Estadarizaci�n por z-scores 

# Convierte los datos z-scores (media = 0), (varianza = 1)
# Com�mente usado para poner variables en igualdad de condiciones
# Esencial cuando las variables tienen diferentes escalas o unidades de medida

std_z_score <- decostand(df, "standardize")

## Estandarizaci�n dividiendo por el total de columnas (especies)
# Com�nmente usado con datos de especies para ajustar abundancias desiguales entre especies.
# Iguala �reas bajo los perfiles de respuesta de curva de especies.

Total_2 <- decostand(df, "total", MARGIN = 2)
colSums(Total)

## Estandarizaci�n dividiendo por el m�ximo de columnas (especies)
# Similar al total de columnas, excepto que: iguala los altos de los picos de las curvas de especies. 
# Basados en valores extremos que pueden introducir ruido.
# Pueden exhacerbar importancia de especies raras.

Max_2 <- decostand(df, "max", MARGIN = 2)

## Estandarizaci�n dividiendo por el total de filas (sitios)
# Com�nmente usado con datos de especies para ajustar abundancias desiguales entre unidades de muestreo.
# Iguala areas bajo las ruvas de perfil de muestras.
# Enfatiza la abundancia relativa dentro de la unidad de muestreo.
# Perfiles de abundancia relativa de las muestras son independientes.

Total_1 <- decostand(df, "total", MARGIN = 1)
rowSums(Total)

## Estandarizaci�n dividiendo por el m�ximo de filas (sitios)
# Similar al total de filas; excepto:
# Iguala la altura de los picos de los perfiles de muestras.
# Basados en valore extremos que pueden reducir el ruido

Max_1 <- decostand(df, "max", MARGIN = 1)



## Estandarizaci�n de Wisconsin

# Palacio 2018 indica que la doble estandarizaci�n de Wisconsin, cada elemento se divide por el m�ximo de su columna y luego por el total de sus filas correspondientes de la nueva matriz (Cottam et al. 1978). Los resultados var�an entre 0 y 1.
# La funci�n metaMDS del paquete vegan, efect�a la estadarizaci�n de Wisconsin cuando la abundancia es mayor a nueve y adem�s la transformaci�n ra�z cuadrada cuando la abundancia es mayor a 50.
# Atractivo, pero tiene el costo de disminuir el significado intuitivo de los valores individuales.

Wis_S <- wisconsin(df)


### REGLAS ##########
# El efecto de la estandarizaci�n depende de la variabilidad en las columnas/filas
# Estandarizaci�n de filas para, datos de especies comunmente:
  # Distancia euclidiana (Normalizaci�n por filas)
  # Distancia chi-cuadrado (Estandarizaci�n Chi cuadrado) para CA y CCA
  # Perfil de distancia de especies (Estandarizaci�n por total de filas)
  # Distancia de Hellinger (Distancia de Hellinger)

# Estandarizaci�n de columnas para:
  # Para igualar variables medidas en diferentes unidades y escalas
  # Estandarizaci�n por columnas z-score
  # Normalizaci�n por columnas
  # Estandarizaci�n total de columnas
  # Estandarizaci�n por rango

# Algunas estandarizaciones pueden ser irrelevantes en estos casos:
  # PCA de datos estandarizados por columna
  # CA est� construido sobre la base de la estandarizaci�n chi-cuadrado

# No hay base te�rica para la selecci�n de un m�todo de estandarizaci�n u otro, deben de ser justificados sobre la base biol�gica o mediante an�lisis de sensibilidad.

