######## Análisis de Coordenadas Principales PCoA #############

#### 0. Nota ########

# El análisis de Coordenadas Principales (PCoA) es una función de paquete phyloseq
# Este paquete requiere de un objeto de clase phyloseq o "phyloseq-class".


#### 1. Paquetes ######

library("phyloseq"); packageVersion("phyloseq")
library("tidyverse"); packageVersion("tidyverse")
library("vegan"); packageVersion("vegan")
library("plyr"); packageVersion("plyr")
library("ggforce"); packageVersion("ggforce")

#### 2. Datos #######

## Directorio de trabajo

setwd(paste0(getwd(), "/Data/vegan_dune"))

# Crear objeto phyloseq

## Matriz de especies (OTUs) 

df <- read.csv("dune.csv", row.names = 1)
otu_mat <- as.matrix(df)
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

## Matriz de taxas 

df_tax <- read.csv("dune.taxon.csv", row.names = 1)
tax_mat <- as.matrix(df_tax)
TAX <- tax_table(tax_mat)

## Matriz de muestras

df_env <- read.csv("dune.env.csv", row.names = 1)
samples <- sample_data(df_env)

## Objeto phyloseq

Phy_data <- phyloseq(OTU, TAX, samples)
# Los datos para trabajar están dentro de Phy_data


#### 3. Calculamos el PCoA ###############

# Cálculo del ordenamiento

Ord <- ordinate(Phy_data, "PCoA", "bray")

Ord1 <- ordinate(Phy_data, "CCA", "bray")

# Valores del modelo
Ordsum <- Ord$values

# Diagrama básico

(ordplot <- plot_ordination(Phy_data, Ord, type = "samples", shape = "Management"))

# Extraemos los datos

orddata <- ordplot$data # para visualizar los resultados del ordenamiento

# Personalizacion de datos

ordplot +
  geom_point(size = 3) +
  scale_fill_brewer("Manejo", palette = "Paired") +
  scale_shape_manual("Manejo", values = c(15,16,17, 18)) +
  geom_mark_ellipse(aes(fill = Management)) + # Función del paquete ggforce
  theme_bw()

