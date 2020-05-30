
############## Curvas de comparación de Abundancia - Biomasa (ABC curves) #########################

######## Preparando datos ###########

### Extraemos los datos crudos:

# Abundancias
#df <- read.csv("~/Ecoanalysis/Data/G-Garroch/clma.csv", row.names=1, check.names=FALSE)
#df <- t(df)

# Biomasa
#df_bm <-  read.csv("~/Ecoanalysis/Data/G-Garroch/clmb.csv", row.names=1, check.names = FALSE)
#df_bm <- t(df_bm)


### Las matrices de datos biológicos usualmente presentan la forma de Sitios (Filas) x Especies (Columnas):

### El paquete forams requiere de dos columnas: N: Abundancias y B: Biomasa por cada especie, las cuales constituyen las filas. ver data(NB) del paquete forams

## La función melt transforma la base de datos en un formato de 3 columnas: Sitio `Site`, Especie `Species` y Abundancia `N` o Biomasa `Biomass` según corresponda.

# Abundancia:
df <- melt(df)

df <- df %>%  
  rename(Site = Var1, Species = Var2, N = value) %>% 
  select(Site, Species, N)

# Biomasa:

df_bm <- melt(df_bm)

df_bm <- df_bm %>% 
  rename(Site = Var1, Species = Var2, Biomass = value) %>% 
  select(Site, Species, Biomass) 

# Combinamos ambas matrices en la matriz `df_abc`:

df_abc <- inner_join(df, df_bm)

df_abc <- df_abc %>% 
  filter(N !=0)

str(df_abc)

#####################################


######## Forma básica: Para una estación ############

# Seleccionamos una estación en este caso la estación S1

abc_S1 <- df_abc %>% 
  filter(Site == "S1") %>% 
  select(N, Biomass)

abc_S1


### La función abc del paquete forams calcula los ranking y las abundancias acumuladas

(abc_S1 <- abc(abc_S1))

# Abundancias y biomasas acumuladas:

abc_S1@abc

# Estadístico W y sus intervalos de confianza

abc_S1@W.Stat

# Diagrama por defecto según las R básico

plot(abc_S1)

### Diagrama personalizable en ggplot2

## Previos

# Guardamos los datos en objetos df_abc y W

df_abc <- abc_S1@abc

W <- abc_S1@W.Stat[2]

# Calculamos los rangos en escala logarítmica `ln` para biomasa y abundancia

x.bio = log(1:length(df_abc$Accum.Biomass))
x.ab = log(1:length(df_abc$Accum.Abund))

# Etiquetas del eje x:
brks <- 0:ceiling(log(length(df_abc$Accum.Abund)))

# Tipos de línea
line_types <- c("Abundancia"=1,"Biomasa"= 3)

## Diagrama

ggplot(data = df_abc) +
  geom_line(aes(x = x.bio, y = sort(Accum.Biomass), linetype= "Biomasa"), lwd = 1) +
  geom_line(aes(x = x.ab, y = sort(Accum.Abund), linetype = "Abundancia"), lwd =1) +
  geom_ribbon(aes( ymin = sort(Accum.Abund), ymax = sort(Accum.Biomass), x = x.bio), fill = "grey12", alpha = 0.1) +
  scale_linetype_manual(name = " ", values=line_types) +
  theme_classic() +
  xlab(expression('Rango de especies'~(Log[e]~Scale))) +
  ylab("Dominancia acumulada %") +
  scale_x_continuous(breaks = brks, 
                     labels = round(exp(0:ceiling(log(length(df_abc$Accum.Abund)))))
  ) +
  scale_y_continuous(limits = c(0,101), 
                     breaks = seq(0,100,10)) +
  geom_text(aes(x = 0, y = 100, label = paste('W =',
                                              W, 
                                              sep=' ')),
            vjust = "inward", hjust = "inward")

###################################

######### Diagramas para estaciones múltiples #######################
str(df_abc)
Sit_level <- c("S1", "S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12")
#### Preparamos datos 

## N de especies por estación
Lt <- df_abc %>% 
  split(.$Site) %>% 
  map(~(length(.$N))) %>% 
  map_dfr(~as.data.frame((log(1:.))), .id = "Sitio_")

Lt <- rename(Lt, Lt = "(log(1:.))")


## Matriz de abundancias acumulada
ABC_abun <- df_abc %>%
  select(Site,N, Biomass) %>% 
  split(.$Site) %>% 
  map( ~ (.x %>% select(-Site))) %>% 
  map(~abc(.)) %>% 
  map(~(.@abc)) %>% 
  map(~sort(.$Accum.Abund)) %>% 
  map_dfr(~(as.data.frame(.)), .id = "Sitio") %>% 
  rename(Accum.Abund =".")

ABC_bio <- df_abc %>%
  select(Site,N, Biomass) %>% 
  split(.$Site) %>% 
  map( ~ (.x %>% select(-Site))) %>% 
  map(~abc(.)) %>% 
  map(~(.@abc)) %>% 
  map(~sort(.$Accum.Biomass)) %>% 
  map_dfr(~(as.data.frame(.)), .id = "Sitio__") %>% 
  rename(Accum.Biomass =".")


## unimos ABC y Lt
ABC <- cbind(ABC_abun, ABC_bio, Lt)

ABC <- select(ABC, -c("Sitio_", "Sitio__"))
ABC$Sitio <- factor(ABC$Sitio, levels = Sit_level)

## Valores de W por estación
Ws <- df_abc %>%
  select(Site,N, Biomass) %>% 
  split(.$Site) %>% 
  map( ~ (.x %>% select(-Site))) %>% 
  map(~abc(.)) %>% 
  map_dfr(~(.@W.Stat[2]))

Ws <- as.data.frame(t(Ws))
Ws <- Ws %>% 
  tibble::rownames_to_column() %>% 
  rename(Sitio = rowname, W = V1)

Ws$Sitio <- factor(Ws$Sitio, levels = Sit_level)
#### Diagrama
max(ABC$Lt)

# Etiquetas del eje x:
brks <- 0:ceiling(max(ABC$Lt))
brks

# Tipos de línea
line_types <- c("Abundancia"=1,"Biomasa"= 3)

## Diagrama

str(ABC)
max(brks)
round(exp(0:max(brks)))



(ABC.plot <- ggplot(data = ABC) + 
  # Separamos los gráficos por sitio
  facet_wrap(~Sitio, ncol = 3) +
  # Colocamos las líneas
  geom_line(aes(x = Lt, y = Accum.Biomass, linetype= "Biomasa"), lwd = 0.5) +
  geom_line(aes(x = Lt, y = Accum.Abund, linetype = "Abundancia"), lwd = 0.5) +
  #Sombreamos el espacio entre curvas
  geom_ribbon(aes( ymin = Accum.Abund, ymax = Accum.Biomass, x = Lt), fill = "grey12", alpha = 0.1) +
  # Leyenda
  scale_linetype_manual(name = " ", values=line_types) +
  # Personalizacion
  theme_ipsum() +
  xlab(expression('Rango de especies'~(Log[e]~Scale))) +
  ylab("Dominancia acumulada %") +
  scale_x_continuous(breaks = brks, 
                     labels = round(exp(0:max(brks))))+
  scale_y_continuous(limits = c(0,101), 
                     breaks = seq(0,100,10)) +
  geom_text(data = Ws, aes(x = max(brks), y = 0, label = paste('W =',
                                              W, 
                                              sep=' ')),
            vjust = "inward", hjust = "inward",
            size = 2)  +
    theme(axis.text.x = element_text(angle =90, size = 8),
                             axis.text.y = element_text(size = 8))  +
    theme(
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8),
      plot.title = element_text(size = 14)
    ) +
    theme(legend.position="bottom")
)

### Guardar gráfico como pdf

extrafont::loadfonts()

pdf("~/Ecoanalysis/Figs/ABC.plot.pdf", width = 15/2.54, height = 20/2.54)
print(ABC.plot)
dev.off()
#########################################


######## Gráfico de los valores de W ###########

str(Ws)

(W.plot <- ggplot(Ws,aes(x = Sitio, y = W,
                    colour = Sitio, group = 1)) +
  geom_line(stat = "identity", color = "grey") +
  geom_point(shape=21, color="black", fill="#69b3a2", size=6) +
  geom_hline(yintercept = 0, linetype = "dotted", lwd =1.5) +
  ylab("W") +
  xlab("Sitios") +
  theme_ipsum() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4)) 
)

### Guardar gráfico como pdf

extrafont::loadfonts()

pdf("~/Ecoanalysis/Figs/W.plot.pdf", width = 15/2.54, height = 9/2.54)
print(W.plot)
dev.off()

