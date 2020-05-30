
########### curva de distribución de abundancias de especies - SAD Preston ##################

########### Forma básica ##############

# Modelo de preston para toda la base de datos (df)

pfit <- prestonfit(colSums(df), tiesplit = FALSE)
names(pfit)


pfit$fitted


# La función autoplot diagrama un histograma de la distribución, así como la curva log-normal ajustada con una distribución quasi-Poisson

autoplot(pfit)

############# Diagrama personalizado ##############

# Extraemos los datos octavas y abundancias por rango de clase
df.pr <- data.frame(Octava = as.numeric(names(pfit$freq)), Abundancia = unclass(pfit$freq))
df.pr

df.pr[['OctaveMinusOne']] <- df.pr[['Octava']] - 0.5
df.pr

line.col <- "red"
size <- 1
coefs <- coef(pfit)
brks <- seq(0, nrow(df.pr))
presfun <- function(x, mode, width, S0)  { S0 * exp(-(x - mode)^2/2/width^2) }

ggplot(df.pr, aes_string(x = "Octava", y = "Abundancia")) +
  geom_point() +
  geom_line() +
  stat_function (fun = presfun,
                 args = list(mode = coefs['mode'],
                             width = coefs['width'],
                             S0 = coefs['S0']),
                 colour = line.col, size = size) +
  scale_x_continuous(breaks = brks, labels = 2^brks) +
  theme_classic()



############ Curvas por Estación / Sitio de muestreo #########################

## Creamos un vector que ordene los nombres de las estaciones, usar sólo si el orden de tus estaciones no se acomoda al orden por defecto de R (alfabético)

Sit_level <- c("S1", "S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12")

### Extraemos los datos observados en la lista A

A <- df %>% 
  split(rownames(.)) %>% 
  map(~prestondistr(., tiesplit = FALSE)) %>% 
  map(~(data.frame(Octave = as.numeric(names(.$freq)),
                   Abundance = unclass(.$freq))))

A <- bind_rows(A, .id = "source") # Las compilamos en un data.frame

A$source <- as.factor(A$source)
A$source <- factor(A$source, levels = Sit_level)

levels(A$source)

### Calculamos los datos ajustados a partir de la función presfun 

# presfun <- function(x, mode, width, S0)  { S0 * exp(-(x - mode)^2/2/width^2) }

C <- seq(0, by = 0.1, nrow(df.pr)) # El vector c contiene los valores de la abscisa (X)

B <- df %>% 
  split(rownames(.)) %>% 
  map(~prestondistr(.)) %>% 
  map(~(coefs = coef(.))) %>% 
  map(~(data.frame( Fit = presfun(C, mode = .['mode'], width= .['width'], S0 = .['S0']),
                    Octave=  C)))

B <- bind_rows(B, .id = "source") # Las compilamos en un data.frame

B$source <- as.factor(B$source)
B$source <- factor(B$source, levels = Sit_level)

levels(B$source)
B
###### Graficamos #########

line.col <- "darkred"
size <- 1
brks <- seq(0, nrow(df.pr))

(Preston.plot <- ggplot() +
  geom_line(data = B, aes(x = Octave, y = Fit), colour = line.col, size = size, alpha = 0.5) +
  facet_wrap(~source, ncol = 3) +
  geom_point(data = A, aes(x = Octave, y = Abundance, color = source)) +
  geom_line(data = A, aes(x = Octave, y = Abundance, color = source)) +
  ggtitle("Curvas de Distribución de Abundancia de Especies") +
  scale_color_viridis(discrete = TRUE) +
  scale_x_continuous(breaks = brks, labels = 2^brks) +
  theme_ipsum() +
  xlab("Abundancia geométrica") +
  ylab("Número de especies") +
  theme(axis.text.x = element_text(angle =90, size = 8),
        axis.text.y = element_text(size = 10))  +
  theme(
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8),
    plot.title = element_text(size = 14)
  ) +
  theme(legend.position = "none")
)

# theme_ipsum() es un tema del paquete hrbrthemes y requiere de la instalación de fuentes adicionales para eso el siguiente código (omitir si ya están instalados):

#extrafont::font_import()
#y

extrafont::loadfonts()

pdf("~/Ecoanalysis/Figs/Preston.plot.pdf", width = 15/2.54, height = 20/2.54)
print(Preston.plot)
dev.off()

#################################################################
