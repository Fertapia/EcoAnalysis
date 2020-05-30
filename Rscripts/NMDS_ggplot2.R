

########## Ordenamiento No Métrico Multi-dimensional NMDS - personalizado con ggplot2 ###########

# Datos

## Para este ejemplo se usan los datos de :


# dune del paquete vegan
library(vegan)
df <- data("dune") # abundancias
df_env <- data("dune.env") # ambiental
df <- dune
df_env <- dune.env

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


################ 1. NMDS básico sitios (unidades de estudio) + etiquetas de sitios ################################

ord <- metaMDS(df)

# Extraemos los scores
df_ord <- as.data.frame(scores(ord, display = "sites", scaling = 1, choices = c(1,2)))

# Elegimos etiquetas de sitios y ejes
axis.labels <- ord_labels(ord)[choices]
xlab <- axis.labels[1]
ylab <- axis.labels[2]
colnames(df_ord) <- c("x", "y")
df_ord$site <- rownames(df_ord)

library(ggplot2)
# Diagrama
(
  ggplot(data = df_ord, aes(x = x, y = y)) +
  geom_point(size = 3)+
  xlab(xlab) +
  ylab(ylab) +
  geom_text(mapping = aes(label = site, x = x, y = y),
            size = 2, vjust=-1) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw()
)

# Diagrama personalizado https://www.r-graph-gallery.com/273-custom-your-scatterplot-ggplot2.html
library(hrbrthemes)

(
  ggplot(data = df_ord, aes(x = x, y = y)) +
    geom_point(
      color = "black", 
      fill = "#69b9a2",
      shape = 21,
      alpha = 0.5, 
      size = 6,
      stroke = 2) +
    xlab(xlab) +
    ylab(ylab) +
    geom_text(mapping = aes(label = site, x = x, y = y),
              size = 2) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_ipsum()
)

############## 2. NMDS básico + agrupaciones (1 factor - Color/forma) ########################

df_ord
var.temp <- df_env$Management

# Añadimos a df_ord la variable Management de df_env

df_ord$Management <- var.temp

# Diagrama
(
  ggplot(data = df_ord, aes(x = x, y = y, colour = Management, shape = Management)) +
    geom_point(size = 3) +
    xlab(xlab) +
    ylab(ylab) +
    geom_text(mapping = aes(label = site, x = x, y = y),
              size = 2, vjust=-1, color = "black") +
    scale_shape_manual(values = c(15 ,17, 16, 18)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
)

# Diagrama personalizado https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html
library(viridis)

(
  ggplot() +
    geom_path(data = df_ord, aes(x = x, y = y), color = "grey") +
    geom_point(data = df_ord, aes(x = x, y = y, fill = Management),
               color = "black",
               shape = 21,
               alpha = 0.5, 
               size = 6,
               stroke = 2) +
    xlab(xlab) +
    ylab(ylab) +
    geom_text(data = df_ord,mapping = aes(label = site, x = x, y = y),
              size = 2, color = "black") +
    scale_fill_viridis(discrete = TRUE) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_vline(xintercept = 0, lty = 2) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    theme_ipsum() 

    )


# Diagrama personalizado + Segmentos-flechas que unen los puntos en orden, útil cuando las unidades de estudio poseen un orden correlativo

(
  ggplot() +
    geom_point(data = df_ord, aes(x = x, y = y, color = Management),
               alpha = 0.9, 
               size = 2,
               stroke = 2) +
    xlab(xlab) +
    ylab(ylab) +
    geom_text(data = df_ord,mapping = aes(label = site, x = x, y = y),
              size = 3, color = "black", vjust = -1) +
    geom_segment(data = df_ord, 
                 color = "#69b3a2",
                 aes(
                   x = x,
                   y = y,
                   xend = c(tail(x, n = -1), NA),
                   yend = c(tail(y, n = -1), NA)
                 ),
                 arrow = arrow(length = unit(0.3, "cm"))) +
    scale_color_viridis(discrete = TRUE) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_ipsum() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
)

######################### 3. NMDS + Polígonos/elipses/spiders #####################################

# Obtenemos las coordenadas de las unidades de estudio

ord <- metaMDS(df)

########## Extraemos los scores df_ord
df_ord <- as.data.frame(scores(ord, display = "sites", scaling = 1, choices = c(1,2)))

# Elegimos etiquetas de sitios y ejes
axis.labels <- ord_labels(ord)[choices]
xlab <- axis.labels[1]
ylab <- axis.labels[2]
colnames(df_ord) <- c("x", "y")

site <- rownames(df_ord)

axis.labels <- ord_labels(ord)[choices]
xlab <- axis.labels[1]
ylab <- axis.labels[2]

# Añadimos a df_ord la variable Management de df_env
groups <- as.factor(df_env$Management)
Group <- groups
df_ord$Group <- Group



# Definimos los grupos (niveles del factor)
show.groups <- as.vector(levels(groups))

########## Elipses

# Obtenemos el centro de los elipses df_mean.ord
df_mean.ord <- aggregate(df_ord[,1:2], by = list(df_ord$Group), mean)
colnames(df_mean.ord) <- c("Group", "x", "y")
df_mean.ord <- df_mean.ord[df_mean.ord$Group %in% show.groups, ]

# Elegir valores para los elipses:
# sd para desviaciones estándar de los puntos
# se para errores estándar
# ehull para elipses que encierran todos los puntos del grupo

kind <- "se" 

# En caso de elegir se, elegir el intervalo de confianza

conf <- 0.95

# Calculamos los elipses en df_ellipse

rslt <- ordiellipse(ord, groups=groups, display = "sites", scaling=scaling, 
                           choices=choices, kind = kind, show.groups = show.groups,
                           draw = "none", conf = conf, label = label)


df_ellipse <- data.frame()
for(g in show.groups) {
  df_ellipse <- rbind(df_ellipse, 
                      cbind(as.data.frame(with(df_ord[df_ord$Group==g,],
                                                           vegan:::veganCovEllipse(rslt[[g]]$cov,rslt[[g]]$center, 
                                                                                   rslt[[g]]$scale))),Group=g))
}
colnames(df_ellipse) <- c("x", "y", "Group")
df_ellipse <- df_ellipse[ , c(3,1,2)]


# Diagrama
(
ggplot() +
  geom_point(data=df_ord, aes(x=x, y=y, color=Group), size = 3) +
  xlab(xlab) + 
  ylab(ylab) +
  geom_path(data = df_ellipse, aes(x=x, y=y, color=Group), show.legend = FALSE) +
  geom_text(data=df_mean.ord, aes(x=x, y=y, label=Group, color=Group), show.legend = FALSE) +
  coord_fixed()
)



########### Polígonos

# Calculamos los polígonos df_hull
rslt.hull <- ordihull(ord, groups = groups, scaling = scaling, choices = choices, 
                             show.groups = show.groups, draw = "none")

df_hull <- data.frame()
df_temp <- data.frame()
for (g in show.groups) {
  x <- rslt.hull[[g]][ , 1]
  y <- rslt.hull[[g]][ , 2]
  Group <- rep(g, length(x))
  df_temp <- data.frame(Group = Group, x=x, y=y)
  df_hull <- rbind(df_hull, df_temp)
}

# Diagrama
(
  ggplot() +
    geom_point(data=df_ord, aes(x=x, y=y, color=Group), size = 3) +
    xlab(xlab) + 
    ylab(ylab) +
    geom_path(data=df_hull, aes(x=x, y=y, color=Group), show.legend = FALSE) +
    geom_text(data=df_mean.ord, aes(x=x, y=y, label=Group, color=Group), show.legend = FALSE) +
    coord_fixed()
)


############ Spiders

# Calculamos los spider df_spiders
df_spiders <- df_ord
df_spiders$cntr.x <- NA
df_spiders$cntr.y <- NA
for (g in show.groups) {
  df_spiders[which(df_spiders$Group==g), 4:5] <- df_mean.ord[which(df_mean.ord==g), 2:3]
}

df_spiders <- df_spiders[ , c(3,4,5,1,2)]
df_spiders <- df_spiders[order(df_spiders$Group), ]
df_spiders <- df_spiders[df_spiders$Group %in% show.groups, ]

# Diagrama
(
  ggplot() +
    geom_point(data=df_ord, aes(x=x, y=y, color=Group), size = 3) +
    xlab(xlab) + 
    ylab(ylab) +
    geom_segment(data=df_spiders, aes(x=cntr.x, xend=x, y=cntr.y, yend=y, color=Group), show.legend = FALSE) +
    geom_text(data=df_mean.ord, aes(x=x, y=y, label=Group, color=Group), show.legend = FALSE) +
    coord_fixed()
)


### Diagrama personalizado

## Pueden combinarse 

df_ord$site <- site

p <-  ggplot() +
    # Plotear unidades de estudio 
    geom_point(data = df_ord, aes(x = x, y = y, color = Group),
               alpha = 0.9, 
               size = 2,
               stroke = 2,
               shape = 21)  +

    # Etiquetas de los sitios
    geom_text(data = df_ord,
              mapping = aes(label = site, x = x, y = y),
              size = 3, 
              color = "black", 
              vjust = -1) +
    # Segmentos que unen unidades de estudio
    geom_segment(data = df_ord, 
                 color = "black",
                 lwd = 1,
                 linetype = 3,
                 aes(
                   x = x,
                   y = y,
                   xend = c(tail(x, n = -1), NA),
                   yend = c(tail(y, n = -1), NA)
                 ),
                 arrow = arrow(length = unit(0.3, "cm"))) +
    # Elipses se 95%
    geom_path(data = df_ellipse, aes(x=x, y=y, color=Group), show.legend = FALSE) +
    # Polígonos
    geom_path(data=df_hull, aes(x=x, y=y, color=Group), show.legend = FALSE) +
    # Spiders
    geom_segment(data=df_spiders, aes(x=cntr.x, xend=x, y=cntr.y, yend=y, color=Group), show.legend = FALSE) +
    # Personalización
    xlab(xlab) +
    ylab(ylab) +
    scale_color_viridis(discrete = TRUE, option = "magma", begin = 0.2, end = 0.8) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_ipsum() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    coord_fixed()

  ### Guardar gráfico como pdf
  
  extrafont::loadfonts()
  
  pdf("~/Ecoanalysis/Figs/NMDS_ordi.pdf", width = 15/2.54, height = 12/2.54)
  print(p)
  dev.off()

############## 4. NMDS básico + variables ambientales numéricas ########################

## 4.1 Burbujas ####

# Extraemos los scores
  
df_ord <- as.data.frame(scores(ord, display = "sites", choices = choices))

axis.labels <- ord_labels(ord)[choices]
  
df_ord$var <- df_env$A1

xlab <- axis.labels[1]
ylab <- axis.labels[2]
colnames(df_ord) <- c("x", "y", "var")

# Diagrama básico

ggplot(data = df_ord,
       aes(x = x, y = y, size = var)) +
  geom_point(alpha = 0.5) +
  scale_size(range = c(.1, 12), name="A1") +
  xlab(xlab) +
  ylab(ylab) +
  labs(size = var.label) +
  coord_fixed(ratio = 1)

# Añadir colores 

df_ord$group <- df_env$Management

(
p <- ggplot(data = df_ord,
       aes(x = x, y = y, size = var, fill = group)) +
  geom_point(alpha = 0.5, shape = 21, color = "black") +
  scale_size(range = c(.1, 16), name="A1") +
  scale_fill_viridis(discrete = TRUE, guide = FALSE, option = "A") +
  theme_ipsum() +
  theme(legend.position = "bottom") +
  xlab(xlab) +
  ylab(ylab) +
  labs(size = var.label)  +
  coord_fixed(ratio = 1)
)
  
extrafont::loadfonts()

pdf("~/Ecoanalysis/Figs/NMDS_bubble.pdf", width = 15/2.54, height = 12/2.54)
print(p)
dev.off() 
  

## 4.2 NMDS + diagrama de contorno ####

# Extraemos las coordenadas de los sitios

df_ord <- as.data.frame(scores(ord, choices = choices, display = "sites"))
colnames(df_ord) <- c("x", "y")

# Variable ambiental

env.var <- df_env$A1

# Extraems los datos de ordisurf
ordi <- ordisurf(ord, env.var, plot = FALSE)
ordi.grid <- ordi$grid
ordi.data <- expand.grid(x = ordi.grid$x, y = ordi.grid$y)
ordi.data$z <- as.vector(ordi.grid$z)
df_surf <- data.frame(na.omit(ordi.data))

# Etiquetas de ejes
axis.labels <- ord_labels(ord)[choices]
xlab <- axis.labels[1]
ylab <- axis.labels[2]

# Calcular ancho de rango
r <- range(env.var)
binwidth <- (r[2] - r[1])/15

# Diagrama básico

ggplot(data=df_ord, aes(x=x, y=y)) + 
  geom_point(size = pt.size) +
  xlab(xlab) + 
  ylab(ylab) +
  stat_contour(data = df_surf, aes(x=x, y=y, z=z, color= ..level..), binwidth=binwidth) +
  labs(color=var.label) + 
  coord_fixed(ratio=1) 

# Diagrama personalizado
library(ggnewscale)
library(metR)

site <- rownames(df_ord)
df_ord$site <- site
df_ord$group <- df_env$Management


p <- ggplot() + 
  geom_point(data = df_ord, aes(x = x, y =y, color = group, shape = group), size = 3) +
  geom_text(data = df_ord, mapping= aes(label = site, x = x, y = y), size =3, vjust = -1) +
  scale_shape_manual("Management", values = c(15,16,17,18)) +
  scale_color_grey("Management") +

  new_scale_color() +
  
  stat_contour(data = df_surf, aes(x=x, y=y, z=z, color= ..level..), binwidth=binwidth) +
  geom_text_contour(data = df_surf, aes(x = x, y = y, z = z), stroke = 0.2, binwidth=binwidth) +
  scale_color_viridis("A1", option = "magma", begin = 0.2, end = 0.8) +
  labs(color=var.label) + 
  coord_fixed(ratio=1) +
  theme_ipsum() +
  #theme(legend.position = "none") +
  xlab(xlab) + 
  ylab(ylab) 


extrafont::loadfonts()

pdf("~/Ecoanalysis/Figs/NMDS_contour.pdf", width = 15/2.54, height = 12/2.54)
print(p)
dev.off() 

### 4.3 Ajuste variables ambientales ####

# Extaer scores
df_ord <- scores(ord, display = "sites", choices = choices, scaling = 1)
df_ord <- as.data.frame(df_ord)

axis.labels <- ord_labels(ord)[choices]
colnames(df_ord) <- c("x", "y")

# Extraer ajustes df_arrows
fit <- envfit(ord, df_env, choices = choices, perm = 999)

df_arrows <- as.data.frame(scores(fit, "vectors"))

# Multiplicador para ajustar las flechas al diagrama
mult <- scale_arrow(df_arrows, df_ord[, c("x", "y")])

df_arrows <- mult*df_arrows
df_arrows$var <- rownames(df_arrows)
df_arrows$p.val <- fit$vectors$pvals

colnames(df_arrows) <- c("x", "y", "var", "p.val")
df_arrows <- df_arrows[df_arrows$pval <- alpha,]

df_arrows

xlab <- axis.labels[1]
ylab <- axis.labels[2]


#Diagrama básico
ggplot(data=df_ord, aes(x=x, y=y)) + 
  geom_point(size = 2) +
  xlab(xlab) + 
  ylab(ylab) +
  geom_segment(data = df_arrows, 
               aes(x = 0, xend = x, 
                   y = 0, yend = y),
               arrow = arrow(angle = angle, 
                             length = unit(len, unit)),
               color = arrow.col) +
  geom_text(data=df_arrows, 
            aes(x=x, y=y, label=var), 
            color=arrow.col, 
            hjust="outward") +
  coord_fixed(ratio = 1)

# Diagrama
ggplot(data=df_ord, aes(x=x, y=y)) + 
  geom_point(size=pt.size) +
  xlab(xlab) + 
  ylab(ylab) +
  geom_segment(data = df_arrows, 
               aes(x = 0, xend = x, 
                   y = 0, yend = y),
               arrow = arrow(angle = angle, 
                             length = unit(len, unit)),
               color = "blue") +
  geom_text(data=df_arrows, 
            aes(x=x, y=y, label=var), 
            color= "blue", 
            hjust="outward") +
  coord_fixed(ratio = 1)


### 4.4 Contorno + flechas ####

# Combinar el diagrama de contorno con las flecha de interés

p <- ggplot() + 
  geom_point(data = df_ord, aes(x = x, y =y, color = group, shape = group), size = 3) +
  geom_text(data = df_ord, mapping= aes(label = site, x = x, y = y), size =3, vjust = -1) +
  scale_shape_manual("Management", values = c(15,16,17,18)) +
  scale_color_grey("Management") +
  
  new_scale_color() +
  
  stat_contour(data = df_surf, aes(x=x, y=y, z=z, color= ..level..), binwidth=binwidth) +
  geom_text_contour(data = df_surf, aes(x = x, y = y, z = z), stroke = 0.2, binwidth=binwidth) +
  scale_color_viridis("A1", option = "magma", begin = 0.2, end = 0.8) +
  labs(color=var.label) + 
  coord_fixed(ratio=1) +
  
  geom_segment(data = df_arrows, 
               aes(x = 0, xend = x, 
                   y = 0, yend = y),
               arrow = arrow(angle = angle, 
                             length = unit(len, unit)),
               color = "blue") +
  
  geom_text(data=df_arrows, 
            aes(x=x, y=y, label=var), 
            color= "blue", 
            hjust="outward") +
  
  theme_ipsum() +
  #theme(legend.position = "none") +
  xlab(xlab) + 
  ylab(ylab) 

p

extrafont::loadfonts()

pdf("~/Ecoanalysis/Figs/NMDS_contour_fit.pdf", width = 15/2.54, height = 12/2.54)
print(p)
dev.off() 
