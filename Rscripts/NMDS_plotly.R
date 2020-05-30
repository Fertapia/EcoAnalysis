

# Plotly básico

# Libraries
library(ggplot2)
library(dplyr)
library(plotly)
library(viridis)
library(hrbrthemes)
library(gapminder)

# data en gapminder

data <- gapminder %>% 
  filter(year == "2007") %>% 
  dplyr::select(-year)

# Interactive version
df <- data %>%
  mutate(gdpPercap=round(gdpPercap,0)) %>% 
  mutate(pop=round(pop/1000000,2)) %>%
  mutate(lifeExp=round(lifeExp,1)) %>% 
  
  # Reorder countries to having big bubbles on top
  arrange(desc(pop)) %>%
  mutate(country = factor(country, country)) %>% 
  
  # prepare text for tooltip
  mutate(text = paste("Country: ", country, "\nPopulation (M): ", pop, "\nLife Expectancy: ", 
                      lifeExp, "\nGdp per capita: ", gdpPercap, sep=""))
  
df

p <- ggplot(data = df, aes(x = gdpPercap, y = lifeExp, size = pop, color = continent, text = text)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(1.4, 19), name = "Population (M)") +
  scale_color_viridis(discrete = TRUE, guide = FALSE) +
  theme_ipsum() +
  theme(legend.position = "none")

pp <- ggplotly(p, tooltip = "text")
pp


###### 1. NMDS - plotly

# dune del paquete vegan
library(vegan)
df <- data("dune") # abundancias
df_env <- data("dune.env") # ambiental
df <- dune
df_env <- dune.env

# Extraemos los scores
ord <- metaMDS(df)
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


# Complementamos la base de datos

str(df_env)
df_ord$Moisture <- df_env$Moisture
df_ord$Use <- df_ord$Use

str(df_ord)

# Etiquetas

df_ord <- df_ord %>% 
  mutate(text = paste("A1: ", var,
                      "\nManagement: ", group,
                      "\nMoisture: ", Moisture, 
                      sep = ""))


p <- ggplot(data = df_ord,
            aes(x = x, y = y, size = var, fill = group, text = text)) +
  geom_point(alpha = 0.5, shape = 21, color = "black") +
  scale_size(range = c(.1, 16), name="A1") +
  scale_fill_viridis(discrete = TRUE, guide = FALSE, option = "A") +
  theme_ipsum() +
  theme(legend.position = "bottom") +
  xlab(xlab) +
  ylab(ylab) +
  labs(size = var.label)  +
  coord_fixed(ratio = 1)
p

pp <- ggplotly(p, tooltip = "text")
pp









