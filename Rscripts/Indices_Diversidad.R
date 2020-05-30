
########## Índices ######################## 
N0 <- specnumber(df) # Riqueza - número de especies
Ab <- rowSums(df) # Abundancia por muestra - número de individuos
H <- diversity(df) # Índice de Shannon- Wiener base e
Hb2 <- diversity(df, base = 2) # Índice de Shannon-Wiener base 2
N1 <- exp(H) # Números de Hill - Hill 1 (N1) - Diversidad verdadera - Exponencial de Shannon base e
N1b2 <- 2^Hb2 # Hill 1 calculado desde la base 2
N2 <- diversity(df, "inv") # Números de hill (N2) - Número de especies dominantes
J <- H/log(N0) # Equidad de Pielou
E10 <- N1/N0 # Equidad de Shannon
E20 <- N2/N0 # Equidad de Simpson
D <- diversity(df, index = "simpson") # Diversidad de Simpson

# Tabla de índices
alpha_div <- cbind(df_env, N0, Ab, H, Hb2, N1, N1b2, N2, J, E10, E20, D)

################################################


############ Plots de diversidad #######################

# Calculamos intervalos de confianza

alpha.ci <- alpha_div %>% 
  group_by(Sitio) %>% 
  summarise(
            prom.N0 = mean(N0),
            se.N0 = sd(N0)/sqrt(n()),
            prom.H = mean(H),
            se.H = sd(H)/sqrt(n()),
            prom.D = mean(D),
            se.D = sd(D)/sqrt(n())
            )

# Diagrama 

N0.plot <- 
  ggplot(alpha.ci,aes(x = Sitio, y = prom.N0,
                              colour = Sitio)) +
  geom_point(size = 3) +
  scale_color_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  geom_errorbar(aes(ymin = prom.N0 - se.N0,
                    ymax = prom.N0 + se.N0), 
                width = 0.1) +
  ylab("Riqueza de especies") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))



H.plot <- 
    ggplot(alpha.ci,aes(x = Sitio, y = prom.H,
                        colour = Sitio)) +
    geom_point(size = 3) +
    scale_color_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
    geom_errorbar(aes(ymin = prom.H - se.H,
                      ymax = prom.H + se.H), 
                  width = 0.1) +
    ylab("Shannon H'") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))


D.plot <- 
    ggplot(alpha.ci,aes(x = Sitio, y = prom.D,
                        colour = Sitio)) +
    geom_point(size = 3) +
    scale_color_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
    geom_errorbar(aes(ymin = prom.D - se.D,
                      ymax = prom.D + se.D), 
                  width = 0.1) +
    ylab("Simpson 1-D") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))


legend <- get_legend(D.plot)

(Div.plot <- plot_grid (N0.plot + theme (legend.position = "none"),
           H.plot +  theme (legend.position = "none"),
           D.plot + theme (legend.position = "none"), 
           ncol = 3))

# PDF 15 X 6 cm
pdf("~/Ecoanalysis/Figs/Div.plot_.pdf", width = 15/2.54, height = 6/2.54)
print(Div.plot)
dev.off()

####################################################

########################## Análisis de varianza (ANOVA) ####################

Div.var <- lm(H ~ Sitio, data = alpha_div)

autoplot(Div.var, smoot.colour = NA) +
  theme_bw()

Div.av <- aov(Div.var)
summary(Div.av)

tukey.test <- TukeyHSD(Div.av)
tukey.test

###########################################