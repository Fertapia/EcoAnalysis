
### FUENTE: 

#Vegan: an introduction to ordination

library (vegan)

data(dune)

ord <- metaMDS(dune)

ord
class(ord)

plot(ord, type ="n")
points(ord, display = "sites", cex = 0.9, pch = 21, col = "red", bg = "yellow")
text(ord, display = "spec", cex = 0.7, col = "blue")

data(dune.env)
attach(dune.env)
?attach

plot(ord, disp = "sites", type = "n")
ordihull(ord, Management, col = 1:4, lwd =3)
ordiellipse(ord, Management, col = 1:4, kind = "ehull", lwd = 3)
ordiellipse(ord, Management, col = 1:4, draw = "polygon")
ordispider(ord, Management, col = 1:4, label = TRUE)
points(ord, disp = "sites", pch = 21, col = "red", bg = "yellow", cex = 1.3)

ord.fit <- envfit(ord ~A1 + Management, data = dune.env, perm = 999)
ord.fit

plot(ord, dis = "site")
plot(ord.fit)
ordisurf(ord, A1, add = TRUE)


### FUENTE
## Package 'vegan'

# Función metaMDS

data(dune)

sol <- metaMDS(dune)
sol

plot(sol)

sol2 <- metaMDS(dune, previous.best = sol)
sol2

## NMDS Local y stress 2 de monoMDS

sol3 <- metaMDS(dune, model = "hybrid", stress = 2)
sol3

## Cambiar la distancia Bray por Arrhenius 'z' como una medida binaria

sol4 <- metaMDS(dune, distfun = betadiver, distance = "z")
sol4

#################### Función monoMDS #####################
stressplot(sol3)

data(dune)
dis <- vegdist(dune)
m <- monoMDS(dis, model = "loc")
m
plot(m)

###################### Función ordiarrows#####################

example(pyrifos)
mod <- rda(pyrifos)

plot(mod, type = "n")

## Annual succession by ditches, colour by dose
ordiarrows(mod, ditch, label = TRUE, col = as.numeric(dose))
legend("topright", levels(dose), lty=1, col=1:5, title="Dose")

plot(mod, type = "n")
ordiarrows(mod, ditch, label = TRUE,
           show.groups = c("2", "3", "5", "11"))
Arrow <- ordiarrows(mod, ditch, label = TRUE, show = c("6"),
           col = 2)

legend("topright", c("Control", "Pyrifos 44"), lty = 1, col = c(1,2))

ordigrid(mod, ditch, col = as.numeric(dose))

plot(mod)

###################### Función ordiArrowTextXY #####################

## Scale arrows by hand to fill 80% of the plot
## Biplot arrows by hand

data(varespec, varechem)
ord <- cca(varespec ~ Al + P + K, varechem)
plot(ord, display = c("species","sites"))

## biplot scores

bip <- scores(ord, choices = 1:2, display = "bp")

## scaling factor for arrows to fill 80% of plot

(mul <- ordiArrowMul(bip, fill = 0.4))

bip.scl <- bip * mul # Scale the biplot scores
labs <- rownames(bip) # Arrow labels

## calculate coordinate of labels for arrows

(bip.lab <- ordiArrowTextXY(bip.scl, rescale = FALSE, labels = labs))
## draw arrows and text labels

arrows(0, 0, bip.scl[,1], bip.scl[,2], length = 0.1)
text(bip.lab, labels = labs)

## Handling of ordination objects directly
mul2 <- ordiArrowMul(ord, display = "bp", fill = 0.4)
stopifnot(all.equal(mul, mul2))

##################### Función ordihull #################

data(dune)
data(dune.env)
mod <- cca(dune ~ Management, dune.env)
plot(mod, type="n", scaling = "symmetric")

## Catch the invisible result of ordihull...
pl <- with(dune.env, ordihull(mod, Management,
                              scaling = "symmetric", label = TRUE))

## ... and find centres and areas of the hulls

summary(pl)

## use more colours and add ellipsoid hulls

plot(mod, type = "n")
pl <- with(dune.env, ordihull(mod, Management,
                              scaling = "symmetric", col = 1:4,
                              draw="polygon", label =TRUE))
with(dune.env, ordiellipse(mod, Management, scaling = "symmetric",
                           kind = "ehull", col = 1:4, lwd=3))

## ordispider to connect WA and LC scores

plot(mod, dis=c("wa","lc"), type="p")
ordispider(mod)

## Other types of plots
plot(mod, type = "p", display="sites")
cl <- hclust(vegdist(dune))
plot(cl)
ordicluster(mod, cl, prune=3, col = cutree(cl, 4))

## confidence ellipse: location of the class centroids
plot(mod, type="n", display = "sites")
with(dune.env, text(mod, display="sites", labels = as.character(Management),
                    col=as.numeric(Management)))
pl <- with(dune.env, ordiellipse(mod, Management, kind="se", conf=0.95, lwd=2,
                                 draw = "polygon", col=1:4, border=1:4,
                                 alpha=63))
summary(pl)


## add confidence bars
with(dune.env, ordibar(mod, Management, kind="se", conf=0.95, lwd=2, col=1:4,
                       label=TRUE))

############ Función Ordilabel ###################

data(dune)
ord <- cca(dune)
plot(ord, type = "n")
ordilabel(ord, dis="sites", cex=1.2, font=3, fill="hotpink", col="blue")
## You may prefer separate plots, but here species as well
ordilabel(ord, dis="sp", font=2, priority=colSums(dune))

############## Función Ordiplot ###############

# Draw a plot for a non-vegan ordination (cmdscale).
data(dune)
dune.dis <- vegdist(wisconsin(dune))
dune.mds <- cmdscale(dune.dis, eig = TRUE)
dune.mds$species <- wascores(dune.mds$points, dune, expand = TRUE)
pl <- ordiplot(dune.mds, type = "none")
points(pl, "sites", pch=21, col="red", bg="yellow")
text(pl, "species", col="blue", cex=0.9)
# Default plot of the previous using identify to label selected points
## Not run:
pl <- ordiplot(dune.mds)
identify(pl, "spec")
## End(Not run)

############## Función Ordipointlabel ###############

data(dune)
ord <- cca(dune)
plt <- ordipointlabel(ord)
plt
## set scaling - should be no warnings!

ordipointlabel(ord, scaling = "sites")
## plot then add

plot(ord, scaling = "symmetric", type = "n")
ordipointlabel(ord, display = "species", scaling = "symm", add = TRUE)
ordipointlabel(ord, display = "sites", scaling = "symm", add = TRUE)

## redraw plot without rerunning SANN optimisation

plot(plt)


############ Función ordiresids ###################
data(varespec)
data(varechem)
mod <- cca(varespec ~ Al + P + K, varechem)
ordiresids(mod)
ordiresids(mod, formula = Residuals ~ Fitted | Species, residuals="standard",
           cex = 0.5)


########### Función oristep ###################

## See add1.cca for another example
### Dune data
data(dune)
data(dune.env)
mod0 <- rda(dune ~ 1, dune.env) # Model with intercept only
mod1 <- rda(dune ~ ., dune.env) # Model with all explanatory variables
## With scope present, the default direction is "both"
mod <- ordistep(mod0, scope = formula(mod1))
mod
## summary table of steps
mod$anova
## Example of ordistep, forward
ordistep(mod0, scope = formula(mod1), direction="forward")
## Example of ordiR2step (always forward)
## stops because R2 of 'mod1' exceeded
ordiR2step(mod0, mod1)

########### Función orisurf ###################

data(varespec)
data(varechem)
vare.dist <- vegdist(varespec)
vare.mds <- monoMDS(vare.dist)
ordisurf(vare.mds ~ Baresoil, varechem, bubble = 5)

## as above but without the extra penalties on smooth terms,
## and using GCV smoothness selection (old behaviour of `ordisurf()`):

ordisurf(vare.mds ~ Baresoil, varechem, col = "blue", add = TRUE,
         select = FALSE, method = "GCV.Cp")

## Cover of Cladina arbuscula

fit <- ordisurf(vare.mds ~ Cladarbu, varespec, family=quasipoisson)

## Get fitted values

calibrate(fit)

## Variable selection via additional shrinkage penalties
## This allows non-significant smooths to be selected out
## of the model not just to a linear surface. There are 2
## options available:
## - option 1: `select = TRUE` --- the *default*


ordisurf(vare.mds ~ Baresoil, varechem, method = "REML", select = TRUE)
## - option 2: use a basis with shrinkage

ordisurf(vare.mds ~ Baresoil, varechem, method = "REML", bs = "ts")
## or bs = "cs" with `isotropic = FALSE`
## Plot method

plot(fit, what = "contour")
## Plotting the "gam" object

plot(fit, what = "gam") ## 'col' and 'cex' not passed on
## or via plot.gam directly

library(mgcv)
plot.gam(fit, cex = 2, pch = 1, col = "blue")

## 'col' effects all objects drawn...
### controlling the basis functions used
## Use Duchon splines

ordisurf(vare.mds ~ Baresoil, varechem, bs = "ds")

## A fixed degrees of freedom smooth, must use 'select = FALSE'

ordisurf(vare.mds ~ Baresoil, varechem, knots = 4,
         fx = TRUE, select = FALSE)

## An anisotropic smoother with cubic regression spline bases

ordisurf(vare.mds ~ Baresoil, varechem, isotropic = FALSE,
         bs = "cr", knots = 4)

## An anisotropic smoother with cubic regression spline with
## shrinkage bases & different degrees of freedom in each dimension

ordisurf(vare.mds ~ Baresoil, varechem, isotropic = FALSE,
         bs = "cs", knots = c(3,4), fx = TRUE,
         select = FALSE)


############ Función orditkplot ####################

## The example needs user interaction and is not executed directly.
## It should work when pasted to the window.
## Not run:
data(varespec)
ord <- cca(varespec)
## Do something with the graph and end by clicking "Dismiss"
orditkplot(ord, mar = c(4,4,1,1)+.1, font=3)
## Use ordipointlabel to produce a plot that has both species and site
## scores in different colors and plotting symbols
pl <- ordipointlabel(ord)
orditkplot(pl)

########## FUnción oridtorp ###################

## A cluttered ordination plot :
data(BCI)
mod <- cca(BCI)
plot(mod, dis="sp", type="t")
# Now with orditorp and abbreviated species names
cnam <- make.cepnames(names(BCI))
plot(mod, dis="sp", type="n")
stems <- colSums(BCI)
orditorp(mod, "sp", label = cnam, priority=stems, pch="+", pcol="grey")
## show select in action
set.seed(1)
take <- sample(ncol(BCI), 50)
plot(mod, dis="sp", type="n")
stems <- colSums(BCI)
orditorp(mod, "sp", label = cnam, priority=stems, select = take,
         pch="+", pcol="grey")


######### Función ordixyplot ##############

data(dune, dune.env)
ord <- cca(dune)
## Pairs plots
ordisplom(ord)
ordisplom(ord, data=dune.env, choices=1:2)
ordisplom(ord, data=dune.env, form = ~ . | Management, groups=Manure)
## Scatter plot with polygons
ordixyplot(ord, data=dune.env, form = CA1 ~ CA2 | Management,
           groups=Manure, type = c("p","polygon"))
## Choose a different scaling
ordixyplot(ord, scaling = "symmetric")
## ... Slices of third axis
ordixyplot(ord, form = CA1 ~ CA2 | equal.count(CA3, 4),
           type = c("g","p", "polygon"))
## Display environmental variables
ordixyplot(ord, envfit = envfit(ord ~ Management + A1, dune.env, choices=1:3))
## 3D Scatter plots
ordicloud(ord, form = CA2 ~ CA3*CA1, groups = Manure, data = dune.env)
ordicloud(ord, form = CA2 ~ CA3*CA1 | Management, groups = Manure,
          data = dune.env, auto.key = TRUE, type = c("p","h"))


##################################################################################################################################
############################# PAQUETE: ggordiplots ###############################################################################
##################################################################################################################################

install.packages("ggordiplots")

remotes::install_github("jfq3/ggordiplots")

library("ggordiplots")
library("ggplot2")

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

################## Función gg_envfit ############################
library(ggordiplots)

data("varespec")
data("varechem")
vare.dist <- vegdist(varespec)
vare.mds <- monoMDS(vare.dist)
Envfit <- gg_envfit(ord=vare.mds, env=varechem)


### Analizando la función:

groups = NA
scaling = 1
choices = c(1,2)
perm = 999
alpha = 0.05
angle = 20
len = 0.5
unit = "cm"
arrow.col = "red"
pt.size = 3
plot = TRUE

scores(vare.mds)


df_ord <- vegan::scores(vare.mds, display = "sites", choices = c(1,2), scaling= 1)

df_ord <- as.data.frame(df_ord)

axis.labels <- ord_labels(vare.mds)[choices]

colnames(df_ord) <- c("x", "y")


fit <- vegan::envfit(vare.mds, varechem, choices = choices, perm = perm)

df_arrows <- as.data.frame(scores(fit, "vectors"))

mult <- scale_arrow(df_arrows, df_ord[, c("x", "y")])

df_arrows <- mult * df_arrows
df_arrows$var <- rownames(df_arrows)
df_arrows$p.val <- fit$vectors$pvals
df_arrows
colnames(df_arrows) <- c("x", "y", "var", "p.val")
df_arrows <- df_arrows[df_arrows$p.val<=alpha, ]

xlab <- axis.labels[1]
ylab <- axis.labels[2]

ggplot(data=df_ord, aes(x=x, y=y)) + 
  geom_point(size=pt.size) +
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
  


################## Función gg_ordibubble ########################

data(dune)
data(dune.env)
dune.bray <- vegdist(dune, method = "bray")
ord <- cmdscale(dune.bray, k=(nrow(dune)-1), eig=TRUE, add=TRUE)
BUBLE <- gg_ordibubble(ord, env.var=dune.env$A1, var.label="A1")


### Analizando la función:

var.label = "Level"
choices = c(1,2)

df_ord <- as.data.frame(vegan::scores(ord, display = "sites", choices = choices))

axis.labels <- ord_labels(ord)[choices]

df_ord$var <- dune.env$A1

xlab <- axis.labels[1]
ylab <- axis.labels[2]
colnames(df_ord) <- c("x", "y", var.label)

ggplot(data = df_ord,
       aes(x = x, y = y, size = dune.env$A1)) +
  geom_point() +
  xlab(xlab) +
  ylab(ylab) +
  labs(size = var.label) +
  coord_fixed(ratio = 1)
 

################ Función gg_ordicluster ########################

data(dune)
data(dune.env)
dune.bray <- vegdist(dune, method="bray")
ord <- metaMDS(dune, k=3)
cluster <- hclust(dune.bray, method="complete")
CLUSTER <-gg_ordicluster(ord, cluster=cl, treatments=dune.env$Management, prune=3, col=cutree(cl, 4))



### Analizando la función:

treatments = dune.env$Management
choices = c(1,2)
prune = 0
col = 1
pt.size = 3
plot = TRUE
display <- "sites"
show.legend <- TRUE
n.trts <- nlevels(as.factor(treatments))


gg_ordicluster <- function (ord, cluster, treatments=NA, choices=c(1,2), prune = 0, col = 1, pt.size = 3, plot=TRUE)
{

if (is.numeric(treatments)) {
  stop("'treatments' cannot be numeric")
}
n.trts <- nlevels(as.factor(treatments))
if (n.trts==0) {
  treatments <- "none"
  n.trts <- 1
}
if (n.trts==1) {
  show.legend <- FALSE
}
else {
  show.legend <- TRUE
}
display <- "sites"
w <- stats::weights(ord, display)
weights.default <- function(object, ...) NULL
w <- eval(w)
mrg <- cluster$merge
ord.scores <- scores(ord, display = display, choices=choices)
if (nrow(mrg) != nrow(ord.scores) - 1)
  stop("Dimensions do not match in 'ord' and 'cluster'")
if ((nrow(ord.scores) != length(treatments)) & (n.trts > 1))
  stop("Dimensions of 'ord' and 'treatments' do not match")
if (length(w) == 1)
  w <- rep(w, nrow(ord.scores))
n <- if (is.null(w))
  rep(1, nrow(ord.scores))
else w
noden <- numeric(nrow(mrg) - prune)
go <- matrix(0, nrow(mrg) - prune, 2)
col <- rep(col, length = nrow(ord.scores))
col <- col2rgb(col)/255
nodecol <- matrix(NA, nrow(mrg) - prune, 3)
for (i in 1:(nrow(mrg) - prune)) {
  a <- mrg[i, 1]
  b <- mrg[i, 2]
  one <- if (a < 0)
    ord.scores[-a, ]
  else go[a, ]
  two <- if (b < 0)
    ord.scores[-b, ]
  else go[b, ]
  n1 <- if (a < 0)
    n[-a]
  else noden[a]
  n2 <- if (b < 0)
    n[-b]
  else noden[b]
  xm <- weighted.mean(c(one[1], two[1]), w = c(n1, n2))
  ym <- weighted.mean(c(one[2], two[2]), w = c(n1, n2))
  go[i, ] <- c(xm, ym)
  noden[i] <- n1 + n2
  colone <- if (a < 0)
    col[, -a]
  else nodecol[a, ]
  coltwo <- if (b < 0)
    col[, -b]
  else nodecol[b, ]
  nodecol[i, ] <- (n1 * colone + n2 * coltwo)/noden[i]
  
  # Rather than plotting the line segments, collect the coordinates and
  # color into a data frame.
  col.a = rgb(t(nodecol[i, ]))
  temp <- c(one[1], one[2], two[1], two[2], col.a)
  
  if (i==1){
    temp2 <- temp
  } else {
    temp2 <- rbind(temp2, temp)
    rownames(temp2) <- NULL # prevents duplicate row names
  }
  
}

colnames(temp2) <- c("x", "y", "xend", "yend", "Group")
temp2 <- as.data.frame(temp2)
j <- sapply(temp2, is.factor)
temp2[j] <- lapply(temp2[j], as.character)
j <- c( rep(TRUE, 4), FALSE)
temp2[j] <- lapply(temp2[j], as.numeric)
df_segments <- temp2

df_ord <- as.data.frame(ord.scores)
axis.labels <- ord_labels(ord)[choices]
df_ord$Treatment <- treatments
colnames(df_ord) <- c("x", "y", "Treatment")

xlab <- axis.labels[1]
ylab <- axis.labels[2]

plt <- ggplot() +
  geom_segment(data=df_segments, aes(x=x, y=y, xend=xend, yend=yend),
               color=df_segments$Group,
               show.legend = FALSE) +
  geom_point(data=df_ord, aes(x=x, y=y, shape=Treatment), size=pt.size,
             show.legend = show.legend) +
  xlab(xlab) + ylab(ylab) + coord_fixed(ratio=1)

if (plot==TRUE) {print(plt)}

invisible(list(df_ord=df_ord, df_segments=df_segments, plot=plt))
}



############## Función gg_ordiplot #############################

data("dune")
data("dune.env")
dune.hel <- decostand(dune, method = "hellinger")
ord <- rda(dune.hel)
gg_ordiplot(ord, groups = dune.env$Management)



### Analizando la función:

ord <- metaMDS(dune)
groups <- dune.env$Management

scaling = 1
choices = c(1,2)
kind = c("sd", "se", "ehull")
conf = NULL
show.groups = "all"
ellipse = TRUE
label = FALSE
hull = FALSE
spiders = FALSE
pt.size = 3


df_ord <- vegan::scores(ord, display = "sites", scaling = scaling, choices)

groups <- as.factor(groups)

show.groups <- as.vector(levels(groups))
show.groups
show.groups[1]=="all"

# Get site coordinates to plot

df_ord <- vegan::scores(ord, display = "sites", scaling=scaling, choices=choices)
axis.labels <- ord_labels(ord)[choices]
df_ord <- data.frame(x=df_ord[ , 1], y=df_ord[ , 2], Group=groups)

# Get ellipse centers to annotate.
df_mean.ord <- aggregate(df_ord[,1:2], by=list(df_ord$Group),mean)
colnames(df_mean.ord) <- c("Group", "x", "y")
df_mean.ord <- df_mean.ord[df_mean.ord$Group %in% show.groups, ]
df_mean.ord

is.null(conf)
rslt <- vegan::ordiellipse(ord, groups=groups, display = "sites", scaling=scaling, 
                           choices=choices, kind = kind, show.groups = show.groups, draw = "none", label = label)


rslt <- vegan::ordiellipse(ord, groups=groups, display = "sites", scaling=scaling, 
                           choices=choices, kind = kind, show.groups = show.groups,
                           draw = "none", conf = 0.95, label = label)



df_ellipse <- data.frame()
for(g in show.groups) {
  df_ellipse <- rbind(df_ellipse, cbind(as.data.frame(with(df_ord[df_ord$Group==g,],
                                                           vegan:::veganCovEllipse(rslt[[g]]$cov,rslt[[g]]$center, 
                                                                                   rslt[[g]]$scale))),Group=g))
}
colnames(df_ellipse) <- c("x", "y", "Group")
df_ellipse <- df_ellipse[ , c(3,1,2)]


rslt.hull <- vegan::ordihull(ord, groups = groups, scaling = scaling, choices = choices, 
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

# Make a data frame for the spiders..
df_spiders <- df_ord
df_spiders$cntr.x <- NA
df_spiders$cntr.y <- NA
for (g in show.groups) {
  df_spiders[which(df_spiders$Group==g), 4:5] <- df_mean.ord[which(df_mean.ord==g), 2:3]
}
df_spiders <- df_spiders[ , c(3,4,5,1,2)]
df_spiders <- df_spiders[order(df_spiders$Group), ]
df_spiders <- df_spiders[df_spiders$Group %in% show.groups, ]

# Make basic ggplot with ellipses.
xlab <- axis.labels[1]
ylab <- axis.labels[2]


ggplot2::ggplot() +
  geom_point(data=df_ord, aes(x=x, y=y, color=Group), size = pt.size) +
  xlab(xlab) + 
  ylab(ylab) +
  # Add ellipses
  geom_path(data = df_ellipse, aes(x=x, y=y, color=Group), show.legend = FALSE) +
  # Add labels
  geom_text(data=df_mean.ord, aes(x=x, y=y, label=Group, color=Group), show.legend = FALSE) +
  # Add hulls
  geom_path(data=df_hull, aes(x=x, y=y, color=Group), show.legend = FALSE) +
  # Add spiders
  geom_segment(data=df_spiders, aes(x=cntr.x, xend=x, y=cntr.y, yend=y, color=Group), show.legend = FALSE) +
  # Detalles
  coord_fixed()


############## Función gg_ordisurf ###############################

data(varespec)
data(varechem)
vare.dist <- vegdist(varespec)
vare.mds <- monoMDS(vare.dist)
gg_ordisurf(vare.mds, env.var = varechem$Baresoil, var.label="Bare Soil")


### Analizando la función:

ord <- metaMDS(varespec)
env.var <- varechem$Baresoil
choices = c(1,2)
var.label = "Level"
pt.size = 3

# Extract ordisurf data for plotting
ordi <- vegan::ordisurf(ord ~ env.var, plot=FALSE) #created the ordisurf object
ordi.grid <- ordi$grid #extracts the ordisurf object
ordi.data <- expand.grid(x = ordi.grid$x, y = ordi.grid$y) #get x and ys
ordi.data$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
df_surf <- data.frame(na.omit(ordi.data)) #gets rid of the nas

# Extract site coordinates for plotting.
df_ord <- as.data.frame(scores(ord, choices = choices, display = "sites"))
colnames(df_ord) <- c("x", "y")

# Make axis labels.
axis.labels <- ord_labels(ord)[choices]
xlab <- axis.labels[1]
ylab <- axis.labels[2]

# Calculate default binwidth
r <- range(env.var)
binwidth <- (r[2]-r[1])/15

## Plotting in ggplot2
ggplot(data=df_ord, aes(x=x, y=y)) + 
  geom_point(size = pt.size) +
  xlab(xlab) + 
  ylab(ylab) +
  stat_contour(data = df_surf, aes(x=x, y=y, z=z, color= ..level..), binwidth=binwidth) +
  labs(color=var.label) + 
  coord_fixed(ratio=1) 






























