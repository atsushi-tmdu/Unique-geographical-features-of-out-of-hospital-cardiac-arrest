library(readxl)
library(dplyr)
library(ggplot2)
library(maptools)
library(INLA)
library(data.table)
library(GISTools)
library(tidyverse)
library(sf)
sf_use_s2()

df<- read_xlsx("../../pre_analysys/all7.xlsx")
jpnmap2<-rgdal::readOGR("./tokyo.shp")

df$longitude= df$longtitude
legionmap <-rgdal::readOGR("../../pre_analysys/tokyoall.shp")

# processing data for analysis
inla_df <- df %>% dplyr::select(id,gender,year,age,byCPR,vf_vt,cpc,time1,time2,time3,time4,witness,latitude,longitude) %>% filter(vf_vt==1)
inla_df <- inla_df %>%  mutate(cpc1_2 =if_else( cpc<3, true=1, false=0))
inla_df <- inla_df %>% mutate(delta_time1=time2-time1, delta_time2=time3-time2,delta_time3=time4-time3)
inla_df <- inla_df %>% mutate(id=row_number())
inla_df<-inla_df %>% dplyr::select(id,gender,year,age,byCPR,cpc1_2,delta_time1,delta_time2,delta_time3,witness,,latitude,longitude)


# make legion map and spatial points dataframe
coords<- inla_df %>%  dplyr::select(longitude,latitude)
coords$admYM<-NULL

pj<-CRS("+proj=longlat +ellps=WGS84 +no_defs")
jpnmap2<-spTransform(jpnmap2,CRS=pj)

colnames(coords) <- c("x","y")
coordinates(coords) <- ~x+y

sp <-SpatialPoints(coords, proj4string = pj)
legionmap<-spTransform(legionmap,CRS=pj)

# plot OHCA cases on the map
g<-ggplot(jpnmap2)+geom_polygon(fill="white", color="black", aes(x=long, y=lat, group=group))+
  geom_point(data=df,aes(x=longitude, y=latitude), size=0.7, color="red", alpha=0.4)+
  annotate('text',x=139.7642509,y=35.701735,label='H',size=4)
ggsave(file="../../supportive_file/results/plot_patientsmap.pdf",dpi=800)
g

# construction of mesh
max.edge = 0.02
bound.outer = 0.02
mesh = inla.mesh.2d(boundary = jpnmap2,
                    loc = sp,
                    max.edge = c(1,3)*max.edge,
                    cutoff = 0.005,
                    offset = c(max.edge, bound.outer))
plot(mesh, asp=1,main ="Our Mesh", lwd=0.2); points(coords, col="red",cex=0.2)
pdf("mesh.pdf")
plot(mesh, asp=1,main ="Our Mesh", lwd=0.2); points(coords, col="red",cex=0.05,add=TRUE)
dev.off()

#spde
spde <- inla.spde2.matern(mesh=mesh, alpha=1.0)

# make stack_est
coordinates_matrix <- as.matrix(inla_df %>%  dplyr::select(latitude,longitude))

A_est <- inla.spde.make.A(mesh=mesh,
                          loc=coordinates_matrix)

s_index <- inla.spde.make.index(name="spatial.field",
                                n.spde = spde$n.spde)

stack_est<- inla.stack(data=list(rate = inla_df$cpc1_2),
                       A = A_est,
                       effects = c(s_index,list(Intercept=1)),
                       tag = "est")

# make grids
cs <- c(3.28084, 3.28084)*0.0001
grdpts<-makegrid(jpnmap2,  cellsize = cs)
plot(jpnmap2)
spgrd <-SpatialPoints(grdpts, proj4string = pj)
spgrdWithin <- SpatialPixels(spgrd[jpnmap2,])
tempgrd<-spgrdWithin
plot(tempgrd<-tempgrd[tempgrd$x1<140],add=TRUE)

# make stack prediction
pred_coords<-tempgrd@coords
pred_coords<-as.data.frame(pred_coords)
n_pred_coords<-nrow(pred_coords)
pred_coords$pred_id<-1:nrow(pred_coords)

A_pred<-inla.spde.make.A(mesh=mesh,
                         loc=as.matrix(pred_coords))

stack_pred <- inla.stack(data=list(rate=NA),
                         A= A_pred,
                         effects = c(s_index,list(Intercept=1)),
                         tag="pred"
                         )
# regression, stack
formula <- rate ~ -1 + Intercept+ age + gender + f(spatial.field,model=spde, group=spatial.field.group)
stack <- inla.stack(stack_est, stack_pred)

### execution  ###
output <- inla(formula, data=inla.stack.data(stack,spde=spde),family='binomial'
               ,control.predictor=list(A=inla.stack.A(stack),compute=TRUE),verbose=FALSE)

# process data for plotting
result<-bind_cols(pred_coords,output$summary.fitted.values[index_pred,])
result$latitude = result$x2
result$longitude = result$x1
pred_df<-result %>%  dplyr::select(latitude,longitude,mean)
pred_sf<-sf::st_as_sf(pred_df,coords=c("longitude","latitude"),crs=4326)

# plot result
g<-ggplot()+
  geom_polygon(data=legionmap,fill="white", color="black", aes(x=long, y=lat, group=group))+
  geom_sf(data=pred_sf,aes(color=mean),size=0.3,alpha=0.01)+
  
  ylim(c(35.5,NA))+
  #xlim(c(139.7,NA))+
  annotate("text",x=139.735,y=35.72,label="x",color='black')+
  scale_color_gradient2(limits=c(-0.0015,0.002))+
  scale_fill_continuous(type = "gradient")+
  theme_bw()
plot(g)

