rm(list = ls())

# install.packages(c('ggforce', 'cowplot', 'sphereplot', 'rworldmap', 'ggplot2', 'sf', 'ncdf4', 'dplyr', 'raster', 'gridExtra', 'genlasso', 'mgcv'))

library(ggforce)
library(cowplot)
library(sphereplot)
library(rworldmap)
library(ggplot2)
library(sf)
library(ncdf4)
library(ggplot2)
library(dplyr)
library(raster)
library(gridExtra)
library(stats)
library(genlasso)
library(mgcv)

source("admm_source.R")
source("splineBox1.6_JJH_240719.R")
source("iBspline.2.R")

#==============================================================================
# Data & Fitting (Integrated)
#==============================================================================
test1 = read.csv("Revision(25.05)_data_ldy/drifter.csv")
t1 = subset(test1,test1$WMO == 5102764)
t2 = subset(test1,test1$WMO == 5102765)
t3 = subset(test1,test1$WMO == 5102766)

n = nrow(t2)
order = 3
dimension = 40
t = seq(0, 1, length = n)

t1 = cbind(t1,t[1:nrow(t1)])
t2 = cbind(t2,t[1:nrow(t2)])
t3 = cbind(t3,t[1:nrow(t3)])
t1 = t1[, -which(names(t1) == "WMO")]
t2 = t2[, -which(names(t2) == "WMO")]
t3 = t3[, -which(names(t3) == "WMO")]
names(t2) = names(t1)
names(t3) = names(t1)


t123 = rbind(t1,t2,t3)
colnames(t123)[3] = 't'
t123_sorted = t123[order(t123$t), ]

y = as.matrix(sapply(t123_sorted[, c("longitude", "latitude")], as.numeric))

n = nrow(t123_sorted)
t = t123_sorted$t

#==============================================================================

knots = knots_quantile(t, dimension, order)
B = bsplines(t, knots, order)
D = bspline_jump(knots, order)
fit = bspline.curve.admm_lambdas(y, D, B,
                                 lambdas = NULL,
                                 lam_max = 100,
                                 lam_min = 1e-10,
                                 n_lambda = 200,
                                 max_iter = 1000,
                                 epsilon = 1e-8,
                                 eta_c = 1)

worldMap = getMap()
worldMap_sf = st_as_sf(worldMap)
mar = 0.2


y_df = as.data.frame(y)
y_df_sf = st_as_sf(y_df, coords = c("longitude", "latitude"), crs = 4326)

t_new = seq(0, 1, length = 5000)
B_new = bsplines(t_new, knots, order)
# best aic
fit_new = B_new %*% fit[[which.min(fit$aic)]]$xi
fit_new_df = as.data.frame(fit_new)
fit_new_df_sf = st_as_sf(fit_new_df, coords = c("longitude", "latitude"), crs = 4326)
# best bic
fit_new2 = B_new %*% fit[[which.min(fit$bic)]]$xi
fit_new_df2 = as.data.frame(fit_new2)
fit_new_df_sf2 = st_as_sf(fit_new_df2, coords = c("longitude", "latitude"), crs = 4326)

#==============================================================================
# [Plot1(1,1)] World map
#==============================================================================
nino4_center = c(-155, 4)  # Niño 4 region center (longitude, latitude)
nino4_radius = 20         # Radius (unit is latitude/longitude)

worldmap = ggplot() +
  geom_sf(data = worldMap_sf, color = "grey", fill = "antiquewhite") +
  geom_sf(data = y_df_sf, size = 0.1) +
  geom_circle(
    aes(x0 = nino4_center[1], y0 = nino4_center[2], r = nino4_radius),
    inherit.aes = FALSE,  # Operates independently of existing data
    color = "red",        # Circle border color
    linetype = "dashed",  # Circle style
    # size = 1.2            # Circle thickness
    linewidth = 1.2
  ) +
  coord_sf(
    xlim = c(-180, -40),  # Longitude range
    ylim = c(-70, 70)     # Latitude range
  ) +
  theme(panel.background = element_rect(fill = "lightblue", color = NA)) # Set background color
worldmap

#==============================================================================
# [Plot3(1,1)] best aic
#==============================================================================
zoommap1 = ggplot() +
  geom_sf(data = worldMap_sf, aes(fill = "World Map"), color = "grey", fill = "antiquewhite") +
  geom_sf(data = y_df_sf, aes(color = "Data Points"), size = 1) +
  geom_sf(data = fit_new_df_sf, aes(color = "Fitted Line"), size = 1) +
  coord_sf(xlim = c(min(y_df$longitude) - mar,
                    max(y_df$longitude) + mar),
           ylim = c(min(y_df$latitude) - mar,
                    max(y_df$latitude) + mar)) +
  theme(panel.background = element_rect(fill = "lightblue", color = NA),
        legend.position = c(0.2, 0.2),
        legend.background = element_rect(fill = "white", color = "white"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18)) +
  scale_fill_manual(values = c("World Map" = "antiquewhite")) +
  scale_color_manual(values = c("Data Points" = "slategrey", "Fitted Line" = "red")) +
  labs(fill = "Map Fill",color = 'Color') # Modify legend title
zoommap1

#==============================================================================
# [Plot3(1,2)] best bic
#==============================================================================
zoommap2 = ggplot() +
  geom_sf(data = worldMap_sf, aes(fill = "World Map"), color = "grey", fill = "antiquewhite") +
  geom_sf(data = y_df_sf, aes(color = "Data Points"), size = 1) +
  geom_sf(data = fit_new_df_sf2, aes(color = "Fitted Line"), size = 1) +
  coord_sf(xlim = c(min(y_df$longitude) - mar,
                    max(y_df$longitude) + mar),
           ylim = c(min(y_df$latitude) - mar,
                    max(y_df$latitude) + mar)) +
  theme(panel.background = element_rect(fill = "lightblue", color = NA),
        legend.position = c(0.2, 0.2), # Move legend to the bottom right inside the plot
        legend.background = element_rect(fill = "white", color = "white"), # Set legend background color and border
        legend.title = element_text(size = 20), # Legend title size
        legend.text = element_text(size = 18)) + # Legend text size
  scale_fill_manual(values = c("World Map" = "antiquewhite")) +
  scale_color_manual(values = c("Data Points" = "slategrey", "Fitted Line" = "red")) +
  labs(fill = "Legend", color = "Color") # Set legend title
zoommap2

#==============================================================================
# Data & Fitting (Individual)
#==============================================================================
b1 = subset(test1,test1$WMO == 5102764)
b2 = subset(test1,test1$WMO == 5102765)
b3 = subset(test1,test1$WMO == 5102766)


n1 = nrow(b1)
n2 = nrow(b2)
n3 = nrow(b3)

time1 = seq(0, 1, length = n1)
time2 = seq(0, 1, length = n2)
time3 = seq(0, 1, length = n3)

b1 = cbind(b1,time1)
b2 = cbind(b2,time2)
b3 = cbind(b3,time3)
b1 = b1[, -which(names(b1) == "WMO")]
b2 = b2[, -which(names(b2) == "WMO")]
b3 = b3[, -which(names(b3) == "WMO")]


y1 = as.matrix(sapply(b1[,c("longitude", "latitude")], as.numeric))
y2 = as.matrix(sapply(b2[,c("longitude", "latitude")], as.numeric))
y3 = as.matrix(sapply(b3[,c("longitude", "latitude")], as.numeric))

order = 3
dimension = 10

knots1 = knots_quantile(time1, dimension, order)
B1 = bsplines(time1, knots1, order)
D1 = bspline_jump(knots1, order)
fit1 = bspline.curve.admm_lambdas(y1, D1, B1,
                                  lambdas = NULL,
                                  lam_max = 100,
                                  lam_min = 1e-10,
                                  n_lambda = 200,
                                  max_iter = 1000,
                                  epsilon = 1e-8,
                                  eta_c = 1)
knots2 = knots_quantile(time2, dimension, order)
B2 = bsplines(time2, knots2, order)
D2 = bspline_jump(knots2, order)
fit2 = bspline.curve.admm_lambdas(y2, D2, B2,
                                  lambdas = NULL,
                                  lam_max = 100,
                                  lam_min = 1e-10,
                                  n_lambda = 200,
                                  max_iter = 1000,
                                  epsilon = 1e-8,
                                  eta_c = 1)
knots3 = knots_quantile(time3, dimension, order)
B3 = bsplines(time3, knots3, order)
D3 = bspline_jump(knots3, order)
fit3 = bspline.curve.admm_lambdas(y3, D3, B3,
                                  lambdas = NULL,
                                  lam_max = 100,
                                  lam_min = 1e-10,
                                  n_lambda = 200,
                                  max_iter = 1000,
                                  epsilon = 1e-8,
                                  eta_c = 1)

worldMap = getMap()
worldMap_sf = st_as_sf(worldMap)

y1_df = as.data.frame(y1)
y2_df = as.data.frame(y2)
y3_df = as.data.frame(y3)
y1_df_sf = st_as_sf(y1_df, coords = c("longitude", "latitude"), crs = 4326)
y2_df_sf = st_as_sf(y2_df, coords = c("longitude", "latitude"), crs = 4326)
y3_df_sf = st_as_sf(y3_df, coords = c("longitude", "latitude"), crs = 4326)

t_new = seq(0, 1, length = 5000)
B1_new = bsplines(t_new, knots1, order)
B2_new = bsplines(t_new, knots2, order)
B3_new = bsplines(t_new, knots3, order)

# best aic
fit1_new = B1_new %*% fit1[[which.min(fit1$aic)]]$xi
fit2_new = B2_new %*% fit2[[which.min(fit2$aic)]]$xi
fit3_new = B3_new %*% fit3[[which.min(fit3$aic)]]$xi
fit1_new_df = as.data.frame(fit1_new)
fit2_new_df = as.data.frame(fit2_new)
fit3_new_df = as.data.frame(fit3_new)
fit1_new_df_sf = st_as_sf(fit1_new_df, coords = c("longitude", "latitude"), crs = 4326)
fit2_new_df_sf = st_as_sf(fit2_new_df, coords = c("longitude", "latitude"), crs = 4326)
fit3_new_df_sf = st_as_sf(fit3_new_df, coords = c("longitude", "latitude"), crs = 4326)

# best bic (for reference)
fit1_new2 = B1_new %*% fit1[[which.min(fit1$bic)]]$xi
fit2_new2 = B2_new %*% fit2[[which.min(fit2$bic)]]$xi
fit3_new2 = B3_new %*% fit3[[which.min(fit3$bic)]]$xi
fit1_new_df2 = as.data.frame(fit1_new2)
fit2_new_df2 = as.data.frame(fit2_new2)
fit3_new_df2 = as.data.frame(fit3_new2)
fit1_new_df_sf2 = st_as_sf(fit1_new_df2, coords = c("longitude", "latitude"), crs = 4326)
fit2_new_df_sf2 = st_as_sf(fit2_new_df2, coords = c("longitude", "latitude"), crs = 4326)
fit3_new_df_sf2 = st_as_sf(fit3_new_df2, coords = c("longitude", "latitude"), crs = 4326)

#==============================================================================
# [Figure1(1,2)] basic point default
#==============================================================================
mar = 0.2

b123basic = ggplot() +
  geom_sf(data = worldMap_sf, aes(fill = "World Map"), color = "grey", fill = "antiquewhite") +
  geom_sf(data = y1_df_sf, aes(color = "5102764"), size = 2) +
  geom_sf(data = y2_df_sf, aes(color = "5102765"), size = 2) +
  geom_sf(data = y3_df_sf, aes(color = "5102766"), size = 2) +
  coord_sf(xlim = c(min(y_df$longitude) - mar,
                    max(y_df$longitude) + mar),
           ylim = c(min(y_df$latitude) - mar,
                    max(y_df$latitude) + mar)) +
  theme(panel.background = element_rect(fill = "lightblue", color = NA),
        legend.position = c(0.2, 0.2), # Move legend to the bottom left inside the plot
        legend.background = element_rect(fill = "white", color = "white"), # Set legend background color and border
        legend.title = element_text(size = 20), # Legend title size
        legend.text = element_text(size = 18)) + # Legend text size
  scale_fill_manual(values = c("World Map" = "antiquewhite")) +
  scale_color_manual(values = c("5102764" = "red",
                                "5102765" = "blue",
                                "5102766" = "black")) +
  labs(fill = "Map Fill", color = "Drifters") # Set legend title
b123basic

#==============================================================================
# [Figure 2(1,1) ~ 2(2,2)] Total 4 plots
# Plot b123 fitted lines together (AIC)
#==============================================================================
b123fitted = ggplot() +
  geom_sf(data = worldMap_sf, aes(fill = "World Map"), color = "grey", fill = "antiquewhite") +
  geom_sf(data = y1_df_sf, aes(color = "5102764"), size = 1) +
  geom_sf(data = y2_df_sf, aes(color = "5102765"), size = 1) +
  geom_sf(data = y3_df_sf, aes(color = "5102766"), size = 1) +
  geom_sf(data = fit1_new_df_sf, aes(color = "5102764"), size = 1) +
  geom_sf(data = fit2_new_df_sf, aes(color = "5102765"), size = 1) +
  geom_sf(data = fit3_new_df_sf, aes(color = "5102766"), size = 1) +
  coord_sf(xlim = c(min(y_df$longitude) - mar,
                    max(y_df$longitude) + mar),
           ylim = c(min(y_df$latitude) - mar,
                    max(y_df$latitude) + mar)) +
  theme(panel.background = element_rect(fill = "lightblue", color = NA),
        legend.position = c(0.2, 0.2), # Move legend to the bottom right inside the plot
        legend.background = element_rect(fill = "white", color = "white"), # Set legend background color and border
        legend.title = element_text(size = 20), # Legend title size
        legend.text = element_text(size = 18)) + # Legend text size
  scale_fill_manual(values = c("World Map" = "antiquewhite")) +
  scale_color_manual(values = c("5102764" = "red",
                                "5102765" = "blue",
                                "5102766" = "black")) +
  labs(color = "Drifters") # Set single legend title
b123fitted


# b1 fitted line (AIC)
b1fitted = ggplot() +
  geom_sf(data = worldMap_sf, aes(fill = "World Map"), color = "grey", fill = "antiquewhite") +
  geom_sf(data = y1_df_sf, aes(color = "Data Points"), size = 1) +
  geom_sf(data = fit1_new_df_sf, aes(color = "Fitted Line"), size = 1) +
  coord_sf(xlim = c(min(y_df$longitude) - mar,
                    max(y_df$longitude) + mar),
           ylim = c(min(y_df$latitude) - mar,
                    max(y_df$latitude) + mar)) +
  theme(panel.background = element_rect(fill = "lightblue", color = NA),
        legend.position = c(0.2, 0.2),
        legend.background = element_rect(fill = "white", color = "white"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18)) +
  scale_fill_manual(values = c("World Map" = "antiquewhite")) +
  scale_color_manual(values = c("Data Points" = "slategrey", "Fitted Line" = "red")) +
  labs(fill = "Map Fill", color = 'Color') # Set legend title
b1fitted


# b2 fitted line (AIC)
b2fitted = ggplot() +
  geom_sf(data = worldMap_sf, aes(fill = "World Map"), color = "grey", fill = "antiquewhite") +
  geom_sf(data = y2_df_sf, aes(color = "Data Points"), size = 1) +
  geom_sf(data = fit2_new_df_sf, aes(color = "Fitted Line"), size = 1) +
  coord_sf(xlim = c(min(y_df$longitude) - mar,
                    max(y_df$longitude) + mar),
           ylim = c(min(y_df$latitude) - mar,
                    max(y_df$latitude) + mar)) +
  theme(panel.background = element_rect(fill = "lightblue", color = NA),
        legend.position = c(0.2, 0.2),
        legend.background = element_rect(fill = "white", color = "white"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18)) +
  scale_fill_manual(values = c("World Map" = "antiquewhite")) +
  scale_color_manual(values = c("Data Points" = "slategrey", "Fitted Line" = "blue")) +
  labs(fill = "Map Fill", color = 'Color') # Set legend title
b2fitted


# b3 fitted line (AIC)
b3fitted = ggplot() +
  geom_sf(data = worldMap_sf, aes(fill = "World Map"), color = "grey", fill = "antiquewhite") +
  geom_sf(data = y3_df_sf, aes(color = "Data Points"), size = 1) +
  geom_sf(data = fit3_new_df_sf, aes(color = "Fitted Line"), size = 1) +
  coord_sf(xlim = c(min(y_df$longitude) - mar,
                    max(y_df$longitude) + mar),
           ylim = c(min(y_df$latitude) - mar,
                    max(y_df$latitude) + mar)) +
  theme(panel.background = element_rect(fill = "lightblue", color = NA),
        legend.position = c(0.2, 0.2),
        legend.background = element_rect(fill = "white", color = "white"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18)) +
  scale_fill_manual(values = c("World Map" = "antiquewhite")) +
  scale_color_manual(values = c("Data Points" = "slategrey", "Fitted Line" = "black")) +
  labs(fill = "Map Fill", color = 'Color') # Set legend title
b3fitted

#==============================================================================
# Figure 2
#==============================================================================
# 200405 Wind & Current
#==============================================================================
nc_file = "Revision(25.05)_data_ldy/cmems_mod_glo_wav_my_0.2deg_PT3H-i_1739711784083.nc"
nc = nc_open(nc_file)

# Get variables
lon = ncvar_get(nc, "longitude")
lat = ncvar_get(nc, "latitude")
wind_dir = ncvar_get(nc, "VMDR_WW")  # Wind Wave Direction
wave_dir = ncvar_get(nc, "VMDR")    # Mean Wave Direction
stokes_x = ncvar_get(nc, "VSDX")    # Stokes Drift X (current)
stokes_y = ncvar_get(nc, "VSDY")    # Stokes Drift Y (current)
wave_height = ncvar_get(nc, "VHM0") # Significant Wave Height

# Close file
nc_close(nc)

# Convert to dataframe
df = expand.grid(lon = lon, lat = lat)
df$wind_dir = as.vector(wind_dir[,,1])   # Wind direction (azimuth)
df$wave_dir = as.vector(wave_dir[,,1])   # Wave direction
df$stokes_x = as.vector(stokes_x[,,1])   # Stokes Drift X component
df$stokes_y = as.vector(stokes_y[,,1])   # Stokes Drift Y component
df$wave_height = as.vector(wave_height[,,1]) # Wave height

# Convert wind direction to vector (u, v)
df$wind_u = cos(df$wind_dir * pi / 180)  # x component
df$wind_v = sin(df$wind_dir * pi / 180)  # y component

# Map visualization
p1 = ggplot(df, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = wave_height), alpha = 0.5) +  # Background color: wave height
  scale_fill_viridis_c(name = "wave_height(m)") +
  geom_segment(aes(xend = lon + 0.2 * wind_u, yend = lat + 0.2 * wind_v),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", alpha = 0.8, size = 0.5) +  # Wind direction (black arrow)
  geom_segment(aes(xend = lon + 0.2 * stokes_x, yend = lat + 0.2 * stokes_y),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "red", alpha = 0.8, size = 0.5) +  # Stokes Drift (current direction, red arrow)
  coord_quickmap() +
  labs(title = "200405 Wind & Current", x = "longitude", y = "latitude") +
  theme_minimal()
p1

#==============================================================================
# 200405 Primary Swell & Current
#==============================================================================
nc_file =  "Revision(25.05)_data_ldy/cmems_mod_glo_wav_my_0.2deg_PT3H-i_1739713922090.nc"
nc = nc_open(nc_file)

# Get variables
lon = ncvar_get(nc, "longitude")
lat = ncvar_get(nc, "latitude")
swell_dir = ncvar_get(nc, "VMDR_SW1")    # Primary swell direction
swell_height = ncvar_get(nc, "VHM0_SW1") # Primary swell height
swell_period = ncvar_get(nc, "VTM01_SW1")# Primary swell period
stokes_x = ncvar_get(nc, "VSDX")         # Stokes Drift X (current)
stokes_y = ncvar_get(nc, "VSDY")         # Stokes Drift Y (current)
wave_height = ncvar_get(nc, "VHM0")      # Total wave height (for reference)

# Close file
nc_close(nc)

# Convert to dataframe
df = expand.grid(lon = lon, lat = lat)
df$swell_dir = as.vector(swell_dir[,,1])    # Primary swell direction (azimuth)
df$swell_height = as.vector(swell_height[,,1]) # Primary swell height
df$swell_period = as.vector(swell_period[,,1]) # Primary swell period
df$stokes_x = as.vector(stokes_x[,,1])      # Stokes Drift X component
df$stokes_y = as.vector(stokes_y[,,1])      # Stokes Drift Y component
df$wave_height = as.vector(wave_height[,,1])  # Total wave height

# Convert swell direction to vector (u, v)
df$swell_u = cos(df$swell_dir * pi / 180)  # x component
df$swell_v = sin(df$swell_dir * pi / 180)  # y component

# Map visualization
p2 = ggplot(df, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = swell_height), alpha = 0.5) +  # Background color: primary swell height
  scale_fill_viridis_c(name = "swell_height(m)") +
  geom_segment(aes(xend = lon + 0.2 * swell_u, yend = lat + 0.2 * swell_v),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "blue", alpha = 0.8, size = 0.5) +  # Primary swell direction (blue arrow)
  geom_segment(aes(xend = lon + 0.2 * stokes_x, yend = lat + 0.2 * stokes_y),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "red", alpha = 0.8, size = 0.5) +  # Stokes Drift (current direction, red arrow)
  coord_quickmap() +
  labs(title = "200405 Primary Swell & Current", x = "longitude", y = "latitude") +
  theme_minimal()
p2

#==============================================================================
# 200405 Secondary Swell & Current
#==============================================================================
nc_file = "Revision(25.05)_data_ldy/cmems_mod_glo_wav_my_0.2deg_PT3H-i_1739714468705.nc"
nc = nc_open(nc_file)

# Get variables
lon = ncvar_get(nc, "longitude")
lat = ncvar_get(nc, "latitude")
swell_dir = ncvar_get(nc, "VMDR_SW2")    # Secondary swell direction
swell_height = ncvar_get(nc, "VHM0_SW2") # Secondary swell height
swell_period = ncvar_get(nc, "VTM01_SW2")# Secondary swell period
stokes_x = ncvar_get(nc, "VSDX")         # Stokes Drift X (current)
stokes_y = ncvar_get(nc, "VSDY")         # Stokes Drift Y (current)
wave_height = ncvar_get(nc, "VHM0")      # Total wave height (for reference)

# Close file
nc_close(nc)

# Convert to dataframe
df = expand.grid(lon = lon, lat = lat)
df$swell_dir = as.vector(swell_dir[,,1])    # Secondary swell direction (azimuth)
df$swell_height = as.vector(swell_height[,,1]) # Secondary swell height
df$swell_period = as.vector(swell_period[,,1]) # Secondary swell period
df$stokes_x = as.vector(stokes_x[,,1])      # Stokes Drift X component
df$stokes_y = as.vector(stokes_y[,,1])      # Stokes Drift Y component
df$wave_height = as.vector(wave_height[,,1])  # Total wave height

# Convert swell direction to vector (u, v)
df$swell_u = cos(df$swell_dir * pi / 180)  # x component
df$swell_v = sin(df$swell_dir * pi / 180)  # y component

# Map visualization
p3 = ggplot(df, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = swell_height), alpha = 0.5) +  # Background color: secondary swell height
  scale_fill_viridis_c(name = "swell_height(m)") +
  geom_segment(aes(xend = lon + 0.2 * swell_u, yend = lat + 0.2 * swell_v),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "purple", alpha = 0.8, size = 0.5) +  # Secondary swell direction (purple arrow)
  geom_segment(aes(xend = lon + 0.2 * stokes_x, yend = lat + 0.2 * stokes_y),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "red", alpha = 0.8, size = 0.5) +  # Stokes Drift (current direction, red arrow)
  coord_quickmap() +
  labs(title = "200405 Secondary Swell & Current", x = "longitude", y = "latitude") +
  theme_minimal()
p3

#==============================================================================
# 200421 Wind & Current
#==============================================================================
nc_file = "Revision(25.05)_data_ldy/cmems_mod_glo_wav_my_0.2deg_PT3H-i_1739805996585.nc"
nc = nc_open(nc_file)

# Get variables
lon = ncvar_get(nc, "longitude")
lat = ncvar_get(nc, "latitude")
wind_dir = ncvar_get(nc, "VMDR_WW")  # Wind Wave Direction
wave_dir = ncvar_get(nc, "VMDR")    # Mean Wave Direction
stokes_x = ncvar_get(nc, "VSDX")    # Stokes Drift X (current)
stokes_y = ncvar_get(nc, "VSDY")    # Stokes Drift Y (current)
wave_height = ncvar_get(nc, "VHM0") # Significant Wave Height

# Close file
nc_close(nc)

# Convert to dataframe
df = expand.grid(lon = lon, lat = lat)
df$wind_dir = as.vector(wind_dir[,,1])   # Wind direction (azimuth)
df$wave_dir = as.vector(wave_dir[,,1])   # Wave direction
df$stokes_x = as.vector(stokes_x[,,1])   # Stokes Drift X component
df$stokes_y = as.vector(stokes_y[,,1])   # Stokes Drift Y component
df$wave_height = as.vector(wave_height[,,1]) # Wave height

# Convert wind direction to vector (u, v)
df$wind_u = cos(df$wind_dir * pi / 180)  # x component
df$wind_v = sin(df$wind_dir * pi / 180)  # y component

# Map visualization
p4 = ggplot(df, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = wave_height), alpha = 0.5) +  # Background color: wave height
  scale_fill_viridis_c(name = "wave_height(m)") +
  geom_segment(aes(xend = lon + 0.2 * wind_u, yend = lat + 0.2 * wind_v),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", alpha = 0.8, size = 0.5) +  # Wind direction (black arrow)
  geom_segment(aes(xend = lon + 0.2 * stokes_x, yend = lat + 0.2 * stokes_y),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "red", alpha = 0.8, size = 0.5) +  # Stokes Drift (current direction, red arrow)
  coord_quickmap() +
  labs(title = "200421 Wind & Current", x = "longitude", y = "latitude") +
  theme_minimal()
p4

#==============================================================================
# 200421 Primary Swell & Current
#==============================================================================
nc_file = "Revision(25.05)_data_ldy/cmems_mod_glo_wav_my_0.2deg_PT3H-i_1739806290060.nc"
nc = nc_open(nc_file)

# Get variables
lon = ncvar_get(nc, "longitude")
lat = ncvar_get(nc, "latitude")
swell_dir = ncvar_get(nc, "VMDR_SW1")    # Primary swell direction
swell_height = ncvar_get(nc, "VHM0_SW1") # Primary swell height
swell_period = ncvar_get(nc, "VTM01_SW1")# Primary swell period
stokes_x = ncvar_get(nc, "VSDX")         # Stokes Drift X (current)
stokes_y = ncvar_get(nc, "VSDY")         # Stokes Drift Y (current)
wave_height = ncvar_get(nc, "VHM0")      # Total wave height (for reference)

# Close file
nc_close(nc)

# Convert to dataframe
df = expand.grid(lon = lon, lat = lat)
df$swell_dir = as.vector(swell_dir[,,1])    # Primary swell direction (azimuth)
df$swell_height = as.vector(swell_height[,,1]) # Primary swell height
df$swell_period = as.vector(swell_period[,,1]) # Primary swell period
df$stokes_x = as.vector(stokes_x[,,1])      # Stokes Drift X component
df$stokes_y = as.vector(stokes_y[,,1])      # Stokes Drift Y component
df$wave_height = as.vector(wave_height[,,1])  # Total wave height

# Convert swell direction to vector (u, v)
df$swell_u = cos(df$swell_dir * pi / 180)  # x component
df$swell_v = sin(df$swell_dir * pi / 180)  # y component

# Map visualization
p5 = ggplot(df, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = swell_height), alpha = 0.5) +  # Background color: primary swell height
  scale_fill_viridis_c(name = "swell_height(m)") +
  geom_segment(aes(xend = lon + 0.2 * swell_u, yend = lat + 0.2 * swell_v),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "blue", alpha = 0.8, size = 0.5) +  # Primary swell direction (blue arrow)
  geom_segment(aes(xend = lon + 0.2 * stokes_x, yend = lat + 0.2 * stokes_y),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "red", alpha = 0.8, size = 0.5) +  # Stokes Drift (current direction, red arrow)
  coord_quickmap() +
  labs(title = "200421 Primary Swell & Current", x = "longitude", y = "latitude") +
  theme_minimal()
p5

#==============================================================================
# 200421 Secondary Swell & Current
#==============================================================================
nc_file = "Revision(25.05)_data_ldy/cmems_mod_glo_wav_my_0.2deg_PT3H-i_1739806388750.nc"
nc = nc_open(nc_file)

# Get variables
lon = ncvar_get(nc, "longitude")
lat = ncvar_get(nc, "latitude")
swell_dir = ncvar_get(nc, "VMDR_SW2")    # Secondary swell direction
swell_height = ncvar_get(nc, "VHM0_SW2") # Secondary swell height
swell_period = ncvar_get(nc, "VTM01_SW2")# Secondary swell period
stokes_x = ncvar_get(nc, "VSDX")         # Stokes Drift X (current)
stokes_y = ncvar_get(nc, "VSDY")         # Stokes Drift Y (current)
wave_height = ncvar_get(nc, "VHM0")      # Total wave height (for reference)

# Close file
nc_close(nc)

# Convert to dataframe
df = expand.grid(lon = lon, lat = lat)
df$swell_dir = as.vector(swell_dir[,,1])    # Secondary swell direction (azimuth)
df$swell_height = as.vector(swell_height[,,1]) # Secondary swell height
df$swell_period = as.vector(swell_period[,,1]) # Secondary swell period
df$stokes_x = as.vector(stokes_x[,,1])      # Stokes Drift X component
df$stokes_y = as.vector(stokes_y[,,1])      # Stokes Drift Y component
df$wave_height = as.vector(wave_height[,,1])  # Total wave height

# Convert swell direction to vector (u, v)
df$swell_u = cos(df$swell_dir * pi / 180)  # x component
df$swell_v = sin(df$swell_dir * pi / 180)  # y component

# Map visualization
p6 = ggplot(df, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = swell_height), alpha = 0.5) +  # Background color: secondary swell height
  scale_fill_viridis_c(name = "swell_height(m)") +
  geom_segment(aes(xend = lon + 0.2 * swell_u, yend = lat + 0.2 * swell_v),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "purple", alpha = 0.8, size = 0.5) +  # Secondary swell direction (purple arrow)
  geom_segment(aes(xend = lon + 0.2 * stokes_x, yend = lat + 0.2 * stokes_y),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "red", alpha = 0.8, size = 0.5) +  # Stokes Drift (current direction, red arrow)
  coord_quickmap() +
  labs(title = "200421 Secondary Swell & Current", x = "longitude", y = "latitude") +
  theme_minimal()
p6

#==============================================================================
# 200509 Wind & Current
#==============================================================================
nc_file = "Revision(25.05)_data_ldy/cmems_mod_glo_wav_my_0.2deg_PT3H-i_1739806673914.nc"
nc = nc_open(nc_file)

# Get variables
lon = ncvar_get(nc, "longitude")
lat = ncvar_get(nc, "latitude")
wind_dir = ncvar_get(nc, "VMDR_WW")  # Wind Wave Direction
wave_dir = ncvar_get(nc, "VMDR")    # Mean Wave Direction
stokes_x = ncvar_get(nc, "VSDX")    # Stokes Drift X (current)
stokes_y = ncvar_get(nc, "VSDY")    # Stokes Drift Y (current)
wave_height = ncvar_get(nc, "VHM0") # Significant Wave Height

# Close file
nc_close(nc)

# Convert to dataframe
df = expand.grid(lon = lon, lat = lat)
df$wind_dir = as.vector(wind_dir[,,1])   # Wind direction (azimuth)
df$wave_dir = as.vector(wave_dir[,,1])   # Wave direction
df$stokes_x = as.vector(stokes_x[,,1])   # Stokes Drift X component
df$stokes_y = as.vector(stokes_y[,,1])   # Stokes Drift Y component
df$wave_height = as.vector(wave_height[,,1]) # Wave height

# Convert wind direction to vector (u, v)
df$wind_u = cos(df$wind_dir * pi / 180)  # x component
df$wind_v = sin(df$wind_dir * pi / 180)  # y component


# Map visualization
p7 = ggplot(df, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = wave_height), alpha = 0.5) +  # Background color: wave height
  scale_fill_viridis_c(name = "wave_height(m)") +
  geom_segment(aes(xend = lon + 0.2 * wind_u, yend = lat + 0.2 * wind_v),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", alpha = 0.8, size = 0.5) +  # Wind direction (black arrow)
  geom_segment(aes(xend = lon + 0.2 * stokes_x, yend = lat + 0.2 * stokes_y),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "red", alpha = 0.8, size = 0.5) +  # Stokes Drift (current direction, red arrow)
  coord_quickmap() +
  labs(title = "200509 Wind & Current", x = "longitude", y = "latitude") +
  theme_minimal()
p7

#==============================================================================
# 200509 Primary Swell & Current
#==============================================================================
nc_file = "Revision(25.05)_data_ldy/cmems_mod_glo_wav_my_0.2deg_PT3H-i_1739807093609.nc"
nc = nc_open(nc_file)

# Get variables
lon = ncvar_get(nc, "longitude")
lat = ncvar_get(nc, "latitude")
swell_dir = ncvar_get(nc, "VMDR_SW1")    # Primary swell direction
swell_height = ncvar_get(nc, "VHM0_SW1") # Primary swell height
swell_period = ncvar_get(nc, "VTM01_SW1")# Primary swell period
stokes_x = ncvar_get(nc, "VSDX")         # Stokes Drift X (current)
stokes_y = ncvar_get(nc, "VSDY")         # Stokes Drift Y (current)
wave_height = ncvar_get(nc, "VHM0")      # Total wave height (for reference)

# Close file
nc_close(nc)

# Convert to dataframe
df = expand.grid(lon = lon, lat = lat)
df$swell_dir = as.vector(swell_dir[,,1])    # Primary swell direction (azimuth)
df$swell_height = as.vector(swell_height[,,1]) # Primary swell height
df$swell_period = as.vector(swell_period[,,1]) # Primary swell period
df$stokes_x = as.vector(stokes_x[,,1])      # Stokes Drift X component
df$stokes_y = as.vector(stokes_y[,,1])      # Stokes Drift Y component
df$wave_height = as.vector(wave_height[,,1])  # Total wave height

# Convert swell direction to vector (u, v)
df$swell_u = cos(df$swell_dir * pi / 180)  # x component
df$swell_v = sin(df$swell_dir * pi / 180)  # y component

# Map visualization
p8 = ggplot(df, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = swell_height), alpha = 0.5) +  # Background color: primary swell height
  scale_fill_viridis_c(name = "swell_height(m)") +
  geom_segment(aes(xend = lon + 0.2 * swell_u, yend = lat + 0.2 * swell_v),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "blue", alpha = 0.8, size = 0.5) +  # Primary swell direction (blue arrow)
  geom_segment(aes(xend = lon + 0.2 * stokes_x, yend = lat + 0.2 * stokes_y),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "red", alpha = 0.8, size = 0.5) +  # Stokes Drift (current direction, red arrow)
  coord_quickmap() +
  labs(title = "200509 Primary Swell & Current", x = "longitude", y = "latitude") +
  theme_minimal()
p8

#==============================================================================
# 200509 Secondary Swell & Current
#==============================================================================
nc_file = "Revision(25.05)_data_ldy/cmems_mod_glo_wav_my_0.2deg_PT3H-i_1739807415115.nc"
nc = nc_open(nc_file)

# Get variables
lon = ncvar_get(nc, "longitude")
lat = ncvar_get(nc, "latitude")
swell_dir = ncvar_get(nc, "VMDR_SW2")    # Secondary swell direction
swell_height = ncvar_get(nc, "VHM0_SW2") # Secondary swell height
swell_period = ncvar_get(nc, "VTM01_SW2")# Secondary swell period
stokes_x = ncvar_get(nc, "VSDX")         # Stokes Drift X (current)
stokes_y = ncvar_get(nc, "VSDY")         # Stokes Drift Y (current)
wave_height = ncvar_get(nc, "VHM0")      # Total wave height (for reference)

# Close file
nc_close(nc)

# Convert to dataframe
df = expand.grid(lon = lon, lat = lat)
df$swell_dir = as.vector(swell_dir[,,1])    # Secondary swell direction (azimuth)
df$swell_height = as.vector(swell_height[,,1]) # Secondary swell height
df$swell_period = as.vector(swell_period[,,1]) # Secondary swell period
df$stokes_x = as.vector(stokes_x[,,1])      # Stokes Drift X component
df$stokes_y = as.vector(stokes_y[,,1])      # Stokes Drift Y component
df$wave_height = as.vector(wave_height[,,1])  # Total wave height

# Convert swell direction to vector (u, v)
df$swell_u = cos(df$swell_dir * pi / 180)  # x component
df$swell_v = sin(df$swell_dir * pi / 180)  # y component

# Map visualization
p9 = ggplot(df, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = swell_height), alpha = 0.5) +  # Background color: secondary swell height
  scale_fill_viridis_c(name = "swell_height(m)") +
  geom_segment(aes(xend = lon + 0.2 * swell_u, yend = lat + 0.2 * swell_v),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "purple", alpha = 0.8, size = 0.5) +  # Secondary swell direction (purple arrow)
  geom_segment(aes(xend = lon + 0.2 * stokes_x, yend = lat + 0.2 * stokes_y),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "red", alpha = 0.8, size = 0.5) +  # Stokes Drift (current direction, red arrow)
  coord_quickmap() +
  labs(title = "200509 Secondary Swell & Current", x = "longitude", y = "latitude") +
  theme_minimal()
p9

#==============================================================================
# Figure 2
#==============================================================================
p0 = grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9, nrow = 3, ncol = 3)

# pdf("output.pdf", width=12, height=8)
# par(mfrow=c(3,3), mar=c(4,4,2,2), oma=c(1,1,1,1))
# plot(p0)

# pdf("output.pdf", width = 12, height = 8, bg = "white")
# par(mfrow = c(3, 3), mar = c(4, 4, 2, 2), oma = c(1, 1, 1, 1))
# grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3)

#==============================================================================
# [Table1] Bsplinecurve_MSE,MAE,MXDV
#==============================================================================
n_sim <- 50
results <- data.frame(seed = integer(n_sim), mse = numeric(n_sim), mae = numeric(n_sim), mxdv = numeric(n_sim),
                      mse_ind = numeric(n_sim), mae_ind = numeric(n_sim), mxdv_ind = numeric(n_sim)
)

for(s in 1:n_sim) {
  set.seed(s)
  
  fixed_indices <- c(1, n)
  remaining_indices <- 2:(n-1)
  
  n_train_target <- floor(0.8 * n)
  n_train_remaining <- n_train_target - length(fixed_indices)
  
  additional_train <- sort(sample(remaining_indices, size = n_train_remaining))
  train_indices <- sort(c(fixed_indices, additional_train))
  test_indices <- setdiff(1:n, train_indices)
  
  y_train <- y[train_indices, ]
  t_train <- t[train_indices]
  
  y_test <- y[test_indices, ]
  t_test <- t[test_indices]
  
  knots_train <- knots_quantile(t_train, dimension, order)
  B_train <- bsplines(t_train, knots_train, order)
  D_train <- bspline_jump(knots_train, order)
  fit_train <- bspline.curve.admm_lambdas(y_train, D_train, B_train, 
                                          lambdas = NULL, 
                                          lam_max = 100,
                                          lam_min = 1e-10, 
                                          n_lambda = 200, 
                                          max_iter = 1000, 
                                          epsilon = 1e-8,
                                          eta_c = 1)
  fit_train_1 = bspline.curve.admm_lambdas(y_train[,1], D_train, B_train, 
                                           lambdas = NULL, 
                                           lam_max = 100,
                                           lam_min = 1e-10, 
                                           n_lambda = 200, 
                                           max_iter = 1000, 
                                           epsilon = 1e-8,
                                           eta_c = 1)
  fit_train_2 = bspline.curve.admm_lambdas(y_train[,2], D_train, B_train, 
                                           lambdas = NULL, 
                                           lam_max = 100,
                                           lam_min = 1e-10, 
                                           n_lambda = 200, 
                                           max_iter = 1000, 
                                           epsilon = 1e-8,
                                           eta_c = 1)
  
  B_test <- bsplines(t_test, knots_train, order)
  fit_test <- B_test %*% fit_train[[which.min(fit_train$aic)]]$xi
  fit_test_1 = B_test %*% fit_train_1[[which.min(fit_train_1$aic)]]$xi
  fit_test_2 = B_test %*% fit_train_2[[which.min(fit_train_2$aic)]]$xi
  
  point_errors <- sqrt(rowSums((y_test - fit_test)^2))
  
  fit_test_combined = cbind(fit_test_1,fit_test_2)
  point_errors_individual <- sqrt(rowSums((y_test - fit_test_combined)^2))
  
  mse <- mean(point_errors^2)
  mae <- mean(point_errors)
  mxdv <- max(point_errors)
  mse_ind <- mean(point_errors_individual^2)
  mae_ind <- mean(point_errors_individual)
  mxdv_ind <- max(point_errors_individual)
  
  
  results[s, ] <- c(s, mse, mae, mxdv, mse_ind, mae_ind, mxdv_ind)
  
  cat("Seed:", s, "MSE:", mse, "MAE:", mae, "MXDV:", mxdv, 
      "MSE_ind", mse_ind, "MAE_ind", mae_ind, "MXDV_ind", mxdv_ind, "\n")
}

mse_mean <- mean(results$mse)
mse_se <- sd(results$mse) / sqrt(n_sim)
mae_mean <- mean(results$mae)
mae_se <- sd(results$mae) / sqrt(n_sim)
mxdv_mean <- mean(results$mxdv)
mxdv_se <- sd(results$mxdv) / sqrt(n_sim)
mse_ind_mean = mean(results$mse_ind)
mse_ind_se = sd(results$mse_ind) / sqrt(n_sim)
mae_ind_mean <- mean(results$mae_ind)
mae_ind_se   <- sd(results$mae_ind) / sqrt(n_sim)
mxdv_ind_mean <- mean(results$mxdv_ind)
mxdv_ind_se   <- sd(results$mxdv_ind) / sqrt(n_sim)

cat("Overall MSE: Mean =", mse_mean, ", SE =", mse_se, "\n")
cat("Overall MAE: Mean =", mae_mean, ", SE =", mae_se, "\n")
cat("Overall MXDV: Mean =", mxdv_mean, ", SE =", mxdv_se, "\n")
cat("Overall MSE_ind: Mean =", mse_ind_mean, ", SE =", mse_ind_se, "\n")
cat("Overall MAE_ind: Mean =", mae_ind_mean, ", SE =", mae_ind_se, "\n")
cat("Overall MXDV_ind: Mean =", mxdv_ind_mean, ", SE =", mxdv_ind_se, "\n")

#==============================================================================
# [Table1] Others_MSE,MAE,MXDV
#==============================================================================
n_sim <- 50
results_others <- data.frame(seed = integer(n_sim),
                             ss_mse = numeric(n_sim), ss_mae = numeric(n_sim), ss_mxdv = numeric(n_sim),  # smoothspline
                             gam_mse = numeric(n_sim), gam_mae = numeric(n_sim), gam_mxdv = numeric(n_sim),  # gam
                             ks_mse = numeric(n_sim), ks_mae = numeric(n_sim), ks_mxdv = numeric(n_sim)     # k-smooth
)

for(s in 1:n_sim){
  set.seed(s)
  
  fixed_indices <- c(1, n)
  remaining_indices <- 2:(n-1)
  
  n_train_target <- floor(0.8 * n)
  n_train_remaining <- n_train_target - length(fixed_indices)
  
  additional_train <- sort(sample(remaining_indices, size = n_train_remaining))
  train_indices <- sort(c(fixed_indices, additional_train))
  test_indices <- setdiff(1:n, train_indices)
  
  y_train <- y[train_indices, ]
  t_train <- t[train_indices]
  
  y_test <- y[test_indices, ]
  t_test <- t[test_indices]
  
  
  ############################# smooth spline : ss ###########################
  ss_lon <- smooth.spline(t_train, y_train[,1])
  ss_lat <- smooth.spline(t_train, y_train[,2])
  
  ss_pred_lon <- predict(ss_lon, x = t_test)
  ss_pred_lat <- predict(ss_lat, x = t_test)
  
  
  ss_pred_combined <- cbind(ss_pred_lon$y, ss_pred_lat$y)
  
  ss_point_errors <- sqrt(rowSums((y_test - ss_pred_combined)^2))
  ss_mse <- mean(ss_point_errors^2)
  ss_mae <- mean(ss_point_errors)
  ss_mxdv <- max(ss_point_errors)
  
  ################################### gam #####################################
  
  train_data <- data.frame(
    t = t_train,
    lon = y_train[,1],
    lat = y_train[,2]
  )
  test_data <- data.frame(t = t_test)
  
  k_candidates <- c(5, 7, 9, 11, 13, 15)
  n_train <- nrow(train_data)
  
  bic_values_lon <- numeric(length(k_candidates))
  bic_values_lat <- numeric(length(k_candidates))
  models_lon <- list()
  models_lat <- list()
  
  for(i in seq_along(k_candidates)) {
    k_val <- k_candidates[i]
    model_lon <- gam(lon ~ s(t, k = k_val), data = train_data, method = "REML")
    model_lat <- gam(lat ~ s(t, k = k_val), data = train_data, method = "REML")
    
    models_lon[[i]] <- model_lon
    models_lat[[i]] <- model_lat
    
    bic_values_lon[i] <- AIC(model_lon, k = log(n_train))
    bic_values_lat[i] <- AIC(model_lat, k = log(n_train))
    
    cat("k =", k_val, "-> BIC (lon):", bic_values_lon[i], ", BIC (lat):", bic_values_lat[i], "\n")
  }
  
  best_index_lon <- which.min(bic_values_lon)
  best_index_lat <- which.min(bic_values_lat)
  
  best_model_lon <- models_lon[[best_index_lon]]
  best_model_lat <- models_lat[[best_index_lat]]
  
  cat("optimal_lon: k =", k_candidates[best_index_lon], "BIC:", bic_values_lon[best_index_lon], "\n")
  cat("optimal_lat: k =", k_candidates[best_index_lat], "BIC:", bic_values_lat[best_index_lat], "\n")
  
  pred_gam_lon <- predict(best_model_lon, newdata = test_data)
  pred_gam_lat <- predict(best_model_lat, newdata = test_data)
  pred_gam <- cbind(pred_gam_lon, pred_gam_lat)
  
  gam_point_errors <- sqrt(rowSums((y_test - pred_gam)^2))
  
  gam_mse <- mean(gam_point_errors^2)
  gam_mae <- mean(gam_point_errors)
  gam_mxdv <- max(gam_point_errors)
  
  
  
  ############################### k-smooth : ks ###############################
  
  ks_lon <- ksmooth(t_train, y_train[,1], kernel = "normal", bandwidth = 0.07, x.points = t_test)
  ks_lat <- ksmooth(t_train, y_train[,2], kernel = "normal", bandwidth = 0.07, x.points = t_test)
  
  pred_ksmooth <- cbind(ks_lon$y, ks_lat$y)
  
  ks_point_errors <- sqrt(rowSums((y_test - pred_ksmooth)^2))
  
  ks_mse = mean(ks_point_errors^2)
  ks_mae = mean(ks_point_errors)
  ks_mxdv = max(ks_point_errors)
  
  
  results_others[s, ] <- c(s, ss_mse, ss_mae, ss_mxdv,
                           gam_mse, gam_mae, gam_mxdv,
                           ks_mse, ks_mae, ks_mxdv)
  
  cat("Seed:", s, "ss_MSE:", ss_mse, "ss_MAE:", ss_mae, "ss_MXDV:", ss_mxdv,
      "gam_MSE:", gam_mse, "gam_MAE:", gam_mae, "gam_MXDV:", gam_mxdv,
      "ks_MSE:", ks_mse, "ks_MAE:", ks_mae, "ks_MXDV:", ks_mxdv) 
}


ss_mse_mean <- mean(results_others$ss_mse)
ss_mse_se    <- sd(results_others$ss_mse) / sqrt(n_sim)
ss_mae_mean <- mean(results_others$ss_mae)
ss_mae_se    <- sd(results_others$ss_mae) / sqrt(n_sim)
ss_mxdv_mean <- mean(results_others$ss_mxdv)
ss_mxdv_se    <- sd(results_others$ss_mxdv) / sqrt(n_sim)

# GAM
gam_mse_mean <- mean(results_others$gam_mse)
gam_mse_se    <- sd(results_others$gam_mse) / sqrt(n_sim)
gam_mae_mean <- mean(results_others$gam_mae)
gam_mae_se    <- sd(results_others$gam_mae) / sqrt(n_sim)
gam_mxdv_mean <- mean(results_others$gam_mxdv)
gam_mxdv_se    <- sd(results_others$gam_mxdv) / sqrt(n_sim)


ks_mse_mean <- mean(results_others$ks_mse)
ks_mse_se    <- sd(results_others$ks_mse) / sqrt(n_sim)
ks_mae_mean <- mean(results_others$ks_mae)
ks_mae_se    <- sd(results_others$ks_mae) / sqrt(n_sim)
ks_mxdv_mean <- mean(results_others$ks_mxdv)
ks_mxdv_se    <- sd(results_others$ks_mxdv) / sqrt(n_sim)

cat("Overall MSE (ss): Mean =", ss_mse_mean, ", SE =", ss_mse_se, "\n")
cat("Overall MAE (ss): Mean =", ss_mae_mean, ", SE =", ss_mae_se, "\n")
cat("Overall MXDV (ss): Mean =", ss_mxdv_mean, ", SE =", ss_mxdv_se, "\n")
cat("Overall MSE (gam): Mean =", gam_mse_mean, ", SE =", gam_mse_se, "\n")
cat("Overall MAE (gam): Mean =", gam_mae_mean, ", SE =", gam_mae_se, "\n")
cat("Overall MXDV (gam): Mean =", gam_mxdv_mean, ", SE =", gam_mxdv_se, "\n")
cat("Overall MSE (ks): Mean =", ks_mse_mean, ", SE =", ks_mse_se, "\n")
cat("Overall MAE (ks): Mean =", ks_mae_mean, ", SE =", ks_mae_se, "\n")
cat("Overall MXDV (ks): Mean =", ks_mxdv_mean, ", SE =", ks_mxdv_se, "\n")