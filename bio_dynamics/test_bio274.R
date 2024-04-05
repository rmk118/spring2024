# testing code for bio 274 HW

library(tidyverse)

y0 <- matrix(c(0, 1, 0, 4, -6, 1, 5, 0.5, 0, -3), 5, 2, byrow = TRUE)

y0_df <- tibble(x=c(0,0,-6,5,0), 
                y=c(1,4,1,0.5,-3),
                IC=c("V1", "V2", "V3", "V4", "V5")) %>% 
  mutate(pt = paste0("(",x, ",", y, ")"))

fun_traj <- function(t, y, parameters) {
  x <- y[1]
  y <- y[2]
  a <- parameters[1]
  b <- parameters[2]
  c <- parameters[3]
  d <- parameters[4]
  dy <- numeric(2)
  dy[1] <- a*x + b*y
  dy[2] <- c*x + d*y  
  return(list(dy))
}

fun_plane <- c(dx ~ a*x+b*y, dy ~c*x+d*y)

make_plane <- function(params) {
  params_unnamed <- unname(params)
  traj <- trajectory(fun_traj, y0 = y0, tlim= c(0, 1), tstep = 0.005,
                     parameters = params_unnamed, add=FALSE)
  
  t <- c(traj$t)
  x <- as_tibble(traj$x) %>% bind_cols("t"=t) %>% 
    pivot_longer(cols = c(1:5), names_to = "IC", values_to = "x")
  y <- as_tibble(traj$y) %>% bind_cols("t"=t) %>% 
    pivot_longer(cols = c(1:5), names_to = "IC", values_to = "y")
  z <- left_join(x,y, by = join_by(t, IC)) %>% filter(x<7 & x>-7 & y<7 & y>-7) %>% 
    left_join(y0_df %>% select(IC, pt), by="IC")
  
  p1 <- phaseplane(fun_plane, x_var="x" ,y_var="y",
                   x_window = c(-7, 7), y_window = c(-7, 7),
                   parameters = params)
  
  A<- matrix(params_unnamed, 2, 2, byrow=TRUE)
  v1_slope <- eigen(A)$vectors[2,1] / eigen(A)$vectors[1,1]
  v2_slope <- eigen(A)$vectors[2,2] / eigen(A)$vectors[1,2]
  
  p1+
    theme_bw()+
    geom_abline(slope=v2_slope, color="gray40", linewidth=1)+
    geom_abline(slope=v1_slope, color="gray40", linewidth=1)+
    geom_path(data=z, aes(x=x, y=y, color=pt), linewidth=1.5)+
    geom_point(data=y0_df, aes(x=x,y=y, color=pt), shape=21, fill="white",size=4, stroke=2)+labs(color="IC")
}

parameters=c(a=1,b=2,c=3,d=2)
make_plane(parameters)

make_plane_complex(c(a=-3,b=10,c=-1,d=3))

A <- matrix(c(-3, 10, -1, 3), 2, 2, byrow=TRUE)


flowField(fun_traj,  xlim = c(-50, 50),
          ylim       = c(-40, 40),
          parameters = c(-3, 10, -1, 3),
          points     = 30,
          add        = FALSE)



# April 5 -----------------------------------------------------------------

# Phase Portrait setup
y0 <- matrix(c(5, 0), 1, 2, byrow = TRUE)
# 
# y0_df <- tibble(x=c(0), 
#                 y=c(4),
#                 IC=c("V1")) %>% 
#   mutate(pt = paste0("(",x, ",", y, ")"))

fun_traj <- function(t, y, parameters) {
  x <- y[1]
  y <- y[2]
  a <- parameters[1]
  b <- parameters[2]
  c <- parameters[3]
  d <- parameters[4]
  dy <- numeric(2)
  dy[1] <- a*x + b*y
  dy[2] <- c*x + d*y  
  return(list(dy))
}

fun_plane <- c(dx ~ a*x+b*y, dy ~c*x+d*y)


  params_unnamed <- unname(c(a=-0.25,b=-0.5,c=0.25,d=0))
  traj <- trajectory(fun_traj, y0 = y0, tlim= c(0, 555), tstep = 0.1,
                     parameters = params_unnamed, add=FALSE)
  
  t <- c(traj$t)
  x <- as_tibble(traj$x) %>% rename("x"="V1") %>% bind_cols("t"=t)
  y <- as_tibble(traj$y) %>% rename("y"="V1")%>% bind_cols("t"=t)
  z <- left_join(x,y, by = join_by(t)) %>% filter(x<7 & x>-7 & y<7 & y>-7)
  
  p1 <- phaseplane(fun_plane, x_var="x" ,y_var="y",
                   x_window = c(-7, 7), y_window = c(-7, 7),
                   parameters = c(a=-0.25,b=-0.5,c=0.25,d=0))
  
  p1+theme_bw()+
    geom_path(data=z, aes(x=x, y=y, group=pt), linewidth=1.5)

