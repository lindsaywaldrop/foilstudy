# Custom R functions

calc_camber <- function(x, y, check = NULL){
  camber <- (max(y) - min(y)) / (max(x) - min(x))
  if(!is.null(check)){
    return(all.equal(camber, check))
  }else{
    return(camber) 
  }
}

mag <- function(x) {
  return(sqrt((x[1])^2 + (x[2])^2))
}

dist <- function(x1, x2, y1, y2) {
  return(sqrt((x2 - x1)^2 + (y2 - y1)^2))
}

calc_angle <- function(p1, p2, p3){
  a <- c(p3[1] - p2[1], p3[2] - p2[2])
  b <- c(p1[1] - p2[1], p1[2] - p2[2])
  alpha <- acos((a[1]*b[1] + a[2]*b[2]) / (mag(a) * mag(b)))
  return(alpha)
}

correct_2D_pts <- function(pts, center_pt, alpha){
  rotm <- matrix(c(cos(alpha), sin(alpha), -sin(alpha), cos(alpha)), ncol = 2)
  if(!is.matrix(pts)) pts <- as.matrix(pts)
  rot_pts <- t(rotm %*% ( t(pts) - c(center_pt[1], center_pt[2])) + 
                 c(center_pt[1], center_pt[2]))
  corr_pts <- rot_pts
  corr_pts[,1] <- rot_pts[,1] - (center_pt[1])
  corr_pts[,2] <- rot_pts[,2] - (center_pt[2])
  return(corr_pts)
}

find_midline <- function(x, y, smoothed = TRUE, plot = TRUE){
  npts <- length(x)
  foil_midline_x <- rep(NA, npts/2)
  foil_midline_y <- rep(NA, npts/2)
  thickness_top_x <- rep(NA, npts/2)
  thickness_top_y <- rep(NA, npts/2)
  thickness_bot_x <- rep(NA, npts/2)
  thickness_bot_y <- rep(NA, npts/2)
  for(i in 1:(npts/2)){
    foil_midline_x[i] <- (x[i] - x[(npts + 1) - i])/2 + 
      x[(npts + 1) - i]
    foil_midline_y[i] <- (y[i] - y[(npts + 1) - i])/2 + 
      y[(npts + 1) - i]
  }
  for(i in 1:(npts/2)){
    thickness_top_x[i] <- x[i] - foil_midline_x[i]
    thickness_top_y[i] <- y[i] - foil_midline_y[i]
    thickness_bot_x[i] <- foil_midline_x[i] - x[(npts + 1) - i]
    thickness_bot_y[i] <- foil_midline_y[i] - y[(npts + 1) - i]
  } 
  if(smoothed) {
    loess_sm <- loess(foil_midline_y ~ foil_midline_x, span = 0.5, parametric = "span")
    foil_midline_y <- predict(loess_sm, foil_midline_x)
  }
  
  if(plot){
    plot(x, y, type = "o", col = viridis::viridis(npts), asp = 1)
     #xlim = c(0.9, 1), ylim = c(0, 0.03))
    lines(x, y)
    lines(c(-10, 10), c(0, 0), lty = 2)
    points(foil_midline_x, foil_midline_y, col = viridis::viridis(npts/2, option = "A"), 
           pch = 19)
    for(i in 1:(npts/2)) {
      lines(x = c(x[i], x[(npts + 1) - i]), 
            y = c(y[i], y[(npts + 1) - i]))
      lines(x = c(foil_midline_x[i], foil_midline_x[i] + thickness_top_x[i]), 
            y = c(foil_midline_y[i], foil_midline_y[i] + thickness_top_y[i]), 
            col = "red")
    }
    for(i in (npts/2):1){
      lines(x = c(foil_midline_x[i], foil_midline_x[i] - thickness_top_x[i]), 
            y = c(foil_midline_y[i], foil_midline_y[i] - thickness_top_y[i]), 
            col = "blue")
    }
  }
  foil_midline <- data.frame("x" = foil_midline_x, "y" = foil_midline_y, 
                             "thickness_top_x" = thickness_top_x, 
                             "thickness_top_y" = thickness_top_y, 
                             "thickness_bot_x" = thickness_bot_x, 
                             "thickness_bot_y" = thickness_bot_y)
  return(foil_midline)
}

adjust_midline <- function(foil_midline, camber_new, plot = FALSE){
  new_midline <- foil_midline
  new_midline$y <- ((foil_midline$y - min(foil_midline$y)) * camber_new ) /
    (max(foil_midline$y) - min(foil_midline$y))
  if(min(new_midline$y) < 0) new_midline$y <- new_midline$y + abs(min(new_midline$y))
  if(plot){
    require(viridis)
    plot(foil_midline$x, foil_midline$y, asp = 1, col = viridis(80))
    points(new_midline$x, new_midline$y)
  }
  return(new_midline)
}

calc_midline_angles <- function(new_midline){
  angles <- rep(NA, nrow(new_midline))
  for(i in 2:(nrow(new_midline))){
    angles[i] <- calc_angle(p1 = c(new_midline$x[i], new_midline$y[i]), 
                            p2 = c(new_midline$x[i-1], new_midline$y[i-1]), 
                            p3 = c(new_midline$x[i], new_midline$y[i-1]))
  }
  angles[1] <- angles[2]
  return(angles)
}

create_new_foil <- function(new_midline, plot = FALSE){
  npts <- nrow(new_midline)*2
  half_thickness <- sqrt((new_midline$thickness_top_x)^2 + (new_midline$thickness_top_y)^2)
  top_angles <- calc_midline_angles(new_midline) - 0.5*pi
  bot_angles <- calc_midline_angles(new_midline) + 0.5*pi
  new_foil <- data.frame("x" = rep(NA, npts), "y" = rep(NA, npts))
  new_foil$x[1:(npts/2)] <- new_midline$x - half_thickness*cos(top_angles)
  new_foil$y[1:(npts/2)] <- new_midline$y - half_thickness*sin(top_angles)
  new_foil$x[npts:(npts/2 + 1)] <- new_midline$x - half_thickness*cos(bot_angles)
  new_foil$y[npts:(npts/2 + 1)] <- new_midline$y - half_thickness*sin(bot_angles)
  new_foil$x[1] <- 1
  new_foil$y[1] <- 0
  new_foil$x[npts] <- new_foil$x[1] - 5e-4
  new_foil$y[npts] <- 0
  if(plot){
    require(viridis)
    plot(new_foil$x, new_foil$y, asp = 1)#, 
         #xlim = c(0.9, 1), ylim = c(0, 0.03))
    lines(new_foil$x, new_foil$y, col = "black")
    points(new_midline$x, new_midline$y, col = viridis(80), pch = 19)
    lines(c(0,1), c(0,0), lty = 2)
    for(i in 1:80){
      lines(x = c(new_midline$x[i], new_midline$x[i] - half_thickness[i]*cos(top_angles[i])), 
            y = c(new_midline$y[i], new_midline$y[i] - half_thickness[i]*sin(top_angles[i])), 
            col = "red")
    }
    for(i in 80:1){
      lines(x = c(new_midline$x[i], 
                  new_midline$x[i] - half_thickness[i]*cos(bot_angles[i])), 
            y = c(new_midline$y[i], 
                  new_midline$y[i] - half_thickness[i]*sin(bot_angles[i])), 
            col = "blue")
    }
  }
  return(new_foil)
}
