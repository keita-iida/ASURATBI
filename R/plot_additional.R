#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Plot elastic trees computed by MERLoT.
#'
#' This function shows elastic trees in low-dimensional spaces.
#'
#' @param ElasticTree An elastic tree object computed by MERLoT.
#' @param labels NULL or a vector of labels of all the samples,
#'   corresponding to colors.
#' @param colors NULL or a vector of colors of all the samples,
#'   corresponding to labels.
#' @param theta Angle of the plot.
#' @param phi Angle of the plot.
#' @param title Title.
#' @param xlabel x-axis label.
#' @param ylabel y-axis label.
#' @param zlabel z-axis label.
#'
#' @return A scatter3D object in plot3D package.
#' @export
#'
plot_elasticTree3D <- function(
  ElasticTree, labels, colors, theta = 30, phi = 30, title = "", xlabel = "",
  ylabel = "", zlabel = ""
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  n_groups <- length(unique(labels))
  df <- as.data.frame(ElasticTree$CellCoords)
  df$branch <- ElasticTree$Cells2Branches
  df$label <- labels
  df$color <- colors
  df <- df[order(df$label), ]
  dg <- as.data.frame(ElasticTree$Nodes)
  dh <- as.data.frame(ElasticTree$EndpointsCoords)
  di <- ElasticTree$BranchpointsCoords
  if(is.null(dim(di))){
    di <- as.data.frame(t(ElasticTree$BranchpointsCoords))
  }else{
    di <- as.data.frame(ElasticTree$BranchpointsCoords)
  }
  dj <- list()
  for(j in seq_len(n_groups)){
    dj[[j]] <- as.data.frame(dg[df$branch[[j]], ])
  }
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  mycolors <- unique(df[order(df$branch), ]$color)
  plot3D::scatter3D(df[, 1], df[, 2], df[, 3], main = title,
                    xlab = xlabel, ylab = ylabel, zlab = zlabel,
                    box = TRUE, bty = "b2", axes = TRUE, nticks = 5,
                    theta = theta, phi = phi, pch = 16, cex = 0.5, alpha = 0.1,
                    col = mycolors[df$branch], colvar = NA, colkey = FALSE)
  for(j in seq_len(n_groups)){
    plot3D::scatter3D(dj[[j]][, 1], dj[[j]][, 2], dj[[j]][, 3],
                      lwd = 2, alpha = 1, col = "black", add = TRUE)
  }
  plot3D::scatter3D(dh[, 1], dh[, 2], dh[, 3], pch = 16, cex = 3, alpha = 0.25,
                    col = "black", add = TRUE)
  plot3D::scatter3D(di[, 1], di[, 2], di[, 3], pch = 16, cex = 3, alpha = 0.40,
                    col = "black", add = TRUE)
  plot3D::scatter3D(dg[, 1], dg[, 2], dg[, 3], pch = 16, cex = 1, alpha = 1,
                    col = "black", add = TRUE)
  graphics::legend("bottomright", legend = unique(sort(df$branch)), pch = 16,
                   col = mycolors, cex = 1.3, inset = c(0.02))
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Plot pseudotime computed by MERLoT.
#'
#' This function shows pseudotimes in low-dimensional spaces.
#'
#' @param ElasticTree An elastic tree object computed by MERLoT.
#' @param Pseudotimes A pseudotime (high-dim) object computed by MERLoT.
#' @param labels NULL or a vector of labels of all the samples,
#'   corresponding to colors.
#' @param colors NULL or a vector of colors of all the samples,
#'   corresponding to labels.
#' @param theta Angle of the plot.
#' @param phi Angle of the plot.
#' @param title Title.
#' @param xlabel x-axis label.
#' @param ylabel y-axis label.
#' @param zlabel z-axis label.
#'
#' @return A scatter3D object in plot3D package.
#' @export
#'
plot_pseudotime3D <- function(
  ElasticTree, Pseudotimes, labels, colors, theta = 30, phi = 30, title = "",
  xlabel = "", ylabel = "", zlabel = ""
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  df <- as.data.frame(ElasticTree$CellCoords)
  df$pseudotime <- Pseudotimes$Times_cells
  df$pseudotime <- (df$pseudotime - min(df$pseudotime)) /
    (max(df$pseudotime) - min(df$pseudotime))
  #--------------------------------------------------
  # Definition of color
  #--------------------------------------------------
  legend.col <- function(col, lev){
    opar <- par ; n <- length(col) ; bx <- par("usr")
    box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
                bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 40)
    box.cy <- c(bx[3], bx[3]) ; box.sy <- (bx[4] - bx[3]) / n
    xx <- rep(box.cx, each = 2) ; par(xpd = TRUE)
    for(i in seq_len(n)){
      yy <- c(box.cy[1] + (box.sy * (i - 1)), box.cy[1] + (box.sy * (i)),
              box.cy[1] + (box.sy * (i)), box.cy[1] + (box.sy * (i - 1)))
      polygon(xx, yy, col = col[i], border = col[i])
    }
    par(new = TRUE)
    plot(0, 0, type = "n", ylim = c(min(lev), max(lev)),
         yaxt = "n", ylab = "", xaxt = "n", xlab = "", frame.plot = FALSE)
    axis(side = 4, las = 2, tick = FALSE, line = .25)
    par <- opar
  }
  mycolor <- viridis::plasma(50)[as.numeric(cut(df$pseudotime, breaks = 50))]
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  par(mar = c(2.0, 2.0, 2.0, 2.0))
  par(oma = c(2.0, 2.0, 2.0, 2.0))
  plot3D::scatter3D(df[,1], df[,2], df[,3], main = title, xlab = "DC_1",
                    ylab = "DC_2", zlab = "DC_3", box = TRUE, bty = "b2",
                    axes = TRUE, nticks = 5, theta = theta, phi = phi, pch = 16,
                    cex = 1.0, alpha = 0.8, col = mycolor, colvar = NA,
                    colkey = FALSE)
  legend.col(col = viridis::plasma(50), lev = df$pseudotime)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Plot pseudotime courses of sign scores with elastic tree.
#'
#' This function shows sign scores and elastic trees along pseudotimes.
#'
#' @param sce A SingleCellExperiment object.
#' @param signID A signID.
#' @param ElasticTree An elastic tree object computed by MERLoT.
#' @param SignsSpaceEmbedding A sign space embedding object computed by MERLoT.
#' @param Pseudotimes A pseudotime (high-dim) object computed by MERLoT.
#' @param labels NULL or a vector of labels of all the samples.
#' @param range_y "cells" or "tree", which determines the range of y-axis.
#'
#' @return A ggplot object.
#' @export
#'
plot_pseudotimecourse_wTree <- function(
  sce, signID, ElasticTree, SignsSpaceEmbedding, Pseudotimes, labels, range_y
){
  mat <- as.matrix(assay(sce, "counts"))
  #--------------------------------------------------
  # Prepare data frames to be output
  #--------------------------------------------------
  Sign <- which(colnames(SignsSpaceEmbedding$CellCoords) == signID)
  SignScore = SignsSpaceEmbedding$Nodes[, Sign]
  SignScoreOriginal = SignsSpaceEmbedding$CellCoords[, Sign]

  if(range_y == "cells"){
    range_y = range(SignScoreOriginal)
  }else if(range_y == "tree"){
    range_y = range(SignScore)
  }

  ptc <- Pseudotimes$Times_cells
  Pseudotimes$Times_cells <- (ptc - min(ptc)) / (max(ptc) - min(ptc))
  pty <- Pseudotimes$Times_yk
  Pseudotimes$Times_yk <- (pty - min(pty)) / (max(pty) - min(pty))
  pty <- Pseudotimes$Times_yk

  # cells
  df <- numeric(0)
  df <- data.frame(pseudotime = Pseudotimes$Times_cells,
                   signscore = SignScoreOriginal,
                   color = labels)

  # branch and end points
  dg <- numeric(0)
  dg <- data.frame(pseudotime = pty[SignsSpaceEmbedding$Topology$Endpoints],
    signscore = SignScore[SignsSpaceEmbedding$Topology$Endpoints],
    color = 0)
  if(length(SignsSpaceEmbedding$Topology$Endpoints) > 2){
    add <- data.frame(
      pseudotime = pty[SignsSpaceEmbedding$Topology$Branchpoints],
      signscore = SignScore[SignsSpaceEmbedding$Topology$Branchpoints],
      color = 0)
    dg <- rbind(dg, add)
  }

  # branches
  dh <- list()
  for(i in seq_along(Pseudotimes$Branches)){
    # If and endpoint or branchpoint was selected as t0
    if(is.null(Pseudotimes$C0)){
      branch_i <- Pseudotimes$Branches[[i]]
      aux <- sort(pty[branch_i], index.return = T)
      dh[[i]] <- data.frame(x = aux$x,
        y = SignsSpaceEmbedding$Nodes[branch_i[aux$ix], Sign],
        color = i)
    }else{
      stop("Set t0 at end point or branch point")
    }
  }
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  mycol <- as.factor(mycolor[dh[[i]]$color])
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = df[, 1], y = df[, 2],
                                     color = as.factor(df$color)),
                        size = 2, alpha = 0.1) +
    ggplot2::geom_point(ggplot2::aes(x = dg[,1], y = dg[,2]),
                        color = "black", size = 5, alpha = 1)
  for(i in seq_along(dh)){
    p <- p + ggplot2::geom_line(
      ggplot2::aes_string(x = dh[[i]][, 1], y = dh[[i]][, 2],
                          color = as.factor(dh[[i]]$color)),
      size = 2, alpha = 1)
  }

  return(p)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Plot pseudotime courses of sign scores without elastic tree.
#'
#' This function shows sign scores along pseudotimes.
#'
#' @param sce A SingleCellExperiment object.
#' @param signID A signID.
#' @param Pseudotimes A pseudotime (high-dim) object computed by MERLoT.
#' @param range_y "cells" or "tree", which determines the range of y-axis.
#'
#' @return A ggplot object.
#' @export
#'
plot_pseudotimecourse_woTree <- function(sce, signID, Pseudotimes, range_y){
  df <- data.frame(x = Pseudotimes$Times_cells)
  df$pseudotime <- (df$x - min(df$x)) / (max(df$x) - min(df$x))
  df$branch <- Pseudotimes$Cells2Branches
  n_branch <- length(unique(df$branch))
  #--------------------------------------------------
  # Prepare a submatrix for the input sign.
  #--------------------------------------------------
  ind <- which(rownames(sce) == signID)
  subsce <- sce[ind, ]
  df$score <- t(as.matrix(assay(subsce, "counts")))
  dg <- numeric(0)
  for(j in seq_len(n_branch)){
    tmp <- df[which(df$branch == j),]
    pt <- unique(sort(tmp$pseudotime))
    for(k in seq_along(pt)){
      dg <- rbind(dg,
        c(pt[k], j, mean(tmp[which(tmp$pseudotime == pt[k]), ]$score),
          ifelse(length(tmp[which(tmp$pseudotime == pt[k]), ]$score) >= 3,
                 sd(tmp[which(tmp$pseudotime == pt[k]), ]$score), 0)))
    }
  }
  dg <- as.data.frame(dg)
  colnames(dg) <- c("pseudotime", "branch", "mean", "sd")
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = df$pseudotime, y = df$score,
                                     color = as.factor(df$branch)),
                        size = 2, alpha = 0.1) +
    ggplot2::geom_line(ggplot2::aes(x = dg$pseudotime, y = dg$mean,
                                    color = as.factor(dg$branch)),
                       size = 2, alpha = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(x = dg$pseudotime, y = dg$mean,
                                      ymin = dg$mean - dg$sd,
                                      ymax = dg$mean + dg$sd,
                                      color = as.factor(dg$branch),
                                      fill = as.factor(dg$branch)),
                         size = 0, alpha = 0.3)

  return(p)
}

