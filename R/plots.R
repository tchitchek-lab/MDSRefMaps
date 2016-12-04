#' @title Generation of ggplot graphical representations for MDS Reference Maps and MDS Projections.
#'
#' @description
#' This function generates a graphical representation for MDS Reference Maps and MDS Projections.
#' Objects can be shaped and colored using the 'color' and 'shape' parameters.
#' Convex hulls for given sets of objects can also be displayed.
#'
#' @param mds a MDS result provided by the 'MDSReferenceMap()', 'MDSProjection()', or 'MDSReferenceMapwithProjection()' functions
#' @param color a vector used to color the points
#' @param shape a vector used for the shape of the points
#' @param polygon a vector used to color the convex hull
#' @param title a character specifying the title of the representation
#' @param display.legend a logical specifying if the graphic legend must be displayed
#'
#' @export 
plotMDS = function(mds, color=NULL, shape=NULL, polygon=NULL, title = "MDS", display.legend = TRUE){
	
	if (is.null(mds$points)) {stop(
        "mds$points must be a matrix of a MDS representation")}
    if (!is.null(color) && (length(color) != nrow(mds$points))) {
        stop("color must be the same length as the number of objects in the mds")}
    if (!is.null(shape) && (length(shape) != nrow(mds$points))) {
        stop("shape must be the same length as the number of objects in the mds")}
    if (!is.null(polygon) && (length(polygon) != nrow(mds$points))) {
        stop("polygon must be the same length as the number of objects in the mds")}

    subtitle  = paste("Kruskal Stress =",format(mds$stress,digits = 4),"; Entourage Score =",format(mds$entourage,digits = 4))

    idcolor   = 0
    idshape   = 0
    idpolygon = 0

    shapes = list(
		"square",
		"circle",
		"triangle",
		"plus",
		"cross",
		"diamond",
		"trianglepointdown",
		"squarecross",
		"star",
		"diamondplus",
		"circleplus",
		"trianglesupanddown",
		"squareplus",
		"circlecross",
		"squareandtriangledown",
		"filledsquare",
		"filledcircle",
		"filledtrianglepointup",
		"filleddiamond",
		"solidcircle",
		"bullet",
		"filledcircleblue",
		"filledsquareblue",
		"filleddiamondblue",
		"filledtrianglepointupblue",
		"filledtrianglepointdownblue")

    dicShapes = list(
		"square"=0,
		"circle"=1,
		"triangle"=2,
		"plus"=3,
		"cross"=4,
		"diamond"=5,
		"trianglepointdown"=6,
		"squarecross"=7,
		"star"=8,
		"diamondplus"=9,
		"circleplus"=10,
		"trianglesupanddown"=11,
		"squareplus"=12,
		"circlecross"=13,
		"squareandtriangledown"=14,
		"filledsquare"=15,
		"filledcircle"=16,
		"filledtrianglepointup"=17,
		"filleddiamond"=18,
		"solidcircle"=19,
		"bullet"=20,
		"filledcircleblue"=21,
		"filledsquareblue"=22,
		"filleddiamondblue"=23,
		"filledtrianglepointupblue"=24,
		"filledtrianglepointdownblue"=25)
    
	nref = NULL
	split = "_"
	
    splitter = function(x, max){
        if (nchar(x) > max){
            paste0(substr(x,1,max),"\n",substr(x,max+1,nchar(x)))
        }else{
            x
        }
    }
    
    points = data.frame(x = mds$points[,1], y = mds$points[,2])
    
    if (!is.null(color)) {
        points$color = color
        if (all(points$color %in% grDevices::colors())) {
            idcolor  = 1
        }else if (is.numeric(split)) {
            lmax = max(nchar(points$color))
            lseuil = split
            while (lmax > lseuil) {
                points$color = sapply(points$color, splitter, lseuil)
                lseuil = lseuil + split
            }
        }else {
            points$color = gsub(split,"\n", points$color)
        }
    }else {
        idcolor      = 1
        points$color = rep("black",length(points$x))
        if (!is.null(nref)) {points$color = c(rep("grey", nref), rep("red", (length(points$x) - nref)))}
    }

    if (!is.null(shape)) {
        points$shape = shape
        if (is.numeric(points$shape)) {
            points$shape = as.factor(cut(points$shape, breaks = 5))
        }else if (any(points$shape %in% shapes)) {
            points$shape = unlist(sapply(points$shape,function(x)dicShapes[[x]]))
            idshape  = 1
            if (length(points$shape) != nrow(mds$points)) {
                stop(paste("A specific shape is used, please ensure to only use recognised shapes :",paste(shapes, collapse = ', ')))
            }
        }else if (is.numeric(split)) {
            lmax = max(nchar(points$shape))
            lseuil = split
            while (lmax > lseuil) {
                points$shape = sapply(points$shape, splitter, lseuil)
                lseuil = lseuil + split
            }
        }else {
            points$shape = gsub(split,"\n", points$shape)
        }
    }else {
        idshape      = 1
        points$shape = rep(19,length(points$x))
        if (!is.null(nref)) {points$shape = c(rep(17, nref), rep(19, (length(points$x) - nref)))}
    }
    
    if (!is.null(polygon)) {
        points$polygon = polygon
        if (all(points$polygon %in% grDevices::colors())) {
            idpolygon = 1
        }else if (is.numeric(split)) {
            lmax = max(nchar(points$polygon))
            lseuil = split
            while (lmax > lseuil) {
                points$polygon = sapply(points$polygon, splitter, lseuil)
                lseuil = lseuil + split
            }
        }else {
            points$polygon = gsub(split,"\n", points$polygon)
        }
        hulls          = plyr::ddply(points, "polygon", function(df){df[grDevices::chull(df$x, df$y), ]})
        hulls$polygoncolor = hulls$polygon
        if (!is.null(nref)) {
            idpolygon          = 1
            hulls$polygoncolor = hulls$color
            hullsproj          = hulls[which(hulls$polygoncolor == "red"),]
            hulls              = hulls[which(hulls$polygoncolor == "grey"),]
        }
    }
    
    if (!is.null(nref)) {
        points$color = c(rep("Reference objects", nref), 
                         rep("Projected objects", (length(points$x) - nref)))
        if (!is.null(points$polygon)) {
            tmp = points$polygon[nref+1]
            points$color = c(rep("Reference objects", nref), 
                             rep(tmp, (length(points$x) - nref)))
        }
        points$color = as.factor(points$color)
        points$color = factor(points$color, levels = rev(levels(points$color)))
    }
    
	title = bquote(atop(.(title), atop(italic(.(subtitle)), "")))
	plot = ggplot2::ggplot(points, ggplot2::aes_string(x = "x", y = "y")) +
           ggplot2::ggtitle(title)
    
    plot = plot + ggplot2::geom_point(ggplot2::aes_string(color = "color",shape = "shape")) +
        ggplot2::labs(colour = names(color), shape = names(shape))
    
    if (!is.null(polygon)) {
        plot = plot + ggplot2::geom_polygon(data = hulls, 
                      ggplot2::aes_string(x = "x", y = "y", group = "polygon", fill = "polygoncolor"), alpha = 0.6) +
                      ggplot2::labs(fill = names(polygon))
        if (!is.null(nref)) {
            plot = plot + ggplot2::geom_polygon(data = hullsproj, 
                          ggplot2::aes_string(x = "x", y = "y", group = "polygon", fill = "polygoncolor"), alpha = 0.6) +
                          ggplot2::labs(fill = names(polygon))
        }
        if (idpolygon) {
            plot = plot + ggplot2::scale_fill_identity()
        }
    }
    
    if (idcolor && is.null(nref)) {
        values=levels(as.factor(points$color))
        plot = plot + ggplot2::scale_colour_manual(values=values)
    }
    
    if (!is.null(nref)) {
        values = c("grey", "red")
        if (idcolor && !is.null(color)) {values = c("grey", color[nref+1])}
        names(values) = c("Reference objects", "Projected objects")
        if (!is.null(points$polygon)) {
            tmp = points$polygon[nref+1]
            names(values) = c("Reference objects", tmp)
        }
        plot = plot + ggplot2::scale_colour_manual(values=values)
    }
    
    if (idshape) {
        plot = plot + ggplot2::scale_shape_identity()
    }
    
    plot = plot + 
        ggplot2::geom_hline(yintercept = (min(points$y) + max(points$y))/2, linetype = "dashed", color = "black", size = 0.5) + 
        ggplot2::geom_vline(xintercept = (min(points$x) + max(points$x))/2, linetype = "dashed", color = "black", size = 0.5)
    
    plot = plot + ggplot2::scale_y_continuous(expand = c(0,0)) +
                  ggplot2::scale_x_continuous(expand = c(0,0)) +
                  ggplot2::xlab("dimension 1") +
                  ggplot2::ylab("dimension 2")
                  
	if(display.legend==TRUE){
		plot = plot + ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'white', colour = 'black'),
			   panel.border     = ggplot2::element_rect(fill = NA),
			   panel.grid.minor = ggplot2::element_blank(),
			   panel.grid.major = ggplot2::element_blank(),
			   axis.ticks       = ggplot2::element_blank(),
			   axis.text.x      = ggplot2::element_blank(),
			   axis.text.y      = ggplot2::element_blank(),
			   legend.text      = ggplot2::element_text(size = 8),
			   aspect.ratio     = 1,
			   plot.title       = ggplot2::element_text(hjust = 0.5))
	}else{	
		plot = plot + ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'white', colour = 'black'),
				   panel.border     = ggplot2::element_rect(fill = NA),
				   panel.grid.minor = ggplot2::element_blank(),
				   panel.grid.major = ggplot2::element_blank(),
				   axis.ticks       = ggplot2::element_blank(),
				   axis.text.x      = ggplot2::element_blank(),
				   axis.text.y      = ggplot2::element_blank(),
				   legend.text      = ggplot2::element_text(size = 8),
				   aspect.ratio     = 1,
				   legend.position  = "none",
				   plot.title       = ggplot2::element_text(hjust = 0.5))
	}
	
    return(plot)
}
