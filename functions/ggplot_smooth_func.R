# https://gist.github.com/kdauria/524eade46135f6348140

stat_smooth_func <- function(mapping = NULL, data = NULL,
                             geom = "smooth", position = "identity",
                             ...,
                             method = "auto",
                             formula = y ~ x,
                             se = TRUE,
                             n = 80,
                             span = 0.75,
                             fullrange = FALSE,
                             level = 0.95,
                             method.args = list(),
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE,
                             xpos = NULL,
                             ypos = NULL) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatSmoothFunc,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      method = method,
      formula = formula,
      se = se,
      n = n,
      fullrange = fullrange,
      level = level,
      na.rm = na.rm,
      method.args = method.args,
      span = span,
      xpos = xpos,
      ypos = ypos,
      ...
    )
  )
}

# get p-value using lmp ()
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

StatSmoothFunc <- ggproto("StatSmooth", Stat,
                          
                          setup_params = function(data, params) {
                            # Figure out what type of smoothing to do: loess for small datasets,
                            # gam with a cubic regression basis for large data
                            # This is based on the size of the _largest_ group.
                            if (identical(params$method, "auto")) {
                              max_group <- max(table(data$group))
                              
                              if (max_group < 1000) {
                                params$method <- "loess"
                              } else {
                                params$method <- "gam"
                                params$formula <- y ~ s(x, bs = "cs")
                              }
                            }
                            if (identical(params$method, "gam")) {
                              params$method <- mgcv::gam
                            }
                            
                            params
                          },
                          
                          compute_group = function(data, scales, method = "auto", formula = y~x,
                                                   se = TRUE, n = 80, span = 0.75, fullrange = FALSE,
                                                   xseq = NULL, level = 0.95, method.args = list(),
                                                   na.rm = FALSE, xpos=NULL, ypos=NULL) {
                            if (length(unique(data$x)) < 2) {
                              # Not enough data to perform fit
                              return(data.frame())
                            }
                            
                            if (is.null(data$weight)) data$weight <- 1
                            
                            if (is.null(xseq)) {
                              if (is.integer(data$x)) {
                                if (fullrange) {
                                  xseq <- scales$x$dimension()
                                } else {
                                  xseq <- sort(unique(data$x))
                                }
                              } else {
                                if (fullrange) {
                                  range <- scales$x$dimension()
                                } else {
                                  range <- range(data$x, na.rm = TRUE)
                                }
                                xseq <- seq(range[1], range[2], length.out = n)
                              }
                            }
                            # Special case span because it's the most commonly used model argument
                            if (identical(method, "loess")) {
                              method.args$span <- span
                            }
                            
                            if (is.character(method)) method <- match.fun(method)
                            
                            base.args <- list(quote(formula), data = quote(data), weights = quote(weight))
                            model <- do.call(method, c(base.args, method.args))
                            
                            m = model
                            #eq <- substitute(italic(y) == a + b %.% italic(x)*"\n"~~italic(r)^2~"="~r2~","~~italic(p)~"="~pval, 
                            eq <- substitute(italic(y) == a + b %.% italic(x)*~~italic(r)^2~"="~r2~","~~italic(p)~"="~pval, 
                                             list(a = format(coef(m)[[1]], digits = 2), 
                                                  b = format(coef(m)[[2]], digits = 2), 
                                                  r2 = format(summary(m)$r.squared, digits = 3),
                                                  pval = ifelse(lmp(m)<1e-4, "<0.0001",format(lmp(m), digits = 3))))
                            # dput(eq)
                            func_string = as.character(as.expression(eq))
                            # dput(func_string)
                            if(is.null(xpos)) xpos = min(data$x)*0.9
                            if(is.null(ypos)) ypos = max(data$y)*0.9
                            data.frame(x=xpos, y=ypos, label=func_string)
                            
                          },
                          
                          required_aes = c("x", "y")
)


stat_smooth_func_with_pval <- function(mapping = NULL, data = NULL,
                                            geom = "smooth", position = "identity",
                                            ...,
                                            method = "auto",
                                            formula = y ~ x,
                                            se = TRUE,
                                            n = 80,
                                            span = 0.75,
                                            fullrange = FALSE,
                                            level = 0.95,
                                            method.args = list(),
                                            na.rm = FALSE,
                                            show.legend = NA,
                                            inherit.aes = TRUE,
                                            xpos = NULL,
                                            ypos = NULL,
                                            xpos2 = NULL,
                                            ypos2 = NULL) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatSmoothFunc,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      method = method,
      formula = formula,
      se = se,
      n = n,
      fullrange = fullrange,
      level = level,
      na.rm = na.rm,
      method.args = method.args,
      span = span,
      xpos = xpos,
      ypos = ypos,
      xpos2 = xpos2,
      ypos2 = ypos2,
      ...
    )
  )
}

StatSmoothFunc <- ggproto("StatSmooth", Stat,
                          
                          setup_params = function(data, params) {
                            # Figure out what type of smoothing to do: loess for small datasets,
                            # gam with a cubic regression basis for large data
                            # This is based on the size of the _largest_ group.
                            if (identical(params$method, "auto")) {
                              max_group <- max(table(data$group))
                              
                              if (max_group < 1000) {
                                params$method <- "loess"
                              } else {
                                params$method <- "gam"
                                params$formula <- y ~ s(x, bs = "cs")
                              }
                            }
                            if (identical(params$method, "gam")) {
                              params$method <- mgcv::gam
                            }
                            
                            params
                          },
                          
                          compute_group = function(data, scales, method = "auto", formula = y~x,
                                                   se = TRUE, n = 80, span = 0.75, fullrange = FALSE,
                                                   xseq = NULL, level = 0.95, method.args = list(),
                                                   na.rm = FALSE, xpos=NULL, ypos=NULL, 
                                                   xpos2=NULL, ypos2=NULL) {
                            if (length(unique(data$x)) < 2) {
                              # Not enough data to perform fit
                              return(data.frame())
                            }
                            
                            if (is.null(data$weight)) data$weight <- 1
                            
                            if (is.null(xseq)) {
                              if (is.integer(data$x)) {
                                if (fullrange) {
                                  xseq <- scales$x$dimension()
                                } else {
                                  xseq <- sort(unique(data$x))
                                }
                              } else {
                                if (fullrange) {
                                  range <- scales$x$dimension()
                                } else {
                                  range <- range(data$x, na.rm = TRUE)
                                }
                                xseq <- seq(range[1], range[2], length.out = n)
                              }
                            }
                            # Special case span because it's the most commonly used model argument
                            if (identical(method, "loess")) {
                              method.args$span <- span
                            }
                            
                            if (is.character(method)) method <- match.fun(method)
                            
                            base.args <- list(quote(formula), data = quote(data), weights = quote(weight))
                            model <- do.call(method, c(base.args, method.args))
                            
                            m = model
                            a <- format(coef(m)[[1]], digits = 3)
                            b <- format(coef(m)[[2]], digits = 3)
                            if (b<0){
                              eq1 <- substitute(italic(y) ==  a - b %.% italic(x), 
                                                list(a = a, b = abs(b)))
                            }
                            else{
                              eq1 <- substitute(italic(y) ==  a + b %.% italic(x), 
                                              list(a = a, b = b))
                            }
                            
                            
                            eq2 <- substitute(italic(r)^2~"="~r2*","~~italic(p)~"="~pval, 
                                              list(r2 = format(summary(m)$r.squared, digits = 3),
                                                   pval = format(summary(m)$coef[2,4], digits = 3)))
                            
                            func_string1 = as.character(as.expression(eq1))
                            func_string2 = as.character(as.expression(eq2))
                            
                            if(is.null(xpos)) xpos = min(data$x)-0.1*(max(data$x)-min(data$x))
                            if(is.null(ypos)) ypos = max(data$y)+0.1*(max(data$y)-min(data$y))
                            if(is.null(xpos2)) xpos2 = xpos
                            if(is.null(ypos2)) ypos2 = max(data$y)-0.1*(max(data$y)-min(data$y))
                            
                            data.frame(x = rbind(xpos, xpos2), 
                                       y = rbind(ypos, ypos2), 
                                       label = rbind(func_string1, func_string2))
                            
                          },
                          
                          required_aes = c("x", "y")
)
