# From: https://github.com/tavisrudd/r_users_group_1/blob/master/calendarHeat.R

##############################################################################
 # Calendar Heatmap #
 # by #
 # Paul Bleicher #
 # an R version of a graphic from: #
 # http://stat-computing.org/dataexpo/2009/posters/wicklin-allison.pdf #
 # requires lattice, chron, grid packages #
 ##############################################################################

## calendarHeat: An R function to display time-series data as a calendar heatmap
## Copyright 2009 Humedica. All rights reserved.

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.

## You can find a copy of the GNU General Public License, Version 2 at:
## http://www.gnu.org/licenses/gpl-2.0.html

################################################################################
# calendarHeat: R function to display time-series data as a calendar heatmap   #
################################################################################
# Author of this slightly modified verison: Mauricio Zambrano-Bigiarini        #
################################################################################
# Started: 2009 ;                                                              #
# Updates: 25-Nov-2023 ; 27-Nov-2023 ; 29-Nov-2023 ; 30-Nov-2023               #
#          17-Jan-2024                                                         #
################################################################################

calendarHeatmap <- function(x, ...) UseMethod("calendarHeatmap")

calendarHeatmap.zoo <- function(x,
                            #ncolors=5,
                            #color="r2g",

                            from,
                            to,
                            date.fmt="%Y-%m-%d",

                            main="Calendar Heat Map",
                            #col= RColorBrewer::brewer.pal(n=8, name="RdYlBu"), # from Mauricio Zambrano-Bigiarini. A color palette, i.e. a vector of n contiguous colors generated by functions like rainbow, heat.colors, topo.colors, bpy.colors or one or your own making, perhaps using colorRampPalette. If none is provided, brewer.pal(n=8, name="RdYlBu") is used
                            col=colorRampPalette(c("red", "orange", "yellow", "white", "lightblue2", "deepskyblue", "blue3"), space = "Lab")(8), # from Mauricio Zambrano-Bigiarini. A color palette, i.e. a vector of n contiguous colors generated by functions like rainbow, heat.colors, topo.colors, bpy.colors or one or your own making, perhaps using colorRampPalette. If none is provided, brewer.pal(n=8, name="RdYlBu") is used
                            cuts,                                              # from Mauricio Zambrano-Bigiarini. Numeric, indicating the values used to divide the range of 'x' in the legend of colours. If not provided, it is automatically selected as a function of 'lenght(col)'
                            cuts.dec=0,                                        # from Mauricio Zambrano-Bigiarini. Number of decimal places used to present the numbers that divide the range of 'x' in the legend of colours
                            cuts.labels,                                       # from Mauricio Zambrano-Bigiarini. Character indicating the label to be used in the ccolour legend for each one of the values defined by 'cuts'. If not provided, as.character(cuts)' is used
                            cuts.style=c("fisher", "equal", "pretty", "fixed", "sd", "quantile", "kmeans", "bclust", "mzb"), # discarded becsue takes too much time or not alway provide the required number of classes: "dpih", "headtails", "hclust",  "jenks"
                            legend.title="",                                   # from Mauricio Zambrano-Bigiarini. text to be displayed above the legned of colours (e.g., showing the measurement units of the raster being displayed)
                            legend.fontsize=15,                                # from Mauricio Zambrano-Bigiarini. The size of text (in points) used in the legend 
                    
                            #varname="Values",
                            #date.fmt = "%Y-%m-%d", 

                            do.png=FALSE,
                              png.fname="mypng.png",
                              png.width=1500,
                              png.height=900,
                              png.pointsize=12,
                              png.res=90,
                            do.pdf=FALSE,
                              pdf.fname="mypdf.pdf",
                              pdf.width=11,
                              pdf.height=8.5,
                              pdf.pointsize=12,
                            ...) {

    # checking 'cuts.style'
    cuts.style <- match.arg(cuts.style)

    #require(lattice)
    #require(grid)
    #require(chron)
    dates  <- stats::time(x)
    values <- zoo::coredata(x)

    #if (class(dates) == "character" | class(dates) == "factor" ) {
    #   dates <- strptime(dates, date.fmt)
    #} # IF end

    # Checking daily time frequency
    if (sfreq(x) != "daily")
      stop("Invalid argument: 'sfreq(x)' must be 'daily', but it is '", sfreq(x), "' !")

    ####################################################################################
    # Lines 95-132 are taken and adpated from izoo2rzoo.R to check 'from' and 'to'
    ####################################################################################
  
    # If the index of 'p' is character, it is converted into a Date object
    if ( class(time(x))[1] %in% c("factor", "character") ) {
      dt <- try(as.Date(time(x), format=date.fmt))
      if("try-error" %in% class(dt) || is.na(dt)) {
        stop("Invalid argument: format of 'time(x)' is not compatible with 'date.fmt' !")
      } else time(x) <- as.Date(time(x), format=date.fmt)
    } # IF end      

    # If 'from' was given as Date, but 'x' is sub-daily
    if (!missing(from)) {
      dt <- try(as.Date(from, format=date.fmt))
      if("try-error" %in% class(dt) || is.na(dt)) {
        stop("Invalid argument: format of 'from' is not compatible with 'date.fmt' !")
      } else from <- as.Date(from, format=date.fmt)

      if (from > to) stop("Invalid argument: 'from > to' !")

      if (from > end(x)) stop("Invalid argument: 'from > end(x)' !")

      # Selecting only those data starting in 'from'
      x <- window(x, start=from)      
    } # IF end

    # If 'to' was given as Date, but 'x' is sub-daily
    if (!missing(to)) {
      dt <- try(as.Date(to, format=date.fmt))
      if("try-error" %in% class(dt) || is.na(dt)) {
        stop("Invalid argument: format of 'to' is not compatible with 'date.fmt' !")
      } else to <- as.Date(to, format=date.fmt)

      if (to < from ) stop("Invalid argument: 'to < from' !")

      if (to < start(x) ) stop("Invalid argument: 'to < start(x)' !")

      # Selecting only those data ending in 'to'
      x <- window(x, end=to)      
    } # IF end    
    ####################################################################################


    # Checking that the maximum amount of years is 6 
    if (length(x) > 2191) # 365 day/year * 6 years + 1 day in a leap year
      stop("Invalid argument: length(x)' must be less than 2191, but it is '", length(x), "' !")

    ######################################################################################
    # Lines 69-97 are from Mauricio Zambrano-Bigiarini
    ######################################################################################
    if (missing(cuts)) {
      ncuts <- length(col) + 1
      temp  <- as.numeric(values)
      temp <-  temp[!is.nan(temp)]

      if (cuts.style=="mzb") {
        probs <- seq( 0, 1, length.out=ncuts )      
        lcuts <- round( stats::quantile( temp, probs=probs, na.rm=TRUE), cuts.dec)
        if ( length(unique(lcuts)) != ncuts )
          lcuts <- .findcuts4plot(x=temp, n=ncuts)
      } else {
          lcuts <- classInt::classIntervals(temp, n=length(col), dataPrecision=cuts.dec, style = cuts.style )[["brks"]]
          lcuts <- round(lcuts, cuts.dec)
          ncuts <- length(lcuts)
          col   <- col[1:(ncuts-1)]
        }
    } else {
        if ( (length(col)+1) != length(cuts) ) {
          stop("Invalid argument: 'length(col)+1 != length(cuts)' ") 
        } else lcuts <- cuts
      } # ELSE end


    if ( missing(cuts.labels) ) {
      cuts.labels <- as.character(lcuts)
    } else {
        if (length(lcuts) != length(cuts.labels) )
          stop("Invalid argument: 'length(cuts) != length(cuts.labels)' ") 
      } # ELSE end
    ######################################################################################


    caldat   <- data.frame(value = values, dates = dates)
    min.date <- as.Date(paste(format(min(dates), "%Y"),
                        "-1-1",sep = ""))
    max.date <- as.Date(paste(format(max(dates), "%Y"),
                         "-12-31", sep = ""))
    dates.f  <- data.frame(date.seq = seq(min.date, max.date, by="days"))

    # Merge moves data by one day, avoid
    caldat <- data.frame(date.seq = seq(min.date, max.date, by="days"), value = NA)
    dates  <- as.Date(dates)
    caldat$value[match(dates, caldat$date.seq)] <- values

    caldat$dotw  <- as.numeric(format(caldat$date.seq, "%w"))
    caldat$woty  <- as.numeric(format(caldat$date.seq, "%U")) + 1
    caldat$yr    <- as.factor(format(caldat$date.seq, "%Y"))
    caldat$month <- as.numeric(format(caldat$date.seq, "%m"))
    yrs          <- as.character(unique(caldat$yr))
    d.loc        <- as.numeric()
    for (m in min(yrs):max(yrs)) {
      d.subset <- which(caldat$yr == m)
      sub.seq <- seq(1,length(d.subset))
      d.loc <- c(d.loc, sub.seq)
    }
    caldat <- cbind(caldat, seq=d.loc)

    #color styles
    r2b <- c("#0571B0", "#92C5DE", "#F7F7F7", "#F4A582", "#CA0020") #red to blue
    r2g <- c("#D61818", "#FFAE63", "#FFFFBD", "#B5E384") #red to green
    w2b <- c("#045A8D", "#2B8CBE", "#74A9CF", "#BDC9E1", "#F1EEF6") #white to blue
                
    #assign("col.sty", get(color))
    #calendar.pal <- colorRampPalette((col.sty), space = "Lab")
    def.theme    <- lattice::lattice.getOption("default.theme")

    cal.theme <- function() {
      theme <- list(
                    strip.background = list(col = "transparent"),
                    strip.border = list(col = "transparent"),
                    axis.line = list(col="transparent"),
                    par.strip.text=list(cex=0.8))
    } # 'cal.theme' END

    lattice::lattice.options(default.theme = cal.theme)

    yrs <- (unique(caldat$yr))
    nyr <- length(yrs)
    
    print(
      cal.plot <- lattice::levelplot(value~woty*dotw | yr, data=caldat,
          as.table=TRUE,
          aspect=.12,
          layout = c(1, nyr%%7),
          between = list(x=0, y=c(1,1)),
          strip=TRUE,
          main = main,
          scales = list(
          x = list(
                   at= c(seq(2.9, 52, by=4.42)),
                   labels = month.abb,
                   alternating = c(1, rep(0, (nyr-1))),
                   tck=0,
                   cex = 0.7),
         y=list(
                at = c(0, 1, 2, 3, 4, 5, 6),
                labels = c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday",
                           "Friday", "Saturday"),
                alternating = 1,
                cex = 0.6,
                tck=0)),
         xlim =c(0.4, 54.6),
         ylim=c(6.6,-0.6),

         #cuts= ncolors - 1,

         #col.regions = (calendar.pal(ncolors)),
         at=lcuts, col.regions=col,
         colorkey = list(title = paste0(legend.title, " \n"),
                         title.gpar=list(fontsize=legend.fontsize, font=1), #font=2 -> bold font
                         space = "right",
                         at=as.numeric(factor(lcuts)),                                     # equally spaced color bins in the legend
                         labels=list(labels=cuts.labels, at=as.numeric(factor(lcuts)), cex=0.9), # equally spaced color bins in the legend
                         col=col,
                         width = 0.6, height = 0.5                                  
                         ),

         xlab="" ,
         ylab=""
         #colorkey= list(col = calendar.pal(ncolors), width = 0.6, height = 0.5),
         #subscripts=TRUE
         ) #'levelplot' END
    ) # 'print' END

    panel.locs <- lattice::trellis.currentLayout()

    for (row in 1:nrow(panel.locs)) {
      for (column in 1:ncol(panel.locs)) {


        if (panel.locs[row, column] > 0) {

            lattice::trellis.focus("panel", row = row, column = column, highlight = FALSE)
            xyetc       <- lattice::trellis.panelArgs()
            subs        <- caldat[xyetc$subscripts,]
            dates.fsubs <- caldat[caldat$yr == unique(subs$yr),]
            y.start     <- dates.fsubs$dotw[1]
            y.end       <- dates.fsubs$dotw[nrow(dates.fsubs)]
            dates.len   <- nrow(dates.fsubs)
            adj.start   <- dates.fsubs$woty[1]

            for (k in 0:6) {
              if (k < y.start) {
                x.start <- adj.start + 0.5
              } else {
                 x.start <- adj.start - 0.5
                }
              if (k > y.end) {
                 x.finis <- dates.fsubs$woty[nrow(dates.fsubs)] - 0.5
              } else {
                  x.finis <- dates.fsubs$woty[nrow(dates.fsubs)] + 0.5
                }
              grid::grid.lines(x = c(x.start, x.finis), y = c(k -0.5, k - 0.5),
                               default.units = "native", 
                               gp=grid::gpar(col = "grey", lwd = 1))
            } # FOR end

            if (adj.start < 2) {
              grid::grid.lines(x = c( 0.5, 0.5), y = c(6.5, y.start-0.5),
                   default.units = "native", gp=grid::gpar(col = "grey", lwd = 1))
              grid::grid.lines(x = c(1.5, 1.5), y = c(6.5, -0.5), default.units = "native",
                   gp=grid::gpar(col = "grey", lwd = 1))
              grid::grid.lines(x = c(x.finis, x.finis),
                   y = c(dates.fsubs$dotw[dates.len] -0.5, -0.5), default.units = "native",
                   gp=grid::gpar(col = "grey", lwd = 1))
              if (dates.fsubs$dotw[dates.len] != 6) {
                  grid::grid.lines(x = c(x.finis + 1, x.finis + 1),
                  y = c(dates.fsubs$dotw[dates.len] -0.5, -0.5), default.units = "native",
                  gp=grid::gpar(col = "grey", lwd = 1))
              }
              grid::grid.lines(x = c(x.finis, x.finis),
                    y = c(dates.fsubs$dotw[dates.len] -0.5, -0.5), default.units = "native",
                    gp=grid::gpar(col = "grey", lwd = 1))
            } # IF end

            for (n in 1:51) {
              grid::grid.lines(x = c(n + 1.5, n + 1.5),
                               y = c(-0.5, 6.5), default.units = "native", 
                               gp=grid::gpar(col = "grey", lwd = 1))
            } # FOR end
            x.start <- adj.start - 0.5

            if (y.start > 0) {
              grid::grid.lines(x = c(x.start, x.start + 1),
                y = c(y.start - 0.5, y.start - 0.5), default.units = "native",
                gp=grid::gpar(col = "black", lwd = 1.75))
              grid::grid.lines(x = c(x.start + 1, x.start + 1),
                y = c(y.start - 0.5 , -0.5), default.units = "native",
                gp=grid::gpar(col = "black", lwd = 1.75))
              grid::grid.lines(x = c(x.start, x.start),
                y = c(y.start - 0.5, 6.5), default.units = "native",
                gp=grid::gpar(col = "black", lwd = 1.75))
             if (y.end < 6 ) {
              grid::grid.lines(x = c(x.start + 1, x.finis + 1),
               y = c(-0.5, -0.5), default.units = "native",
               gp=grid::gpar(col = "black", lwd = 1.75))
              grid::grid.lines(x = c(x.start, x.finis),
               y = c(6.5, 6.5), default.units = "native",
               gp=grid::gpar(col = "black", lwd = 1.75))
               } else {
                  grid::grid.lines(x = c(x.start + 1, x.finis),
                   y = c(-0.5, -0.5), default.units = "native",
                   gp=grid::gpar(col = "black", lwd = 1.75))
                  grid::grid.lines(x = c(x.start, x.finis),
                   y = c(6.5, 6.5), default.units = "native",
                   gp=grid::gpar(col = "black", lwd = 1.75))
                 }
             } else {
                  grid::grid.lines(x = c(x.start, x.start),
                                   y = c( - 0.5, 6.5), default.units = "native",
                                   gp=grid::gpar(col = "black", lwd = 1.75))
               } # ELSE end

             if (y.start == 0 ) {
              if (y.end < 6 ) {
              grid::grid.lines(x = c(x.start, x.finis + 1),
               y = c(-0.5, -0.5), default.units = "native",
               gp=grid::gpar(col = "black", lwd = 1.75))
              grid::grid.lines(x = c(x.start, x.finis),
               y = c(6.5, 6.5), default.units = "native",
               gp=grid::gpar(col = "black", lwd = 1.75))
               } else {
                  grid::grid.lines(x = c(x.start + 1, x.finis),
                   y = c(-0.5, -0.5), default.units = "native",
                   gp=grid::gpar(col = "black", lwd = 1.75))
                  grid::grid.lines(x = c(x.start, x.finis),
                   y = c(6.5, 6.5), default.units = "native",
                   gp=grid::gpar(col = "black", lwd = 1.75))
                 } # ELSE end
             } # IF end

            for (j in 1:12) {
               last.month <- max(dates.fsubs$seq[dates.fsubs$month == j])
               x.last.m <- dates.fsubs$woty[last.month] + 0.5
               y.last.m <- dates.fsubs$dotw[last.month] + 0.5
               grid::grid.lines(x = c(x.last.m, x.last.m), y = c(-0.5, y.last.m),
                 default.units = "native", gp=grid::gpar(col = "black", lwd = 1.75))
               if ((y.last.m) < 6) {
                  grid::grid.lines(x = c(x.last.m, x.last.m - 1), y = c(y.last.m, y.last.m),
                   default.units = "native", gp=grid::gpar(col = "black", lwd = 1.75))
                 grid::grid.lines(x = c(x.last.m - 1, x.last.m - 1), y = c(y.last.m, 6.5),
                   default.units = "native", gp=grid::gpar(col = "black", lwd = 1.75))
               } else {
                  grid::grid.lines(x = c(x.last.m, x.last.m), y = c(- 0.5, 6.5),
                   default.units = "native", gp=grid::gpar(col = "black", lwd = 1.75))
                 } # ELSE end
             } # FOR end
        } # IF end

      } # FOR end

      lattice::trellis.unfocus()

    } # FOR end

  lattice::lattice.options(default.theme = def.theme)

} # 'calendarHeatmap' END


################################################################################
###                             .findcuts4plot                                ###
################################################################################
### Author     : Mauricio Zambrano-Bigiarini                                 ###
################################################################################
### Started    :  20-Dic-2017                                                ###
### Updates    :  12-May-2018                                                ###
###               01-Aug-2019                                                ###
###               11-Jan-2020                                                ###
###               09-May-2021 ; 19-Jun-2021                                  ###
###               10-May-2022                                                ###
################################################################################

# Function used to compute suitable cutting points to be used as separators of
# different colour classes used with rasterVIS::levelplot.
# This function is usually necessary when the distribution of x is highly skewed 
# (too many zeros) to find the 'cuts' values using the normal quantiles
# dec: decimals used to round each one of the output quantiles
.findcuts4plot <- function(x, n=8, thr=0, dec=0) {

  # Finding the minimum value in 'x'
  lmin <- round( min(x, na.rm = TRUE), dec)
  
  # Finding the quantile corresponding to the first value higher than 'thr'
  probs <- seq(0.01, 1, by=0.01)
  q     <- stats::quantile(x, probs=probs, na.rm=TRUE)
  index <- which(q > thr)
  qmin  <- index[1]/100
  
  # Finding n-1 quantiles that are higher than 'thr', only for 'x > thr'
  probs <- seq(qmin, 1, length.out = n-1)
  out   <- unique( round( stats::quantile(x[ x > thr ], probs=probs, na.rm=TRUE), dec) )

  # Adding 'lmin'
  out <- c( lmin, out )

  out <- unique(out)

  # cheking length of the output
  if (length(out) < n)
    out <- seq(from=lmin, to=max(out, na.rm=TRUE), length.out=n)
  
  # Creating the output object
  return( out )
  
} # '.findcuts4plot' END


# ## Example of use: Plot financial data

# plot <- FALSE
# ## This code is not run.
# if (plot) {

# #create faux data; skip this to use data from a file or stock data
# #ndays <- 1500 #set number of days
# #dates <- as.POSIXlt(seq(Sys.Date()- ndays, Sys.Date() - 1, by="days"))
# #vals <- runif(ndays, -100, 100)

# #stock data:
# stock      <- "MSFT"
# start.date <- "2010-01-12"
# end.date   <- Sys.Date()
# quote      <- paste("http://ichart.finance.yahoo.com/table.csv?s=",
#                 stock,
#                 "&a=", substr(start.date,6,7),
#                 "&b=", substr(start.date, 9, 10),
#                 "&c=", substr(start.date, 1,4),
#                 "&d=", substr(end.date,6,7),
#                 "&e=", substr(end.date, 9, 10),
#                 "&f=", substr(end.date, 1,4),
#                 "&g=d&ignore=.csv", sep="")
# stock.data <- read.csv(quote, as.is=TRUE)

# # Plot as calendar heatmap
# calendarHeat(stock.data$Date, stock.data$Adj.Close, varname="MSFT Adjusted Close")


# }
