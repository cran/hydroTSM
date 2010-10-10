#####################################################################
#                'lib_TSM_in_HydrologicalModelling'                 #
#####################################################################
#  Library for Time Series Management in Hydrological Modelling     #
#  Author      : Mauricio Zambrano Bigiarini                        #
#  Stared      : July 2008                                          #
#  Version     : 0.1.0 : 07-Sep-2009                                #
#  Version     : 0.1.1 : 15-Sep-2009                                #
#  Version     : 0.1.2 : 01-Oct-2009                                #
#  Version     : 0.1.3 : 13-Oct-2009                                #
#  Version     : 0.1.4 : 11-Nov-2009                                #
#  Version     : 0.1.5 : 13-Nov-2009                                #
#  Version     : 0.1.6 : 30-Nov-2009                                #
#  Version     : 0.2.0 : 07-Oct-2010                                #
#  Version     : 0.2.1 : ongoing...                                 #
#  Last Update : 06-Oct-2010                                        #
#####################################################################

# The function 'matrixplot' requires the 'lattice' package


########################################################################
# Included functions:
# 01) ma.default          : Moving (sliding) average function
# 02) ma.zoo              : Moving (sliding) average function for zoo objects

# 03) dwi.che             : Days with Information in a CHE file
# 04) dwi.default         : Days with Information in a zoo object
# 29) dwi.data.frame      : Days with information in each station stored in a data frame
# 14) .CHE_daily_analysis : Given a directory with CHE files, it computes:       
#                           1) Number of Days with Information in each CHE file#
#                           2) Matrix with the Time series for each station, OR# 
#                              within the time period comprised between 'from' and 'to'                
# 15) plot.daily.analysis : Plot Analysis of the  number of gauging stations vs 
#                           the percent of days with information with respect to the total amount of days in the period
# 18) .DailyTS2DB         : Database (.csv) generation from all the CHE files (time series) in a directory
#                           The fields in the resulting db file are: ID, ID_Station, Date, Value 
# 35) matrixplot          : Plots a color matrix representing the amount of days with information in a set of gauging stations
# 22) .Qcsv2zoo           : qXXX.csv -> ZOO
# 05) .che2zoo            : CHE -> ZOO
# 06) izoo2rzoo           : Irregular Zoo -> Regular Zoo	
# 20) vector2zoo          : Transforms a numerical vector 'x' with a corresponding 'Date' vector (any format) into a 'zoo' object
# 33) sfreq               : String with the sampling frequency of a ts/zoo object ("Daily", "Monthly", "Annual")
# 28) time2season         : Transforms a vector of dates into a vector of seasons (summer, winter, autumm, spring)
# 25) sname2ts            : Given a data.frame whose columns contains the time series of several gauging stations, 
#                           it takes the name of one gauging station and extracts a time series with daily, monthly or annual time step  

# 07) daily2monthly       : daily time series into monthly ones
# 08) daily2annual        : daily time series into annual ones
# 09) monthly2annual      : monthly time series into annual ones
#     dm2seasonal         : Generic function for computing seasonal values for every year of a daily/monthly zoo object  

# 10) .hydroplotts        : Given 3 zoo ts (daily, monthly and annual), it plots 3 ts (lines) with the daily, monthly, and annual values
# 11) .hydroplotboxplot   : Given 3 zoo ts (daily, monthly and annual), it plots 3 boxplots with the daily, monthly, and annual ts
# 12) .hydroplothist      : Given 3 zoo ts (daily, monthly and annual), it plots 3 histograms with the daily, monthly, and annual ts
# 13) hydroplot           : Given 1 'zoo' daily ts, it plots 3 ts, 3 boxplots and 3 histograms with the daily, monthly, and annual ts
# 21) sname2plot          : Given a data.frame whose columns contains the time series of several gauging stations, it  
#                           takes the name of one gauging station and Plots 9 graphs (see 'hydroplot' description) 
    
# 16) diy                  : Given a numeric value of a year, it generates:       
#                            1) a vector with all the days (dates) within the year, OR 
#		                     2) the amount of days in the year  
# 17) dip                  : Given any starting and ending dates, it generates:  
#                            1) a vector with all the days between the two dates, OR #
#		                     2) the amount of days between the two dates     
# 19) infillxy             : Fills in ALL the 'NA' values in 'obs' with the corrensponding values in 'sim'
# 23) interpol1            : Spatial interpolation using a modified IDW, where the Pearson's product-moment coefficient of 
#                            correlation between the time series of all teh stations is used instead of the spatial distance
# 24) rm1stchar            : deletes the first characther of each element of 'x' 
# 26) stdx                 : standarizes a vector or matrix, i.e., scales all the values, in in a way that the 
#                            transformed values will be within the range [0,1]. z = scale(x) = [ x - xmin ] / [ xmax - xmin ]
# 27) istdx                : Transforms back a standarized vector/matrix into their original values, 
#                            i.e., re-scales all the values in the [0,1] interval to the original range of values. x = re-scale(z) = z*[ zmax - zmin ] + xmin
# 30) extractzoo           : Extracts from a zoo object, all the values belonging to a
#                            given month, year or weather season   
# 37) fdc                  : Flow Duration Curve, computation and plot      
# 38) dmc                  : Monthly Double-Mass Curve for daily precipitation or streamflow data (Homogeneity test)
# 39) drawxaxis            : It draws an X axies with daily, monthly, or annual time marks
# 41) hydropairs           : Visualization of a Correlation Matrix 
# 42) smry                 : 13 summary statistics of numeric objects
# 43) hypsometric          : Plot of the hypsometric curve of a given DEM
# 44) monthlyfunction      : Generic function for applying any R function to ALL the values in 'x' belonging to a given month 
# 45) seasonalfunction     : Generic function for applying any R function to summarize the seasonal values of 'x' 
# 46) zoo2RHtest           : creates the input file to the 'RHtest_dlyPrcp.r' script, for homeneity testing of climatological time series
# 47) dwdays               : Generic function for computing the average amount of wet days by month in a daily time series of precipitation
# 48) fdcu                 : Flow Duration Curve with Uncertainty Bounds
# 49) mip                  : Given any starting and ending dates, it generates:  
#                            1) a vector with all the months between the two dates, OR #
#		                     2) the amount of months between the two dates  
# 50) yip                  : Given any starting and ending dates, it generates:  
#                            1) a vector with all the years between the two dates, OR #
#		                     2) the amount of years between the two dates   


########################### START of ex-lib_Plot.R #####################


################################################################################
# 'drawxaxis': It draws an X axies with daily, monthly, or annual time marks   #
################################################################################
# Started on 2008          #
# Modified: March 2009     #
############################ 
# 'x'         : time series that will be plotted using the X axis that will be draw
#               class(x) must be 'ts' or 'zoo'
# 'tick.tstep': Character indicating the time step that have to be used for 
#               putting the ticks ont he time axis. 
#               Possible values are: 'days', 'months', 'years' 
# 'lab.tstep' : Character indicating the time step that have to be used for 
#               putting the labels ont he time axis. 
#               Possible values are: 'days', 'months', 'years' 
# 'cex.axis'  : magnification of axis annotation relative to cex 
drawxaxis <- function(x, tick.tstep="months", lab.tstep="years", cex.axis=1, ... ) {

 # Checking that the user provied a valid argument for 'tick.tstep'       
 if (is.na(match(tick.tstep, c("days", "months", "years") ) ) ) 
         stop("Invalid argument: 'tick.tstep' must be in c('days', 'months', 'years')")   
         
 # Checking that the user provied a valid argument for 'lab.tstep'       
 if (is.na(match(lab.tstep, c("days", "months", "years") ) ) ) 
         stop("Invalid argument: 'lab.tstep' must be in c('days', 'months', 'years')") 
         
 # Computing the position of the tick on the Time axis
 if ( (tick.tstep == "days") | (tick.tstep == "months") ) {
 
   start <- range( time(x) )[1]
   end   <- range( time(x) )[2] 
 
   # Computes daily or monthly ticks on the X axis, without labels
   tt <- seq(from=start, to=end, by = tick.tstep)    
  
 } else if (tick.tstep == "years") {
 
      if ( class(time(x)) == "Date" ) {
          start <- as.Date( range( time(x) )[1])
          end   <- as.Date( range( time(x) )[2])
          # Computes the ticks for Annual series in the x axis, without labels
          tt <- seq(from=start, to=end, by = tick.tstep) 
      } else if ( class(time(x)) == "character" ) {  
           start <- as.numeric( range( time(x) )[1])
           end   <- as.numeric( range( time(x) )[2]) 
           # Computes the ticks for Annual series in the x axis, without labels       
           tt   <- start:end 
        } # ELSE END    
          
   } # ELSE END   
 
 # Draws the ticks for the time series in the x axis, without labels
 Axis(side = 1, at = tt, labels = FALSE, ...)    
 
 # Computes the string of the labels of the X axis
 if (lab.tstep == "days" | lab.tstep == "months" ) { 
 
   # Computes the position of the labels of the X axis
   tt <- seq(from=start, to=end, by = lab.tstep) 
     
   if (lab.tstep == "days") { 
     labs <- format(tt, "%Y-%m-%d")   
   } else labs <- format(tt, "%b")   
   
 } else if (lab.tstep == "years") {
 
   if ( class(time(x)) == "Date" ) {
          # Computes the position of the labels of the X axis
          tt <- seq(from=start, to=end, by = "months")
          
          tt <-  subset(tt, format(tt, "%m")=="01")
          
          labs <- format(tt, "%Y") 
      } else if ( class(time(x)) == "character" ) {          
           start <- as.numeric( range( time(x) )[1])
           end   <- as.numeric( range( time(x) )[2])  
           # Computes the position of the labels of the X axis                  
           tt   <- start:end                         
           labs <- as.character(tt)
        } # ELSE END
    
 }  # IF END
  
 # Draws the labels corresponding to the selected ticks in the X axis
 Axis(side = 1, at = tt, labels = labs, tcl = -0.7, cex.axis = cex.axis, ...)
 
} # 'drawxaxis' END


###########################END of ex-lib_Plot.R#########################


#####################################################
# Generic Moving (sliding) Average function         #
#####################################################
ma <-function(x, ...) UseMethod("ma")


#####################################################
# Default Moving (sliding) Average function         #
#####################################################
# A vector with the moving average computed using a window of 'win.len' elements

# x      : time series with n elements
# win.len: number of terms that will be considered in the mean. It have to be odd

# result : a vector with the moving average termns. The length of the resulting vector
#          is the same of 'x', but the first and last (win.len-1)/2 elements will
#          be NA.

ma.default <- function (x, win.len, FUN=mean,...) {

  if (ceiling(win.len)/2 != win.len/2)
      stop("Invalid argument: 'win.len' must be of odd")

  return( filter(x, rep(1/win.len, win.len), method="convolution", sides=2) ) 

} # 'ma.default' end


####################################################################
# 	ma.zoo:	Moving Average of a DAILY regular time series,         #  
#           by default using a window width =365 (Annual Average)  #
####################################################################

# 'x'       : zoo variable
# 'win.len' : window width. It have to be odd
# 'FUN'     : Function that have to be applied for computing the moving average. Usually, FUN MUST be "mean"

# 'result'  : a vector with the moving average termns. The length of the resulting vector
#             is less than the length of 'x', with (win.len-1)/2 missing elements at the begining and end of 'x'
ma.zoo <- function(x, win.len, FUN=mean,... ) {  

  if (ceiling(win.len)/2 != win.len/2)
      stop("Invalid argument: 'win.len' must be of odd")

  # Requiring the zoo library
  require(zoo)  

  # Generating an Moving Average time series, with a window width win.len1
  return ( rollapply(x, width=win.len, FUN, by.colum=FALSE) )
 
} # 'ma.zoo' end
 
 
#####################################################
# Generic Days With Information function            #
#####################################################
dwi <-function(x, ...) UseMethod("dwi")
               
 
###################################################
#            Zoo Days with Information            #
###################################################
# This function generates a table indicating the number of days 
# with information (<>NA's) within a zoo object,
# aggregated by: Year, Month or Month by Year

# x        : variable of type 'zoo'
# out.unit : aggregation time for the computation of the amount of days with info.
#	         Valid values are:
#            -) "month": monthly; 
#            -) "year" : annual; 
#            -) "mpy"  : month per year
# from     : Character indicating the starting date for the values stored in all the files that 
#            will be read. It HAs to be in the format indicated by 'date.fmt'
# to       : Character indicating the starting date for the values stored in all the files that 
#            will be read. It HAs to be in the format indicated by 'date.fmt'
# date.fmt : Character indicating the date format in which you provide 'from' and 'to', e.g. "%d-%m-%Y"

dwi.default <- function(x, out.unit="years", from=range(time(x))[1], 
                        to=range(time(x))[2], date.fmt="%Y-%m-%d", tstep="days", ...) {
                        
    # Checking that 'class(x)==zoo'
    if (is.na(match(class(x), c("zoo") ) ) ) 
     stop("Invalid argument: 'x' must be of class 'zoo'")
                               
    # Checking the validity of the 'out.unit' argument
    if ( is.na( match(out.unit, c("years", "months", "mpy") ) ) ) {
         stop("Invalid argument value: 'out.unit' must be in c('years', 'months', 'mpy')" ) } 
 
    # Sequence of dates within the time period between 'from' and 'to'
    DateSeq <- seq( from=as.Date(from, format=date.fmt), 
                    to=as.Date(to, format=date.fmt), by=tstep )
                        
    # Selecting only those data that are within the time period between 'from' and 'to'
    x.sel <- x[ as.Date( time(x), format=date.fmt) %in% DateSeq]	
    # Also is possible to use: x.sel <- window(x, start=as.Date(from, format=date.fmt), end=as.Date(to, format=date.fmt) )  
    
    # Computing the Starting and Ending Year of the analysis
    Starting.Year <- as.numeric(format(as.Date(from, format=date.fmt), "%Y"))
    Ending.Year   <- as.numeric(format(as.Date(to, format=date.fmt), "%Y"))
  
    # Amount of Years belonging to the desired period
    nyears <- Ending.Year - Starting.Year + 1
     
    if (out.unit == "months")   {
    
         a <- numeric(12)
         
         a[1:12] <- sapply(1:12, function(j,y) {         
                              tmp         <- extractzoo(y, trgt= j)                              
                              nona.index  <- which(!is.na(tmp))                              
                              a[j] <- length( nona.index )                              
                              }, y = x.sel) 
         names(a) <- month.abb                   
         return(a)
                
     } # IF end
     
    else if (out.unit == "years") {
         
         a <- numeric(nyears)
          
         a[1:nyears] <- sapply(Starting.Year:Ending.Year, function(j,y) {                                   
                               tmp         <- extractzoo(y, trgt= j)                              
                               nona.index  <- which(!is.na(tmp))                              
                               a[j] <- length( nona.index )                              
                               }, y = x.sel)  
          
         names(a) <- as.character(Starting.Year:Ending.Year)
         
         return(a) 
         
     } # ELSE IF end
     
     else if (out.unit == "mpy") {  
         
         a <- matrix(data=NA,nrow=nyears, ncol=12)
          
         a <- sapply(Starting.Year:Ending.Year, function(i,y) {    
                                                    
                               tmp         <- extractzoo(y, trgt= i)  
                               
                               a[i-Starting.Year+1,1:12] <- sapply(1:12, function(j,y) {         
                                             tmp2        <- extractzoo(y, trgt= j)                              
                                             nona.index  <- which(!is.na(tmp2))  
                                             
                                             a[i-Starting.Year+1,j] <- length( nona.index )                              
                                         }, y = tmp)                                                               
                        }, y = x.sel)
                        
         a <- t(a)
         
         #Change the names of the columns of the matrix
         colnames(a) <- month.abb  
         #Change the names of the rows of the matrix   
         rownames(a) <- as.character(Starting.Year:Ending.Year)    
              
         return(a)
     
    }  # ELSE IF end
 
 } # 'dwi.default' end
 

#########################################################################
# dwi.data.frame: days with info in each station stored in a data frame #
#########################################################################
#                             March 21th, 2009                          #
#########################################################################
# This function generates a table indicating the number of days 
# with information (<>NA's) within a data.frame

# 'x'         : variable of type 'data.frame'
# out.unit    : aggregation time for the computation of the amount of days with info.
#	            Valid values are:
#               -) "month": monthly; 
#               -) "year" : annual; 
# from        : starting date for detection of days with inormation
# to          : date format that will be used in the output variable
# date.fmt    : date format of "from" and "to". For CHE files, the format must be "%d-%m-%Y"
# 'dates'     : "numeric", "factor", "Date" indicating how to obtain the 
#               dates for correponding to the 'sname' station
#               If 'dates' is a number, it indicates the index of the column in 
#                 'x' that stores the dates
#               If 'dates' is a factor, it have to be converted into 'Date' class, 
#                 using the date format  specified by 'date.fmt'
#               If 'dates' is already of Date class, the following line verifies that
#                 the number of days in 'dates' be equal to the number of element in the 
#                 time series corresponding to the 'st.name' station
# 'verbose'  : logical; if TRUE, progress messages are printed 
                        
dwi.data.frame <- function(x, out.unit="years", from, to, 
                           date.fmt="%Y-%m-%d", tstep="days", dates=1, verbose=TRUE,...) {
     
  # Checking the validity of the 'out.unit' argument
  if ( is.na( match(out.unit, c("years", "months") ) ) ) {
       stop("Invalid argument value: 'out.unit' must be in c('years', 'months')" ) }
       
  # Checking that the user provied a valid argument for 'dates'       
  if (missing(dates)) {
      stop("Missing argument: 'dates' must be provided") 
  } else  
    {    
     # Checking that the user provied a valid argument for 'dates'  
     if (is.na(match(class(dates), c("numeric", "factor", "Date")))) 
         stop("Invalid argument: 'dates' must be of class 'numeric', 'factor', 'Date'")
         
     # Verification that the number of days in 'dates' be equal to 
     # the number of elements in 'x'
     if ( ( class(dates) == "Date") & (length(dates) != nrow(x) ) ) 
          stop("Invalid argument: 'length(dates)' must be equal to 'nrow(x)'") 
    } # ELSE end
        
  # If 'dates' is a number, it indicates the index of the column of 'x' that stores the dates
  # The column with dates is then substracted form 'x' for easening the further computations
  if ( class(dates) == "numeric" ) {    
    tmp   <- dates
    dates <- as.Date(x[, dates], format= date.fmt)
    x     <- x[-tmp]
  }  else 
      # If 'dates' is a factor, it have to be converted into 'Date' class, 
      # using the date format  specified by 'date.fmt'
      if ( class(dates) == "factor" ) {
	    dates <- as.Date(dates, format= date.fmt)
	  } # IF end
		  
  # Checking the validity of the 'from' argument
  if (missing(from)) { from <- dates[1]
  } else if ( is.na( match(class(from), c("Date", "character") ) ) ) {
            stop("Invalid argument value: 'class(from)' must be in c('Date', 'character')" ) }
		 
  # Checking the validity of the 'to' argument
  if (missing(to)) { to <- dates[length(dates)]
  } else if ( is.na( match(class(to), c("Date", "character") ) ) ) {
            stop("Invalid argument value: 'class(to)' must be in c('Date', 'character')" ) }     

  # Sequence of dates within the time period between 'from' and 'to'
  DateSeq <- seq( from=as.Date(from, format=date.fmt), 
                  to=as.Date(to, format=date.fmt), by=tstep )
                        
  # Selecting only those data that are within the time period between 'from' and 'to'
  x.sel <- x[dates %in% DateSeq, ]	
  # Also is possible to use: x.sel <- window(x, start=as.Date(from, format=date.fmt), end=as.Date(to, format=date.fmt) ) 
  
  # Computing the Starting and Ending Year of the analysis
  Starting.Year <- as.numeric(format(as.Date(from, format=date.fmt), "%Y"))
  Ending.Year   <- as.numeric(format(as.Date(to, format=date.fmt), "%Y"))
  
  # Amount of Years belonging to the desired period
  nyears <- Ending.Year - Starting.Year + 1
  
  # Amount of stations in 'x.sel'
  nstations <- ncol(x.sel)
  
  # NAme of the stations in 'x'
  snames <- colnames(x)
     
  if (out.unit == "years") {
    z <- matrix(data=NA, nrow= nstations, ncol=nyears) 
  } else if (out.unit == "months") {
    z <- matrix(data=NA, nrow=nstations, ncol=12) 
    }  # ELSE end
    
   
  z <- sapply(1:ncol(x.sel), function(j, y) {
              #y[j] <- length( subset(y[,j], !is.na(y[,j]) ) )
                  
              if (verbose) print( paste("Station: ", format(snames[j], width=6, justify="left"),
				                        " : ", format(j, width=3, justify="left"), "/", 
					                    nstations, " => ", 
					                    format(round(100*j/nstations,2), width=6, justify="left"), 
					                    "%", sep=""), quote=FALSE )
                                  
              tmp  <- vector2zoo(x=y[,j], dates=dates, date.fmt=date.fmt)
                  
              z[j,] <- dwi.default(x=tmp, out.unit=out.unit, from=from, to=to, date.fmt=date.fmt)
                 
              }, y = x.sel) 
  
  colnames(z) <- snames
  
  return(z)
  
} # 'dwi.data.frame' END



#####################################################
#         Irregular Zoo -> Regular Zoo              #
#####################################################
# This function takes a time series of (very likely) irregular (with 
# missing dates) daily time series and then transforms it into a variable 
# regulary spaced, filling the voids with some value (by default: NA) 

# x         : time series of type zoo (very likely with some missing days)
# date.fmt  : character indicating the format in which the dates are stored in 'dates', e.g. "%Y-%m-%d"
#             ONLY required when class(dates)=="factor" or "numeric"
# from      : starting date for the merged output
# to	    : ending date for the merged output
# tstep     : time step in which are stored the values of 'x'

izoo2rzoo <- function(x, from= range(time(x))[1],
                      to= range(time(x))[2], date.fmt="%Y-%m-%d", tstep ="days" ) { 
 
     if (is.na(match(class(x), c("zoo")))) 
            stop("Invalid argument: 'x' must be of class 'zoo'") 
            
     # Requiring the Zoo Library (Zoo's ordered observations)
     require(zoo) 
     
     # Generating a daily-regular time series of NA's, 
     # just for being merged with the real daily-irregular time series
     z <- seq( from=as.Date(from, format=date.fmt), 
               to=as.Date(to, format=date.fmt), by= tstep ) # the default name of the column with the dates is "X_data"
     
     z.zoo       <- as.zoo( rep(NA, length(z)) )
     time(z.zoo) <- z 	#class(z) = Date
     
     # Selecting only those data that are within the time period between 'from' and 'to'
     x.sel <- x[ time(x) %in% z, ]
     
     # Creating a daily-regular time series with the read Precipitatoin values and NA's in those days without information
     x.merged <- merge(x.sel,z.zoo)

     # Returning as result only the column containing the Daily-Regular Time Series of Precipitation with NA's in the empy days
     return( x.merged[,1] )
 
} # 'izoo2rzoo' end


########################################################################
#  'vector2zoo': Transforms a numerical vector 'x' with a corresponding#
#             'Date' vector (any format) into a 'zoo' object           # 
#                       17-Dic-2008, 01-Oct-2009, 06-Oct-2010          #
########################################################################
#  Transform a numericl vectorial  and  its corresponding dates into 
#  a 'zoo' variable, for being used by other procedures of this library                                

# 'x':      : numeric vector
# 'dates'   : vector with a complete series of dates, in the same order of 'ts'
# 'date.fmt': format in which the dates are stored in 'dates'

# example:
# > x <- read.csv2("Ebro-Daily_TS_by_station-PP-Thr70-349stations-HistoricalPeriod.csv")
# > d <- vector2zoo(x[,2], dates=as.Date(x[,1]) )
# > summary(d)
vector2zoo <- function(x, dates, date.fmt="%Y-%m-%d") { 
  
  # Requiring the Zoo Library
  require(zoo)
  
  if (is.na(match(class(dates), c("Date", "character", "factor")))) 
        stop("Invalid argument: 'class(dates)' must be in c('Date', 'character', 'factor')")   
  
   if (is.na(match(class(dates), c("Date")))) 
      dates <- as.Date(dates, format= date.fmt)
 
  # Transforming into a 'zoo' type the values in the time series
  b <- as.zoo(x)
 
  # Setting as the date of the 'zoo' object the dates provided by 'dates'
  time(b) <- as.Date(dates, format= date.fmt)
 
  return( b )      
 
}  # 'vector2zoo' END



#################################################
# sfreq: Sampling frequency of a ts/zoo object  # 
#################################################
#              May 13th, 2009                   #
#################################################
# This function generates a table indicating the number of days 
# with information (<>NA's) within a data.frame

# 'x'        : variable of type 'zoo' or 'ts', with AT LEAST 2 elements, AND 
#              with a Daily, Monthly or Annual sampling frequency.
# 'min.year' : integer used for a correct identification of the sampling fequency 
#              when 'x' is an annual zoo (e.g.: time(x) = "1961") => the minimum possible years starts 
#              in 'min.year' 

# Result     : Possible values are:
#              -) 'daily'   : indicating that the sampling freqeuncy in 'x' is daily
#              -) 'monthly' : indicating that the sampling freqeuncy in 'x' is monthly
#              -) 'annual'  : indicating that the sampling freqeuncy in 'x' is annual
sfreq <- function(x, min.year=1800) {

  # Checking that 'class(x)==Date'
  if (is.na(match(class(x), c("zoo", "ts") ) ) ) 
     stop("Invalid argument: 'x' must be in c('zoo', 'ts')" )
     
  if ( length(x) < 2) stop("Invalid argument: 'length(x)' must be larger than 2 for finding its sampling frequency" )
  
  t1 <- time(x[1])
  t2 <- time(x[2])
  
  if ( ( class(t1) == "character" ) & ( as.numeric(t1) > min.year ) ) {
    sfreq <- "annual"
  } else       
      if ( ( t2 - t1 ) == 1)  {
        sfreq <- "daily"
      } else if ( ( t2 - t1 ) %in% c(28,29,30,31) )  {
        sfreq <- "monthly"
        } else if (  ( t2 - t1 ) %in% c(365, 366) )  {
          sfreq <- "annual"
        } else 
          stop("Invalid Argument: the sampling frequency of 'x' is not in c('daily', 'monthly', 'annual') " )
  
  return(sfreq)
  
} # 'sfreq' END


#####################################
#          daily2monthly            #
#####################################
# This function transform a DAILY regular (without missing days, 
#  but possible with days without information) time series into a MONTHLY one

# 'x'   : daily values that will be converted into annual ones.
#         class(x) must be 'zoo'
# 'FUN' : Function that have to be applied for transforming from daily to monthly time step
#         For precipitation FUN MUST be "sum"
#         For temperature and flow time series, FUN MUST be "mean"
# 'na.rm': Logical. Should missing values be removed?
#          TRUE : the monthly and annual values  are computed considering only those values different from NA
#          FALSE: if there is AT LEAST one NA within a year, the monthly and annual values are NA 

daily2monthly <-function(x, ...) UseMethod("daily2monthly")
      
daily2monthly.default <- function(x, FUN, na.rm=TRUE, ... ) { 

  # Checking the user provide a valid value for 'x'
  if (is.na(match(class(x), c("zoo")))) 
        stop("Invalid argument: 'x' must be of class 'zoo'") 
		
  # Checking the user provide a valid value for 'FUN'
  if (missing(FUN)) {
     stop("Missing argument value: 'FUN' must contain a valid function for aggregating the daily values") }
		  
  # Requiring the Zoo Library (Z’s ordered observations)
  require(zoo)
 
  # Checking the user provide a valid value for the sampling frequency of 'x'
  if (sfreq(x) != "daily") {
      stop(paste("Invalid argument: 'x' is not a daily ts, it is a ", sfreq(x), " ts", sep="") ) }		
   
  # Generating a Monthly time series of Total Monthly Precipitation (Monthly sum of daily values)
  tmp <-aggregate( x, by=as.Date( as.yearmon( time(x) ) ), FUN, na.rm= na.rm )
  
  # Getting the position of all the years in which there were no values
  # mean(NA:NA, na.rm=TRUE) == NaN
  nan.index <- which(is.nan(tmp))
  
  # Changing all the NaN's by NA's
  if ( length(nan.index) > 0 ) { tmp[nan.index] <- NA }
 
  return(tmp)         
 
} # 'daily2monthly.default' end



# 'dates'   : "numeric", "factor", "Date" indicating how to obtain the 
#             dates for correponding to the 'sname' station
#             If 'dates' is a number, it indicates the index of the column in 
#                'x' that stores the dates
#             If 'dates' is a factor, it have to be converted into 'Date' class, 
#                using the date format  specified by 'date.fmt'
#             If 'dates' is already of Date class, the following line verifies that
#                the number of days in 'dates' be equal to the number of element in the 
#                time series corresponding to the 'st.name' station
# 'date.fmt': format in which the dates are stored in 'dates'.
#             ONLY required when class(dates)=="factor" or "numeric"
# 'out.type': string that define the desired type of output. Possible values are
#             -) "data.frame": a data.frame, with as many columns as stations 
#                              are included in 'x'
#             -) "db"        : a data.frame, with 4 colums will be produced.
#                              The first column stores the ID of the station, 
#                              The second column stores the Year 
#                              The third column stores the month 
#                              The fourth colum stores the numerical values corresponding to the year and month specified in the two previous fields.
# 'out.fmt' : character, for selecting if the result will be 'numeric' or 'zoo'. Valid values are: c('numeric', 'zoo') 
# 'verbose'      : logical; if TRUE, progress messages are printed 
daily2monthly.data.frame <- function(x, FUN, na.rm=TRUE,
                                     dates, date.fmt="%Y-%m-%d", 
				     out.type="data.frame", 
				     out.fmt="numeric",
				     verbose=TRUE,...) {
    
  # Checking that the user provide a valid value for 'FUN'
  if (missing(FUN)) {
		 stop("Missing argument value: 'FUN' must contain a valid function for aggregating the daily values") }
         
    # Checking that the user provied a valid argument for 'out.type'  
  if (is.na(match( out.type, c("data.frame", "db") ) ) ) 
      stop("Invalid argument: 'out.type' must be in c('data.frame', 'db')")
		 
  # Checking that the user provied a valid argument for 'out.fmt'  
  if (is.na(match( out.fmt, c("numeric", "zoo") ) ) ) 
      stop("Invalid argument: 'out.fmt' must be in c('numeric', 'zoo')")
	  
  # Checking that the user provied a valid argument for 'dates'       
  if (missing(dates)) {
      stop("Missing argument: 'dates' must be provided") 
  } else      
     # Checking that the user provied a valid argument for 'dates'  
     if (is.na(match(class(dates), c("numeric", "factor", "Date")))) 
         stop("Invalid argument: 'dates' must be of class 'numeric', 'factor', 'Date'")
        
  # If 'dates' is a number, it indicates the index of the column of 'x' that stores the dates
  # The column with dates is then substracted form 'x' for easening the further computations
  if ( class(dates) == "numeric" ) {    
    tmp   <- dates
    dates <- as.Date(x[, dates], format= date.fmt)
    x     <- x[-tmp]
  }  # IF end 
  
  # If 'dates' is a factor, it have to be converted into 'Date' class, 
  # using the date format  specified by 'date.fmt'
  if ( class(dates) == "factor" ) dates <- as.Date(dates, format= date.fmt)
  
  # If 'dates' is already of Date class, the following line verifies that
  # the number of days in 'dates' be equal to the number of element in the 
  # time series corresponding to the 'st.name' station
  if ( ( class(dates) == "Date") & (length(dates) != nrow(x) ) ) 
     stop("Invalid argument: 'length(dates)' must be equal to 'nrow(x)'")  
     
  # Amount of stations in 'x'
  nstations <- ncol(x)

  # ID of all the stations in 'x'
  snames <- colnames(x)  
  
  # Computing the Starting and Ending Year of the analysis
  Starting.Year <- as.numeric(format(range(dates)[1], "%Y"))
  Ending.Year   <- as.numeric(format(range(dates)[2], "%Y"))
  
  # Amount of Years belonging to the desired period
  nyears <- Ending.Year - Starting.Year + 1
  
  # Computing the amount of months with data within the desired period
  ndays   <- length(dates) # number of days in the period
  tmp     <- vector2zoo(rep(0,ndays), dates) 
  tmp     <- daily2monthly.default(x= tmp, FUN=FUN, na.rm=na.rm)
  nmonths <- length(tmp)
  
  # Requiring the Zoo Library (Z’s ordered observations)
  require(zoo)
  
  if (verbose) print("Starting the computations...", quote=FALSE ) 
  
  
  if (out.type=="data.frame") {
  
        # Creating a vector with the names of the field that will be used for storing the results
        field.names <- snames

        # Creating the data.frame that will store the computed averages for each subcatchment
        z <- as.data.frame(matrix(data = NA, nrow = nmonths, ncol = nstations, 
                           byrow = TRUE, dimnames = NULL) )
        
	    colnames(z) <- field.names 
		#rownames(z) <- format(as.Date(mip(from=range(dates)[1], to=range(dates)[2])), "%m-%Y")
        rownames(z) <- format(time(tmp), "%m-%Y")
        
        y = x
        
        for (j in 1:nstations) {
            
            if (verbose) print( paste("Station: ", format(snames[j], width=10, justify="left"),
                                      " : ",format(j, width=3, justify="left"), "/", 
                                      nstations, " => ", 
                                      format(round(100*j/nstations,2), width=6, justify="left"), 
                                      "%", sep=""), quote=FALSE )
                            
            # Transforming the column of 'x' into a zoo object, 
    		# using the dates provided by the user
    		tmp <- vector2zoo(x=y[,j], dates=dates, date.fmt=date.fmt)
    		
    		# Computing the monthly values
    		m <- daily2monthly.default(x= tmp, FUN=FUN, na.rm=na.rm)
    		    
    		if (out.fmt == "numeric") {
    		        z[,j] <- as.numeric(m)
            } else if (out.fmt == "zoo") {
    		              z[,j] <- m
              } # IF end          
                              
        } # FOR end
		
  } else  if (out.type=="db") {
  
        # Creating a vector with the names of the field that will be used for storing the results
        field.names <- c("StationID", "Year", "Month", "Value" )  

        # Creating the data.frame that will store the computed averages for each subcatchment
        z <- as.data.frame(matrix(data = NA, nrow = nmonths*nstations, ncol = 4, 
                            byrow = TRUE, dimnames = NULL) )
        colnames(z) <- field.names      
        
        y = x
        
        for (j in 1:nstations) {
            
            if (verbose) print( paste("Station: ", format(snames[j], width=10, justify="left"),
                                      " : ",format(j, width=3, justify="left"), "/", 
                                      nstations, " => ", 
                                      format(round(100*j/nstations,2), width=6, justify="left"), 
                                      "%", sep=""), quote=FALSE )
                            
        # Transforming the column of 'x' into a zoo object, 
		    # using the dates provided by the user
		    tmp <- vector2zoo(x=x[,j], dates=dates, date.fmt=date.fmt)
		
		    # Computing the monthly values
		    m <- daily2monthly.default(x= tmp, FUN=FUN, na.rm=na.rm)
            
		    if (out.fmt == "numeric") {
    		        m <- as.numeric(m)
            } # IF end
			
		# Putting the annual/monthly values in the output data.frame
        # The first column of 'x' corresponds to the Year
        row.ini <- (j-1)*nmonths + 1
        row.fin <-  j*nmonths
            
        z[row.ini:row.fin, 1] <- snames[j] # it is automatically repeated 'nmonths' times
        z[row.ini:row.fin, 2] <- format(as.Date(time(m)), "%Y")
        z[row.ini:row.fin, 3] <- format(as.Date(time(m)), "%b")              
        z[row.ini:row.fin, 4] <- m
                              
        } # FOR end
  }         
   
  return( z )
  
 } #'daily2monthly.data.frame' END
 
 
daily2monthly.matrix  <- function(x, FUN, na.rm=TRUE,
                                  dates, date.fmt="%Y-%m-%d", 
				  out.type="data.frame", 
				  out.fmt="numeric",
                                  verbose=TRUE,...) {
                  
   x <- as.data.frame(x)   
   #NextMethod("daily2annual")  # I don't know why is redirecting to 'daily2monthly.default' instead of 'daily2monthly.data.frame'....
   daily2monthly.data.frame(x=x, FUN=FUN, na.rm=na.rm,
                            dates=dates, date.fmt=date.fmt, 
			    out.type=out.type, 
			    out.fmt=out.fmt,
                            verbose=verbose,...)
                  
} # 'daily2monthly.matrix  ' END



#####################################
#          daily2annual             #
#####################################
# Generic function for transforming a DAILY regular time series into an ANNUAL one

# 'x'      : Daily zoo object which values will be converted into annual one.
# 'FUN'    : Function that have to be applied for transforming from Daily to Annual time step
#            For Precipitation FUN MUST be 'sum'
#            For Temperature and Flow time series, FUN MUST be 'mean'
# 'na.rm'  : TRUE : the annual mean  value is computed considering only those values different from NA
#            FALSE: if there is AT LEAST one NA within a year, the monthly mean value is NA 
# 'out.fmt': character indicating the format for the output time series. Possible values are:
#            -) "%Y"      : only the year will be used for the time. Default option. (e.g., "1961" "1962"...)
#            -) "%Y-%m-%d": a complete date format will be used for the time. Default option. (e.g., "1961" "1962"...)

daily2annual <-function(x, ...) UseMethod("daily2annual")

daily2annual.default <- function(x, FUN, na.rm=TRUE, out.fmt="%Y",...) { 

	 # Checking that 'x' is a zoo object
	 if (is.na(match(class(x), c("zoo")))) 
			stop("Invalid argument: 'x' must be of class 'zoo'") 
			
	 # Checking that the user provide a valid value for 'FUN'
	 if (missing(FUN)) {
		 stop("Missing argument value: 'FUN' must contain a valid function for aggregating the daily values") }
         
     # Checking that 'x' is a zoo object
	 if ( is.na(match(out.fmt, c("%Y", "%Y-%m-%d") ) ) ) 
			stop("Invalid argument: 'out.fmt' must be in c('%Y', '%Y-%m-%d')" ) 
			
	 # Requiring the Zoo Library (Z’s ordered observations)
	 require(zoo) 
	 
	 # Checking the user provide a valid value for 'x'
	 if (is.na(match(sfreq(x), c("daily", "monthly")))) {
		 stop(paste("Invalid argument: 'x' is not a daily or mothly ts, it is a ", sfreq(x), " ts", sep="") ) }
			  
	 # Generating an Annual time series of Total Annual Precipitation (Annual sum of daily values)
	 tmp <- aggregate(x, by=format( time(x), "%Y" ), FUN, na.rm= na.rm)
	 
	 # Getting the position of all the years in which there were no values
     # mean(NA:NA, na.rm=TRUE) == NaN
     nan.index <- which(is.nan(tmp))
  
     # Changing all the NaN's by NA's
     if ( length(nan.index) > 0 ) { tmp[nan.index] <- NA }

	 # If the user wants a complete data format for the output annual series:
	 if (out.fmt == "%Y-%m-%d") { 
	   time(tmp) <- as.Date(paste( time(tmp), "-01-01", sep=""))
	 } # IF END
	 
	 return(tmp)

} # 'daily2annual.default' end



# 'dates'   : "numeric", "factor", "Date" indicating how to obtain the 
#             dates for correponding to the 'sname' station
#             If 'dates' is a number, it indicates the index of the column in 
#                'x' that stores the dates
#             If 'dates' is a factor, it have to be converted into 'Date' class, 
#                using the date format  specified by 'date.fmt'
#             If 'dates' is already of Date class, the following line verifies that
#                the number of days in 'dates' be equal to the number of element in the 
#                time series corresponding to the 'st.name' station
# 'date.fmt': character indicating the format in which the dates are stored in 'dates'.
#             ONLY required when class(dates)=="factor" or "numeric"
# 'out.type': string that define the desired type of output. Possible values are
#             -) "data.frame": a data.frame, with as many columns as stations 
#                              are included in 'x', and an additional column indicating the Year
#             -) "db"        : a data.frame, with 3 colums will be produced.
#                              The first column will store the Year, 
#                              The second column will store the ID of the station,
#                              The third column will contain the seasonal 
#                                value corresponding to that year and that station.
# 'verbose' : logical; if TRUE, progress messages are printed 
daily2annual.data.frame <- function(x, FUN, na.rm=TRUE, out.fmt="%Y",
                                    dates, date.fmt="%Y-%m-%d", 
								    out.type="data.frame", 
									verbose=TRUE,...) {
  # Checking that the user provide a valid value for 'FUN'
  if (missing(FUN)) {
		 stop("Missing argument value: 'FUN' must contain a valid function for aggregating the daily values") }
		 
  # Checking that the user provied a valid argument for 'out.type'  
  if (is.na(match( out.type, c("data.frame", "db") ) ) ) 
      stop("Invalid argument: 'out.type' must be in c('data.frame', 'db'")
	  
  # Checking that the user provied a valid argument for 'dates'       
  if (missing(dates)) {
      stop("Missing argument: 'dates' must be provided") 
  } else      
     # Checking that the user provied a valid argument for 'dates'  
     if (is.na(match(class(dates), c("numeric", "factor", "Date")))) 
         stop("Invalid argument: 'dates' must be of class 'numeric', 'factor', 'Date'")
        
  # If 'dates' is a number, it indicates the index of the column of 'x' that stores the dates
  # The column with dates is then substracted form 'x' for easening the further computations
  if ( class(dates) == "numeric" ) {    
    tmp   <- dates
    dates <- as.Date(x[, dates], format= date.fmt)
    x     <- x[-tmp]
  }  # IF end 
  
  # If 'dates' is a factor, it have to be converted into 'Date' class, 
  # using the date format  specified by 'date.fmt'
  if ( class(dates) == "factor" ) dates <- as.Date(dates, format= date.fmt)
  
  # If 'dates' is already of Date class, the following line verifies that
  # the number of days in 'dates' be equal to the number of element in the 
  # time series corresponding to the 'st.name' station
  if ( ( class(dates) == "Date") & (length(dates) != nrow(x) ) ) 
     stop("Invalid argument: 'length(dates)' must be equal to 'nrow(x)'")  
     
  # Amount of stations in 'x'
  nstations <- ncol(x)

  # ID of all the stations in 'x'
  snames <- colnames(x)  
  
  # Computing the Starting and Ending Year of the analysis
  Starting.Year <- as.numeric(format(range(dates)[1], "%Y"))
  Ending.Year   <- as.numeric(format(range(dates)[2], "%Y"))
  
  # Amount of Years belonging to the desired period
  #nyears <- Ending.Year - Starting.Year + 1
  
  # Computing the amount of years with data within 'x'
  ndays    <- length(dates) # number of days in the period
  tmp      <- vector2zoo(rep(0, ndays), dates) 
  tmp      <- daily2annual.default(x= tmp, FUN=FUN, na.rm=na.rm)  
  nyears   <- length(tmp) #number of years in the period
  
  # Generating a string vector with the years effectively within 'x'
  if (out.fmt == "%Y") {
    chryears <- time(tmp) 
  } else chryears <- format(time(tmp), "%Y")
  
  
  # Requiring the Zoo Library (Z’s ordered observations)
  require(zoo)
  
  if (verbose) print("Starting the computations...", quote=FALSE )
  
  if (out.type == "data.frame") {
  
	# Vector with the names of the field that will be used for storing the results
	field.names <- snames

	# Creating the data.frame that will store the computed averages for each station
	z <- as.data.frame(matrix(data = NA, nrow = nyears, ncol = nstations, 
						byrow = TRUE, dimnames = NULL) )
						
	colnames(z) <- field.names	
	rownames(z) <- chryears        
	  
	z[1:nstations] <- sapply(1:nstations, function(j,y) {

		
		if (verbose) print( paste("Station: ", format(snames[j], width=10, justify="left"),
					              " : ",format(j, width=3, justify="left"), "/", 
					              nstations, " => ", 
					              format(round(100*j/nstations,2), width=6, justify="left"), 
					              "%", sep=""), quote=FALSE )
						
		# Transforming the column of 'x' into a zoo object, 
		# using the dates provided by the user
		tmp <- vector2zoo(x=x[,j], dates=dates, date.fmt=date.fmt)
		
		# Computing the annual values
		z[,j] <- daily2annual.default(x= tmp, FUN=FUN, na.rm=na.rm)
						 
	}, y = x) # sapply END
  
  } else if (out.type == "db") {
  
        # Creating a vector with the names of the field that will be used for storing the results
        field.names <- c("StationID", "Year", "Value" )  

        # Creating the data.frame that will store the computed averages for each subcatchment
        z <- as.data.frame(matrix(data = NA, nrow = nyears*nstations, ncol = 3, 
                           byrow = TRUE, dimnames = NULL) )
        colnames(z) <- field.names      
        
        y = x
        
        for (j in 1:nstations) {
            
            if (verbose) print( paste("Station: ", format(snames[j], width=10, justify="left"),
                                      " : ",format(j, width=3, justify="left"), "/", 
                                      nstations, " => ", 
                                      format(round(100*j/nstations,2), width=6, justify="left"), 
                                      "%", sep=""), quote=FALSE )
                            
            # Transforming the column of 'x' into a zoo object, 
		    # using the dates provided by the user
		    tmp <- vector2zoo(x=x[,j], dates=dates, date.fmt=date.fmt)
		
		    # Computing the annual values
		    a <- daily2annual.default(x= tmp, FUN=FUN, na.rm=na.rm, out.fmt="%Y-%m-%d")
			
			# Putting the annual/monthly values in the output data.frame
            # The first column of 'x' corresponds to the Year
            row.ini <- (j-1)*nyears + 1
            row.fin <-  j*nyears
            
            z[row.ini:row.fin, 1] <- snames[j] # it is automatically repeted 'nmonths' times
            z[row.ini:row.fin, 2] <- format(as.Date(time(a)), "%Y")              
            z[row.ini:row.fin, 3] <- a
                              
        } # FOR end
  
    } # IF end
   
  return( z )
  
 } #'daily2annual.data.frame' END
 
 
daily2annual.matrix  <- function(x, FUN, na.rm=TRUE, out.fmt="%Y",
                                 dates, date.fmt="%Y-%m-%d", 
								 out.type="data.frame", 
								 verbose=TRUE,...) {
                  
   x <- as.data.frame(x)    
   #NextMethod("daily2annual")
   daily2annual.data.frame(x=x, FUN=FUN, na.rm=na.rm, 
                           out.fmt=out.fmt,
                           dates=dates, date.fmt=date.fmt, 
                           out.type=out.type, 
                           verbose=verbose,...)
                  
} # 'daily2annual.matrix  ' END


#######################################
#          monthly2annual             #
#######################################
#           May 13th, 2009            #
#######################################
# This function transform a MONTHLY regular time series into an ANNUAL one .
# It is only a wrapper to 'daily2annual' and it was not merged as unique function just for clarity

# 'x'      : Monthly zoo object that will be converted into annual ones.
#            class(x) must be 'zoo'
# 'FUN'    : Function that have to be applied for transforming from daily to monthly time step
#            For precipitation FUN MUST be "sum"
#            For temperature and flow time series, FUN MUST be "mean"
# 'na.rm'  : Logical. Should missing values be removed?
#            TRUE : the monthly and annual values  are computed considering only those values different from NA
#            FALSE: if there is AT LEAST one NA within a year, the monthly and annual values are NA 
# 'out.fmt': date format for the output time series. Possible values are:
#            -) "%Y"      : only the year will be used for the time. Default option. (e.g., "1961" "1962"...)
#            -) "%Y-%m-%d": a complete date format will be used for the time. Default option. (e.g., "1961" "1962"...)

monthly2annual <-function(x,...) UseMethod("daily2annual")

#monthly2annual.default <- function(x, FUN, na.rm=TRUE, out.fmt="%Y",...) { 

   
#daily2annual.default(x, FUN, na.rm=na.rm, out.fmt=out.fmt )

#} # 'monthly2annual.default' end


#monthly2annual.data.frame <- function(x, FUN, na.rm=TRUE, out.fmt="%Y",
                                      #dates, date.fmt="%Y-%m-%d", 
				      #out.type="data.frame", 
				      #verbose=TRUE,...) {
                                      
#daily2annual.data.frame(x, FUN, na.rm=na.rm, out.fmt=out.fmt, 
                          #dates=dates, date.fmt=date.fmt, out.type=out.type, 
                          #verbose=verbose )
                                      
#} # 'monthly2annual.data.frame' END



#######################################################
# .hydroplotts Daily, Monthly and Annual Time Series  #
#######################################################
# It requires the function'drawxaxis' that is stored in the 'lib_Plot.R' library
# 'x'		 : daily time series of type 'zoo'
# 'x.monthly : monthly time series of type 'zoo'
# 'x.annual' : annual time series of type 'zoo'
# 'win.len1' : number of days for being used in the computation of the first moving average
# 'win.len2' : number of days for being used in the computation of the second moving average
# 'var.type' : string representing the type of variable being plotted (e.g., "Precipitation", "Temperature" or "Flow"). 
#              ONLY used for labelling the y axis and the title of the plot (in case it is missing) 
# 'var.unit' : string representing the measurement unit of the variable being plotted ("mm" for precipitation, "C" for temperature, and "m3/s" for flow). 
#              ONLY used for labelling the y axis
# 'main'     : string representing the main title of the plot
# 'pfreq'    : string indicating how many plots are desired by the user. 
#              Valid values are:                     
#              -) 'dma': Daily, Monthly and Annual values are plotted 
#              -) 'ma' : Monthly and Annual values are plotted 
#              -) 'dm' : Daily and Monthly values are plotted 
# 'tick.tstep': string indicating the time step that have to be used for 
#               putting the ticks ont he time axis. 
#               Possible values are: 'days', 'months', 'years' 
# 'lab.tstep' : string indicating the time step that have to be used for 
#               putting the labels ont he time axis. 

.hydroplotts <- function(x, x.monthly, x.annual, win.len1, win.len2, 
                         var.type="Precipitation", var.unit="mm", main,
			 pfreq="dma", tick.tstep= "months", lab.tstep= "years" ) {
                
      # Checking that 'x' is a zoo object
      if (is.na(match(class(x), c("zoo")))) 
            stop("Invalid argument: 'x' must be of class 'zoo'")
            
      # Checking that 'x.monthly' is a zoo object
      if (is.na(match(class(x.monthly), c("zoo")))) 
            stop("Invalid argument: 'x.monthly' must be of class 'zoo'")
            
      # Checking that 'x.annual' is a zoo object
      if (is.na(match(class(x.annual), c("zoo")))) 
            stop("Invalid argument: 'x.annual' must be of class 'zoo'")
             
      # Checking that the user provied a valid argument for 'pfreq' 
      if (is.na(match(pfreq, c("dma", "ma", "dm")))) 
          stop("Invalid argument: 'pfreq' must be in c('dma', 'ma', 'dm')")
          
      # Checking that the user provided a valid argument for 'tick.tstep'       
      if (is.na(match(tick.tstep, c("days", "months", "years") ) ) ) 
         stop("Invalid argument: 'tick.tstep' must be in c('days', 'months', 'years')")
         
      # Checking that the user provided a valid argument for 'lab.tstep'       
      if (is.na(match(lab.tstep, c("days", "months", "years") ) ) ) 
         stop("Invalid argument: 'lab.tstep' must be in c('days', 'months', 'years')")
            
      # If the user did not provide a title for the plots, this is created automatically
      if (missing(main)) { main= paste(var.type, var.unit, sep=" ") }
      
      # Requiring the Zoo Library (Zoo's ordered observations)
      require(zoo)
      
      # Booleans indicating if the moving averages for the dayly and monthly 
      # time series can be computed and ploted. By default they  are FALSE, 
      # and only if the lenght(x) is large enough they are changed into TRUE
      d.ma1 <- FALSE
      d.ma2 <- FALSE
      m.ma1 <- FALSE
      m.ma2 <- FALSE
      
      # Generating a Moving Average of the Daily time series, with a window width 'win.len1'
      if (win.len1 > 0 ) { # If 'win.len1==0', the moving average is not computed
          win.len <- win.len1  
          if (length(x) >= win.len) { 
            d.ma1 <- TRUE
            daily.ma1 <- ma.zoo(x, win.len) }
      } # IF end
      
      # Generating a Moving Average of the Daily time series, with a window width 'win.len2'
      if (win.len2 > 0 ) {  # If 'win.len2==0', the moving average is not computed
          win.len <- win.len2
          if (length(x) >= win.len) { 
            d.ma2 <- TRUE
            daily.ma2 <- ma.zoo(x, win.len) }
      } # IF end
     
      # Generating a Moving Average of the Monthly time series, with a window width 'win.len1'
      win.len <- round(win.len1/365,1)*12
      if (length(x.monthly) >= win.len) { 
        m.ma1 <- TRUE
        monthly.ma1 <- ma.zoo( x.monthly, win.len )  }
        
      # Generating a Moving Average of the Monthly time series, with a window width 'win.len2'
      win.len <- round(win.len2/365,1)*12
      if (length(x.monthly) >= win.len) { 
        m.ma2 <- TRUE
        monthly.ma2 <- ma.zoo( x.monthly, win.len ) }    
     
      # Checking if the Daily ts have to be plotted 
      if ( pfreq %in% c("dma", "dm") ) { 
          # Generating the title of the Daily Time Series plot
          title <- paste("Daily", main, sep= " ")
          # Plotting the Daily Time Series
          # xaxt = "n": is for avoiding drawing the x axis
          plot.zoo(x, xaxt = "n", type="o", lwd=1, lty=1, col="blue", cex = .5, 
                  main=title, xlab="Time", ylab=paste(var.type," [", var.unit,"/day]", sep="") )
          # Draws monthly ticks in the X axis, but labels only in years
          drawxaxis(x, tick.tstep=tick.tstep, lab.tstep=lab.tstep) 
         
          if (d.ma1) {
            # Plotting the 1st Moving Average of the Daily time series. If win.len1=365*1 => "Annual Moving Average"
            lines(daily.ma1, type="o", lty=2, lwd=1, col="green", cex = .5) } 
          if (d.ma2) {
            # Plotting the 2nd Moving Average of the Daily time series. If win.len2=365*3 => "Moving Average of 3 Years"
            lines(daily.ma2, type="o", lty=3, lwd=1, col="red", cex = .5) }
          # Drawing a legend. y.intersp=0.5, is for the vertical spacin in the legend
          legend("topleft", c("Daily series", paste("MA(", round(win.len1/365,2), " years)", sep=""), 
                 paste("MA(", round(win.len2/365,2), " years)", sep="") ), y.intersp=0.5,
                 bty="n", cex =0.9, col = c("blue","green","red"), lwd= c(1,1,1), lty=c(1,2,3) ) #bty="n" => no box around the legend 
      } # IF end          
     
     
      # Checking if the Monthly ts have to be plotted 
      if ( pfreq %in% c("dma", "dm", "ma") ) { 
        # Generating the title of the Monthly Time Series plot
        title <- paste("Monthly", main, sep= " ")     
        # Plotting the Monthly time series
        plot.zoo(x.monthly, xaxt = "n", type="o", lwd=1, col="blue", cex = .5, 
                 main=title, xlab="Time", ylab=paste(var.type," [", var.unit,"/month]", sep="") )
               
        
        # Draws monthly ticks in the X axis, but labels only in years
        drawxaxis(x.monthly, tick.tstep=tick.tstep, lab.tstep=lab.tstep)
        if (m.ma1) {
        # Plotting the 1st Moving Average of the Daily time series. If win.len1=365*1 => "Annual Moving Average"
        lines(monthly.ma1, type="o", lty=2, lwd=1, col="green", cex = .5) } 
        if (m.ma2) {
        # Plotting the 2nd Moving Average of the Daily time series. If win.len2=365*3 => "Moving Average of 3 Years"
        lines(monthly.ma2, type="o", lty=3, lwd=1, col="red", cex = .5) } 
        # Drawing a legend
        legend("topleft", c("Monthly series", paste("MA(", round(win.len1/365,1), " years)", sep=""), 
             paste("MA(", round(win.len2/365,1), " years)", sep="") ), y.intersp=0.5, 
             bty="n", cex =0.9, col = c("blue","green","red"), lwd= c(1,1,1), lty=c(1,2,3) ) #bty="n" => no box around the legend    
      } # IF end    

      # Checking if the Annual ts have to be plotted 
      if ( pfreq %in% c("dma", "ma") ) { 
          # Generating the title of the Annual Time Series plot
           title <- paste("Annual", main, sep= " ")
          # Plotting the Annual time series
          plot.zoo(x.annual, xaxt = "n", type="o", lwd=1, col="blue", cex = .5, 
                   main=title, xlab="Time", ylab=paste(var.type," [", var.unit,"/year]", sep="") )
          # Draws monthly ticks in the X axis, but labels only in years
          drawxaxis(x.annual, tick.tstep=tick.tstep, lab.tstep=lab.tstep)  
      } # IF end
     
} # '.hydroplotts' end



#####################################################
# BoxPlot of Daily, Monthly and Annual Time Serires #
#####################################################
# 'x'		 : daily time series of type 'zoo'
# 'x.monthly : monthly time series of type 'zoo'
# 'x.annual' : annual time series of type 'zoo'
# 'var.type' : string representing the type of variable being plotted 
#              (e.g., "Precipitation", "Temperature" or "Flow"). 
#              ONLY used for labelling the y axis and the title of the plot (in case it is missing) 
# 'var.unit' : string representing the measurement unit of the variable 
#              being plotted ("mm" for precipitation, "C" for temperature, and "m3/s" for flow). 
#              ONLY used for labelling the y axis
# 'main'     : string representing the main title of the plot
# 'pfreq'    : string indicating how many plots are desired by the user. 
#              Valid values are:                     
#              -) 'dma': Daily, Monthly and Annual values are plotted 
#              -) 'dm' : Daily and Monthly values are plotted 
#              -) 'ma' : Monthly and Annual values are plotted 
#              
.hydroplotboxplot <- function(x, x.monthly, x.annual, 
                              var.type="Precipitation", var.unit="mm", main,
			      pfreq="dma") {
                              
  # Checking that 'x' is a zoo object
  if (is.na(match(class(x), c("zoo")))) 
        stop("Invalid argument: 'x' must be of class 'zoo'")
        
  # Checking that 'x.monthly' is a zoo object
  if (is.na(match(class(x.monthly), c("zoo")))) 
        stop("Invalid argument: 'x.monthly' must be of class 'zoo'")
        
  # Checking that 'x.annual' is a zoo object
  if (is.na(match(class(x.annual), c("zoo")))) 
        stop("Invalid argument: 'x.annual' must be of class 'zoo'")
		 
  # Checking that the user provied a valid argument for 'pfreq' 
  if (is.na(match(pfreq, c("dma", "ma", "dm")))) 
      stop("Invalid argument: 'pfreq' must be in c('dma', 'ma', 'dm')")
 
 # Requiring the Zoo Library (Z’s ordered observations)
 require(zoo)
 
 # Checking if the Daily ts have to be plotted 
 if ( pfreq %in% c("dma", "dm") ) { 
   # Generating a factor based on the year in which each daily date falls
   cyear <- format(time(x), "%Y")
   years <- factor(cyear, levels=unique(cyear), ordered=TRUE)
   # Generating the title of the Daily plot
   title <- paste("Boxplot of Daily", main, sep= " ")
   # Drawing boxplot of Daily values against Year
   boxplot( coredata(x)~years, main=title, ylab=paste(var.type," [", var.unit, "/day]", sep=""), col="lightblue")
 } # IF end
 
 # Checking if the Monthly ts have to be plotted 
 if ( pfreq %in% c("dma", "dm", "ma") ) { 
   # Generating a factor based on the month in which each monthly date falls
   cmonth <- format(time(x.monthly), "%b")
   months <- factor(cmonth, levels=unique(cmonth), ordered=TRUE)
   # Generating the title of the Monthly plot
   title <- paste("Boxplot of Monthly", main, sep= " ")
   # Drawing boxplot of Monthly values against Year
   boxplot( coredata(x.monthly)~months, main=title, ylab=paste(var.type," [", var.unit,"/month]", sep=""), col="lightblue")
 } # IF end
 
 # Checking if the Annual ts have to be plotted 
 if ( pfreq %in% c("dma", "ma") ) { 
   # Generating the title of the Annual plot
   title <- paste("Boxplot of Annual", main, sep= " ")
   # Drawing boxplot of Annual values against Year
   boxplot( coredata(x.annual), main=title, ylab=paste(var.type," [", var.unit, "/year]", sep=""), col="lightblue") 
 } # IF end
 
} #'.hydroplotboxplot' end


#######################################################
# Histogram of Daily, Monthly and Annual Time Serires #
#######################################################
# 'x'		 : daily time series of type 'zoo'
# 'x.monthly : monthly time series of type 'zoo'
# 'x.annual' : annual time series of type 'zoo'
# 'var.type' : string representing the type of variable being plotted (e.g., "Precipitation", "Temperature" or "Flow"). 
#              ONLY used for labelling the y axis and the title of the plot (in case it is missing) 
# 'var.unit' : string representing the measurement unit of the variable being plotted ("mm" for precipitation, "C" for temperature, and "m3/s" for flow). 
#              ONLY used for labelling the x axis
# 'main'     : string representing the main title of the plot
# 'pfreq'    : string indicating how many plots are desired by the user. 
#              Valid values are:                     
#              -) 'dma': Daily, Monthly and Annual values are plotted 
#              -) 'ma' : Monthly and Annual values are plotted 
#              -) 'dm' : Daily and Monthly values are plotted 
.hydroplothist <- function(x, x.monthly, x.annual, 
                           var.type="Precipitation", var.unit="mm", main,
			   pfreq="dma") {
                                
      # Checking that 'x' is a zoo object
      if (is.na(match(class(x), c("zoo")))) 
            stop("Invalid argument: 'x' must be of class 'zoo'")
            
      # Checking that 'x.monthly' is a zoo object
      if (is.na(match(class(x.monthly), c("zoo")))) 
            stop("Invalid argument: 'x.monthly' must be of class 'zoo'")
            
      # Checking that 'x.annual' is a zoo object
      if (is.na(match(class(x.annual), c("zoo")))) 
            stop("Invalid argument: 'x.annual' must be of class 'zoo'")
             
      # Checking that the user provied a valid argument for 'pfreq' 
      if (is.na(match(pfreq, c("dma", "ma", "dm")))) 
          stop("Invalid argument: 'pfreq' must be in c('dma', 'ma', 'dm')")
     
     # Requiring the Zoo Library (Z’s ordered observations)
     require(zoo)
     
     # Checking if the Daily ts have to be plotted 
     if ( pfreq %in% c("dma", "dm") ) { 
       # Generating the title of the Daily plot
       title <- paste("Histogram of Monthly", main, sep= " ")
       # Drawing an histogram of Daily Precipitation
       hist(x, br=100, freq=FALSE, main=title, xlab=paste(var.type," [", var.unit, "/day]", sep=""), ylab="Pbb", col="lightblue")
     } # IF end
     
     # Checking if the Monthly ts have to be plotted 
     if ( pfreq %in% c("dma", "dm", "ma") ) { 
       # Generating the title of the Monthly plot
       title <- paste("Histogram of Monthly", main, sep= " ")
       # Drawing an histogram of Monthly Precipitation
       hist(x.monthly, br=10, freq=FALSE, main=title, xlab=paste(var.type," [", var.unit, "/month]", sep=""), ylab="Pbb", col="lightblue")
     } # IF end
     
     # Checking if the Annual ts have to be plotted 
     if ( pfreq %in% c("dma", "ma") ) { 
       # Generating the title of the Annual plot
       title <- paste("Histogram of Annual", main, sep= " ")
       # Drawing an histogram of Annual Precipitation
       hist(x.annual, br=5, freq=FALSE, main=title, xlab=paste(var.type," [", var.unit, "/year]", sep=""), ylab="Pbb", col="lightblue")
     } # IF end
 
} # '.hydroplothist' end


#########################################################################
# hydroplot: Daily, Monthly and Annual plots of hydrological time series#
#########################################################################
# 9 plots:
# 1: Line plot with Daily time series, with 2 moving averages, specified by 'win.len1' and 'win.len2'
# 2: Line plot with Monthly time series, with 2 moving averages, specified by 'win.len1' and 'win.len2'
# 3: Line plot with Annual time series
# 4: Boxplot with daily time series
# 5: Boxplot with monthly time series
# 6: Boxplot with annual time series
# 7: Histogram of the daily time series
# 8: Histogram of the monthly time series
# 9: Histogram of the annual time series
    
# 'x'	     : variable of type 'zoo', with daily ts
# 'same'     : Character representing the name of the meteorological station
#              ONLY used for labelling the title
# 'var.unit' : Character representing the measurement unit of the variable being plotted. 
#              ONLY used for labelling the axes
#              e.g., "mm" for precipitation, "C" for temperature, and "m3/s" for flow. 
# 'main'     : Character representing the main title of the plot. If the user did not provide a title, this is
#              created automatically as: main= paste(var.type, "at", st.name, sep=" "), 
# 'win.len1' : number of days for being used in the computation of the first moving average
# 'win.len2' : number of days for being used in the computation of the second moving average
# 'ptype'    : Character indicating the type of plot that will be plotted. Valid values are
#              -) ptype= "ts"              => only time series
#              -) ptype= "ts+boxplot"      => only time series + boxplot               
#              -) ptype= "ts+hist"         => only time series + histogram 
#              -) ptype= "ts+boxplot+hist" => time series + boxplot + histogram    
# 'pfreq'    : Character indicating how many plots are desired by the user. 
#              Valid values are:                     
#              -) 'dma': Daily, Monthly and Annual values are plotted 
#              -) 'ma' : Monthly and Annual values are plotted 
#              -) 'dm' : Daily and Monthly values are plotted  
# 'var.type' : string representing the type of variable being plotted 
#              Used for determining the function used for computing the 
#              Monthly and Annual values when 'FUN' is missing
#              Valid values are:
#              -) "Precipitation" => FUN = sum
#              -) "Temperature"   => FUN = mean
#              -) "Flow"          => FUN = mean 
# 'FUN'      : ONLY required when 'var.type' is missing
#              Function that have to be applied for transforming from daily to monthly or annual time step
#              For precipitation FUN MUST be "sum"
#              For temperature and flow time series, FUN MUST be "mean"#             
# 'na.rm'    : Logical. Should missing values be removed?
#              TRUE : the monthly and annual values  are computed considering only those values different from NA
#              FALSE: if there is AT LEAST one NA within a year, the monthly and annual values are NA
# 'tick.tstep': string indicating the time step that have to be used for 
#               putting the ticks ont he time axis. 
#               Possible values are: 'days', 'months', 'years' 
# 'lab.tstep' : string indicating the time step that have to be used for 
#               putting the labels ont he time axis.  
 
hydroplot <- function(x, sname="X", 
                      var.type="", 
                      FUN, 
                      na.rm=TRUE, 
                      var.unit="mm", 
                      main, 
                      win.len1=365*1, 
                      win.len2=365*3,
                      ptype="ts+boxplot+hist",
		      pfreq="dma",
                      tick.tstep= "months", 
                      lab.tstep= "years") {
                       
     # Checking that 'x' is a zoo object
     if (is.na(match(class(x), c("zoo")))) 
            stop("Invalid argument: 'x' must be of class 'zoo'")
             
     # Checking that the user provied a valid argument for 'ptype' 
     if (is.na(match(ptype, c("ts", "ts+boxplot", "ts+hist", "ts+boxplot+hist")))) 
            stop("Invalid argument: 'ptype' must be in c('ts', 'ts+boxplot', 'ts+hist', 'ts+boxplot+hist')")
            
     # Checking that the user provied a valid argument for 'pfreq' 
     if ( sfreq(x) == "daily" ) {
       if (is.na(match(pfreq, c("dma", "dm", "ma")))) 
          stop("Invalid argument: 'pfreq' must be in c('dma', 'ma', 'dm')")
     } else if ( sfreq(x) == "monthly" ) {
         if ( pfreq != "ma" ) 
            print("Warning: 'x' is a monthly object, so 'pfreq' has been changed to 'ma'")
            pfreq="ma"
       } # ELSE end
            
     # If the user did not provide a title for the plots, this is created automatically
     if (missing(main)) {
            main= paste(var.type, "at", sname, sep=" ") }
            
     # Checking that the user provied a valid argument for 'var.type'       
     if (missing(FUN)) {
        # If the user did not provide a title for the plots, this is created automatically
        if (missing(var.type)) { 
          stop("Missing argument: 'var.type' OR 'FUN' must be provided") 
        } else # If 'var.type' is provided
             # Checking that the user provied a valid argument for 'var.type'       
             if (is.na(match(var.type, c("Precipitation", "Temperature", "Flow") ) ) ) {
                   stop("Invalid argument: 'var.type' must be in c('Precipitation', 'Temperature', 'Flow')") 
             } else {
                      if (var.type=="Precipitation") { 
                          FUN <- sum 
                          if (missing(var.unit)) { var.unit <- "mm"   }
                      } else if (var.type=="Temperature") { 
                               FUN <- mean
                               if (missing(var.unit)) { var.unit <- "dC" }
                             } else if (var.type=="Flow") { 
                                    FUN <- mean
                                    if (missing(var.unit)) { var.unit= "m3/s" }
                               }
                      } #ELSE end 
     } # IF end 
     
     def.par <- par(no.readonly = TRUE) # save default, for resetting... 
     on.exit(par(def.par)) 
            
     # Requiring the Zoo Library
     require(zoo)    
      
     # Generating a Monthly time series   
     if ( sfreq(x) == "daily" ) {         
       x.monthly <- daily2monthly(x, FUN=FUN, na.rm=na.rm)
     } else if ( sfreq(x) == "monthly" ) { 
        x.monthly <- x
        } else x.monthly <- NA
     
     # Generating an Annual time series
     if ( !is.na( match( sfreq(x), c("daily", "monthly") ) ) ) {         
       x.annual <- daily2annual(x, FUN=FUN, na.rm=na.rm, out.fmt="%Y-%m-%d")
     } else if ( sfreq(x) == "annual" ) { 
        x.annual <- x
        } else x.annual <- NA
     
     # If the user didn't provided a value for 'tick.tstep' and 
     # the lenght of the daily ts is less than 1 year, the ticks in 
     # the 'x' axis are placed by day; if larger than a year, they are placed by month
     if ( missing(tick.tstep) ) {
       #if (length(x) <= 366) {
       if ( ( (sfreq(x) == "daily") & ( length(x) <= 366 ) ) |   
          ( (sfreq(x) == "monthly") & ( length(x) <= 12 ) ) ) {
         tick.tstep <- "days"
       } else tick.tstep <- "months"
     } # IF end
        
     # If the user didn't provided a value for 'lab.tstep' and 
     # the lenght of the daily ts is less than 1 year, the labels in 
     # the 'x' axis are placed by month; if larger than a year, they are placed by year
     if ( missing(lab.tstep) ) {
       #if (length(x) <=366) {
       if ( ( (sfreq(x) == "daily") & ( length(x) <= 366 ) ) |   
          ( (sfreq(x) == "monthly") & ( length(x) <= 12 ) ) ) {
         lab.tstep <- "months"
       } else lab.tstep <- "years"
     } # IF end
     
     # If 'x' is too short for plotting annual values, 'pfreq' is automatically changed
     if ( ( (sfreq(x) == "daily") & ( length(x) <= 366 ) ) |   
          ( (sfreq(x) == "monthly") & ( length(x) <= 12 ) ) ) {
         if ( match(pfreq, c("dma", "ma") ) ) {         
           if (pfreq == "dma") pfreq <- "dm"
           if (pfreq == "ma") pfreq <- "m"
           print(paste("warning: your ts is too short for plotting annual time series => 'pfreq'= ", pfreq, sep="") )
         }    
     } # IF end       
       
     if (ptype=="ts") {
       # Setting up the screen with 3 rows and 3 columns
       if (pfreq == "dma") { par(mfcol=c(3,1)) 
       } else if (pfreq %in% c("dm", "ma")) { par(mfcol=c(2,1)) 
         } # ELSE end	   
       # Drawing the daily, monthly and annual time series of the variable against time
       .hydroplotts(x, x.monthly, x.annual, win.len1, win.len2, var.type, var.unit, main, pfreq, tick.tstep= tick.tstep, lab.tstep= lab.tstep)
     }
     
     else if (ptype=="ts+boxplot") {
       # Setting up the screen with 3 rows and 3 columns
       if (pfreq == "dma") { par(mfcol=c(3,2)) 
       } else if (pfreq %in% c("dm", "ma")) { par(mfcol=c(2,2)) 
         } # ELSE end 
       # Drawing the daily, monthly and annual time series of the variable against time
       .hydroplotts(x, x.monthly, x.annual, win.len1, win.len2, var.type, var.unit, main, pfreq, tick.tstep= tick.tstep, lab.tstep= lab.tstep)
       # Drawing a boxplot of the daily, monthly and annual time series of the variable    
       .hydroplotboxplot(x, x.monthly, x.annual, var.type, var.unit, main, pfreq)
     } 
     
     else if (ptype=="ts+hist") {
       # Setting up the screen with 3 rows and 3 columns
       if (pfreq == "dma") { par(mfcol=c(3,2)) 
       } else if (pfreq %in% c("dm", "ma")) { par(mfcol=c(2,2)) 
         } # ELSE end 
       # Drawing the daily, monthly and annual time series of the variable against time
       .hydroplotts(x, x.monthly, x.annual, win.len1, win.len2, var.type, var.unit, main, pfreq, tick.tstep= tick.tstep, lab.tstep= lab.tstep)
       # Drawing an histogram of the daily, monthly and annual time series of the variable
       .hydroplothist(x, x.monthly, x.annual, var.type, var.unit, main, pfreq)
     }
     else if (ptype=="ts+boxplot+hist") {
       # Setting up the screen with 3 rows and 3 columns
       if (pfreq == "dma") { par(mfcol=c(3,3)) 
       } else if (pfreq %in% c("dm", "ma")) { par(mfcol=c(2,3)) 
         } # ELSE end  
       # Drawing the daily, monthly and annual time series of the variable against time
       .hydroplotts(x, x.monthly, x.annual, win.len1, win.len2, var.type, var.unit, main, pfreq, tick.tstep= tick.tstep, lab.tstep= lab.tstep)
       # Drawing a boxplot of the daily, monthly and annual time series of the variable    
       .hydroplotboxplot(x, x.monthly, x.annual, var.type, var.unit, main, pfreq)
       # Drawing an histogram of the daily, monthly and annual time series of the variable
       .hydroplothist(x, x.monthly, x.annual, var.type, var.unit, main, pfreq)
     }

 
 } # 'hydroplot end
        

#####################################################################
# 'diy' : Given a numeric value of a year, it generates:            #
#         1) a vector with all the days (dates) within the year, OR #
#		  2) the amount of days in the year                         #
#####################################################################
# year    : numeric, the year for which the sequence of days will be generated
# out.type: Character indicating the type of result that is given by this function.
#           Valid values are: 
#		    -) type= "seq"  => a vectorial sequence with all the days within the given year
#		    -) type= "nmbr" => the number of days in the vectorial sequence with all the days within the given year
diy <- function(year, out.type="seq") {

   if (is.na(match(out.type, c("seq", "nmbr")))) 
        stop("Invalid argument: 'out.type' must be of class 'seq' or 'nmbr'")   

   # Generating a Daily-regular time series of Dates, 
   # just for being counted as the maximum amount of possible daily data
   vec <- seq( from=as.Date( paste(year,"-01-01", sep="") ), to=as.Date( paste(year,"-12-31", sep="") ), by= "days" )
    
   if (out.type=="seq") return(vec)
   else if (out.type=="nmbr") return ( length(vec) )


} # 'diy' END


##################################################################
# 'dip': Given any starting and ending dates, it generates:      #
#        1) a vector with all the days between the two dates, OR #
#		 2) the amount of days between the two dates             #
##################################################################

# 'from'	: Character indicating the starting date for computing the number of dyas. 
#             It MUST have the date format specified by 'date.fmt'
# 'to'		: Character indicating the ending date for computing the number of dyas. 
#             It MUST have the date format specified by 'date.fmt'
# 'date.fmt': Character indicating the date format in which you provide 'from' and 'to'. (e.g., "%d-%m-%Y")
# 'out.type': Character indicating the type of result that is given by this function.
#             Valid values are: 
#		      -) type= "seq"  => a vectorial sequence with all the days within the given time period
#		      -) type= "nmbr" => the number of days in the vectorial sequence with all the days within the given time period
dip <- function(from, to, date.fmt="%Y-%m-%d", out.type="seq") {

     # Generating an Annual-regular time series of Dates. 
     vec.days <- seq( from=as.Date(from, format=date.fmt), to=as.Date(to, format=date.fmt), by="days" )
     
     if (out.type=="seq") return(vec.days)
     else if (out.type=="nmbr") return ( length(vec.days) )
 
} # 'dip' END



####################################################################
# 'mip': Given any starting and ending dates, it generates:        #
#        1) a vector with all the months between the two dates, OR #
#		 2) the amount of months between the two dates             #
####################################################################

# 'from'	: Starting date for computing the number of dyas. MUST have the date format specified by 'date.fmt'
# 'to'		: Ending date for computing the number of dyas. MUST have the date format specified by 'date.fmt'
# 'date.fmt': Format of the dates (e.g., "%d-%m-%Y")
# out.type  : type of result that is given by this function
#		      -) type= "seq"  => a vectorial sequence with all the months within the given year
#		      -) type= "nmbr" => the number of days in the vectorial sequence with all the months within the given year
mip <- function(from, to, date.fmt="%Y-%m-%d", out.type="seq") {

     # Generating an Annual-regular time series of Dates. 
     vec.months <- seq( from=as.Date(from, format=date.fmt), to=as.Date(to, format=date.fmt), by="months" )
     
     if (out.type=="seq") return(vec.months)
     else if (out.type=="nmbr") return ( length(vec.months) )
 
} # 'mip' END



####################################################################
# 'yip': Given any starting and ending dates, it generates:        #
#        1) a vector of dates with all the years between the two dates, OR  #
#		 2) the amount of years between the two dates              #
####################################################################
# 18-May-2010

# 'from'	: Starting date for computing the number of years. MUST have the date format specified by 'date.fmt'
# 'to'		: Ending date for computing the number of years. MUST have the date format specified by 'date.fmt'
# 'date.fmt': Format of the dates (e.g., "%d-%m-%Y")
# out.type  : type of result that is given by this function
#		      -) type= "seq"  => a vectorial sequence with all the months within the given year
#		      -) type= "nmbr" => the number of days in the vectorial sequence with all the months within the given year
yip <- function(from, to, date.fmt="%Y-%m-%d", out.type="seq") {

     # Generating an Annual-regular time series of Dates. 
     vec.years <- seq( from=as.Date(from, format=date.fmt), to=as.Date(to, format=date.fmt), by="years" )
     
     if (out.type=="seq") return(vec.years)
     else if (out.type=="nmbr") return ( length(vec.years) )
 
} # 'yip' END








############################################################################
#  'infillixy': Fills in ALL the 'NA' values in 'x' with the corrensponding#
#               values in 'sim'.                                           #
#               16-Dic-2008, 04-Sep-2009                                   #
############################################################################
#  
# 'x'   : 'data.frame' or 'matrix' in which some (observed) values are 'NA'
# 'sim' : 'data.frame' or 'matrix', with the same dimension of 'obs', 
#          that contains the values that will be used for filling in 
#          the 'NA' values in 'obs'
# Result: a 'data.frame' or 'matrix', with the same dimension of 'obs',
#         without 'NA' values.

infillxy <-function(x,...) UseMethod("infillxy")

infillxy.default <- function(x, sim,...) {
  
  if (length(x) != length(sim))
    stop("'x' and 'sim' does not have the same dimension !!")
    
  # vector with the index of all the elements in 'x' that are 'NA'
  na.index <- which( is.na(x) )

  # Replacing the 'NA' values in 'filled' by the correponding values in 'sim'
  x[na.index] <- sim[na.index]
  
  return(x)

} # 'infillxy.default' END


infillxy.matrix <- function(x, sim, ...) {

    if ( !identical(dim(x), dim(sim) ) )
      stop("'x' and 'sim' does not have the same dimension !!")
    
    if (is.na(match(class(sim), c("matrix")))) 
        stop("Invalid argument: 'sim' must be of class 'matrix'") 
        
    # Creating a copy of the original observed values
	z <- x   
	  
	z[,1:ncol(z)] <- sapply(1:ncol(z), function(j,y) {
 
		# Putting the monthly values in the output data.frame
		# The first column of 'x' corresponds to the Year
		z[,j] <- infillxy.default(x= y[,j], sim=sim[, j])
						 
	}, y = x) # sapply END
    
 return(z)
    
} # 'infillxy.matrix' END


infillxy.data.frame <- function(x, sim, ...) {

    x   <- as.matrix(x)
    sim <- as.matrix(sim)
    
    NextMethod("infillxy.matrix")
    
} # 'infillxy.data.frame' END
  

  


#############################################################################
# 'interpol1' : Given a data.frame whose columns contains the time series   #
#               (without missing dates) of several gauging stations, it     #
#               interpolates the value at the station "s" on the day "i",   #
#               using all the other gauging stations.                       # 
#               The interpolation method is a modified IDW, where the       #
#               Pearson's product-momenet coefficient of correlation between#
#               the time series of all teh stations is used instead of the  #
#               spatial distance, following the paper                       #       
#               
#               Two methods can be used for carrying out interpolations:    # 
#               1) "cc-normal": normal coefficient correlation method,      #
#                   where all the stations are used                         # 
#               2) "cc-neighs": modified coefficient correlation method,    #
#                   where only the stations with the highest coefficient    #
#                   of correlation with the target station are used         #   
#               The overall performance of this method was better than the  #
#               overall performance of the traditional IDW, considering 146 #
#               stations of temperature and 349 stations of precipitation   #
#               with daily data during 30 years                             # 
#               This was made in a way of carrying out a posterior          #
#               cross-validation for each station (comparing the observed   #
#               with the interpolated values)                               #
#############################################################################
#                 Dec 2008, Jan 2009                #
#####################################################
# 'x.ts.catch': data.frame that contains the time series of all the stations involved in the computations 
#               The name of each column in 'x.ts.catch' have to correspond to the names of the gauging station  
#               'x.ts.catch' doesn't need to have a column with Dates or any other thing, ONLY the time series (with some missing values)
# 'cc'        : matrix with the coefficient of correlation among all the time series in 'x.ts.catch'.
#               This value can be computed within this procedure, but it is a waste of time, 
#               because it is a unique value for all the iterations, so it only needs to be computed ONE time
# 'i'         : counter corresponding to the day that it is being interpolated, so, it corresponds to
#               the position of the row of 'x.ts.catch' that is is being used for the computation
# 's'         : counter corresponding to the station in which the interpolation is being computed, so, it corresponds to
#               the position of the column of 'x.ts.catch' that is receiving the computations
# 'method'    : string with the name of the method that will be used for the interpolations. Possible values are:
#               "cc-normal": normal coefficient correlation method, where all the stations 
#                            with values are used for computing the interoplated value in the target station
#               "cc-neighs": modified coefficient correlation method, where only a number of stations (provided by the user)
#                            with the highest coefficient of correlation with the target station are used for computing the interoplated value in the target station
# 'n.neighs'  : number of neighbors, with valid values (NON-'NA'), that
#               will be considered on the computation of the interpolated values.
#               ONLY required when 'mehtod'= "cc-neighs"

# Reference:
# Teegavarapu R.S.V., Chandramouli V. 2005. Improved weighting methods,
# deterministic and stochastic data-driven models for estimation of missing precipitation records.
# Journal of Hydrology, 312 (1-4), pp. 191-206.
.interpol1 <- function(x.ts.catch, cc, i, s, method="cc-neighs", n.neighs) {

    # Checking that the user provides a valid values for 'method'
    if (is.na(match( method, c("cc-normal", "cc-neighs") ) ) ) 
      stop("Invalid argument: 'method' must be of in c('cc-normal', 'cc-neighs')")
    
    # Checking that the user provides the number of neghbouring stations 
    # when he selected the "cc-neighs" method
    if ( ( method == "cc-neighs") & missing(n.neighs) ) 
        stop("Missing argument: 'n.neighs' is required when 'method' is 'cc-neighs'") 

    # index of all the stations that have a valid measure value for day = 'tdate' (NON- 'NA')
    index.valid.values <- which(!is.na( x.ts.catch[i, ] ) )
    
    # if the current station has a valid value, it is substracted from 'index.valid.values'
    # for avoiding that the current station be used for computing its own value  
    # if the current station has NOT a valid value, 'index.valid.values' remains unchanged
    index.valid.values <- setdiff(index.valid.values, s) 
    
    # number of stations with valid values
    no.na <- length( index.valid.values ) 
    
    if (method == "cc-neighs") {
    
        # All the stations with valid values ( NON-'NA') are ordered in decreasing 
        # order according to their correlation with the current station
        index.similar.st <-  sort(cc[s, index.valid.values], decreasing=TRUE )
        
        # column index in 'x.ts.catch' of the stations ordered in decreasing 
        # order according their linear correlation with the current station
        index.cc.ordered <- pmatch( names(index.similar.st), names(cc[s,]) )
        
        # auxiliar variable for computing the real minimum number of neighbors,
        # because 'length(index.cc.ordered)' is the number of stations with valid values (NON- 'NA'),
        # which can be samller than 'n.neighs' (the desired number of neghbors that have to be considered 
        # for the computation of the interpolated value) 
        n.neighs.tmp <- min(length(index.cc.ordered), n.neighs)
        
        # Index in 'x.ts.catch' of the 'n.neighs.tmp' stations with valid values (NON-'NA')
        # and better correlation with the current station
        index.neighs <-  index.cc.ordered[1:n.neighs.tmp] 
        
    } else if (method == "cc-normal") { index.neighs <- index.valid.values   }
    
    # Computation of the interpolated value
    x.idw <- sum( x.ts.catch[i, index.neighs] * cc[s, index.neighs] ) / sum( cc[s, index.neighs] )
    
    return( round(x.idw, 3) )
    
} # 'interpol1' END


#########################################################################
# 'rm1stchar' : deletes the first n characther(s) of 'x'                #
# Last updated: 01-Oct-2010                                             #
#########################################################################

# 'x'   : Character, e.g, each element may represent the name of a single gauging station
# 'n'   : numeric, indicating the number of characters that have to be removed from the beginning of 'x'
rm1stchar <- function(x, n=1) {

   if (n<0) stop("'n' must be a positive integer")

  L <- length(x)

  start.col <- n + 1
  
  x[start.col:length(x)] <- sapply(start.col:length(x), function(j,x) {
                      x[j] <- substr(x[j], start=start.col, stop=nchar(x[j]))
                  }, x = x)
  
  return(x)
  
} # 'rm1stchar' END



########################################################################
#  'sname2plot': Given a data.frame whose columns contains the ts      #
#             (without missing dates) of several gauging stations, it  #
#             takes the name of one gauging station and plots 9 graphs #
#             (see 'hydroplot' description)                            #
#                             17-Dic-2008                              #
########################################################################                          

# 'x':               : data.frame whose columns contains the time series 
#                      (without missing values) for several gauging stations.
# 'sname'            : character with the name of the station whose values will be ploted. 
#                      This name MUST eixst as column name in 'x'
#                      Additinonal columns are allowed, e.g. one column with Dates an other with IDs
# 'dates'            : "numeric", "factor", "Date" indicating how to obtain the 
#                      dates for correponding to the 'st.name' station
#                      If 'dates' is a number, it indicates the index of the column in 'x' that stores the dates
#                      If 'dates' is a factor, it have to be converted into 'Date' class, 
#                          using the date format  specified by 'date.fmt'
#                      If 'dates' is already of Date class, the following line verifies that
#                          the number of days in 'dates' be equal to the number of element in the 
#                          time series corresponding to the 'st.name' station
# 'date.fmt'         : format in which the dates are stored in 'dates'.
#                      ONLY required when class(dates)== "character", "factor" or "numeric"
# 'var.type'         : character representing the type of variable being plotted 
#                      Used for determining the function used for computing the 
#                      Monthly and Annual values when 'FUN' is missing
#                      Valid values are:
#                      -) "Precipitation" => FUN = sum
#                      -) "Temperature"   => FUN = mean
#                      -) "Flow"          => FUN = mean 
# 'FUN'             : ONLY required when 'var.type' is missing
#                     Function that have to be applied for transforming from daily to monthly or annual time step
#                     For precipitation FUN MUST be "sum"
#                     For temperature and flow time series, FUN MUST be "mean"#             
# 'na.rm'           : Logical. Should missing values be removed?
#                     TRUE : the monthly and annual values  are computed considering only those values different from NA
#                     FALSE: if there is AT LEAST one NA within a year, the monthly and annual values are NA  
# 'var.unit'		 : string repreenting the measurement unit of the variable being plotted ("mm" for precipitation, "C" for temperature, and "m3/s" for flow) 
# 'main'             : string repreenting the main title of the plot. If the user did not provide a title, this is
#                      created automatically as: main= paste(var.type, "at", st.name, sep=" "), 
# 'win.len1'		 : number of days for being used in the computation of the first moving average
# 'win.len2'		 : number of days for being used in the computation of the second moving average
# 'ptype'            : type of plot that will be plotted
#                    : ptype= "ts" => only time series
#                    : ptype= "ts+boxplot" => only time series + boxplot               
#                    : ptype= "ts+histogram" => only time series + histogram 
#                    : ptype= "ts+boxplot+histogram" => time series + boxplot + histogram    
# 'tick.tstep'       : string indicating the time step that have to be used for 
#                      putting the ticks ont he time axis. 
#                      Possible values are: 'days', 'months', 'years' 
# 'lab.tstep'        : string indicating the time step that have to be used for 
#                      putting the labels ont he time axis. 
# 'pfreq'            : Passed to the 'hydroplot' function.
#                      Character indicating how many plots are desired by the user. 
#                      Valid values are:                     
#                      -) 'dma': Daily, Monthly and Annual values are plotted 
#                      -) 'ma' : Monthly and Annual values are plotted 
#                      -) 'dm' : Daily and Monthly values are plotted 
sname2plot <- function(x, sname, var.type = "", 
                   FUN, na.rm = TRUE, var.unit = "mm", main, 
                   win.len1 = 365 * 1, win.len2 = 365 * 3, 
                   ptype = "ts+boxplot+hist", pfreq = "dma", 
                   tick.tstep= "months", lab.tstep= "years",
                   dates, date.fmt = "%Y-%m-%d") { 
                   
  # Checking the user provides 'sname'
  if (missing(sname)) { stop("Missing argument: 'sname' must be provided")  
  } else
    # Checking the the station provided for the user exist within 'x'
    if ( !(sname %in% colnames(x) ) )
      stop(paste("Invalid argument: ' The station '", sname, "' is not a column name in 'x'", sep="") )
      
  # If monthly or annual values are required, 'FUN' or 'var.type' must be provided
  if (pfreq %in% c("dma", "ma", "dm")) {     
    # Checking that the user provied a valid argument for 'var.type'       
    if (missing(FUN)) {
       # If the user did not provide a title for the plots, this is created automatically
       if (missing(var.type)) { 
          stop("Missing argument: 'var.type' OR 'FUN' must be provided")  
       } else # If 'var.type' is provided
            # Checking that the user provied a valid argument for 'var.type'       
            if (is.na(match(var.type, c("Precipitation", "Temperature", "Flow") ) ) ) {
                  stop("Invalid argument: 'var.type' must be in c('Precipitation', 'Temperature', 'Flow')") 
              } else {
                     if (var.type=="Precipitation") { 
                         FUN <- sum 
                         if (missing(var.unit)) { var.unit <- "mm"   }
                     } else if (var.type=="Temperature") { 
                              FUN <- mean
                              if (missing(var.unit)) { var.unit <- "dC" }
                            } else if (var.type=="Flow") { 
                                   FUN <- mean
                                   if (missing(var.unit)) { var.unit= "m3/s" }
                              }
                     } #ELSE end 
    } # IF end
  } # IF end
   
  # Checking the user provides the dates
  if (missing(dates)) { stop("Missing argument: 'dates' must be provided")  
  } else            
      if (is.na(match(class(dates), c("numeric", "factor", "Date")))) 
          stop("Invalid argument: 'dates' must be of class 'numeric', 'factor', 'Date'")
      
  if (is.na(match(ptype, c("ts", "ts+boxplot", "ts+hist", "ts+boxplot+hist")))) 
        stop("'ptype' valid values are: 'ts', 'ts+boxplot', 'ts+hist', 'ts+boxplot+hist'")
        
  # If 'dates' is a number, it indicates the index of the column of 'x' that stores the dates
  if ( class(dates) == "numeric" ) dates <- as.Date(x[, dates], format= date.fmt)
  
  # If 'dates' is a factor, it have to be converted into 'Date' class, 
  # using the date format  specified by 'date.fmt'
  if ( class(dates) == "factor" ) dates <- as.Date(dates, format= date.fmt)
  
  # If 'dates' is already of Date class, the following line verifies that
  # the number of days in 'dates' be equal to the number of element in the 
  # time series corresponding to the 'sname' station
  if ( ( class(dates) == "Date") & (length(dates) != nrow(x) ) ) 
     stop("Invalid argument: 'length(dates)' must be equal to 'nrow(x)'") 
     
  # If the user did not provide a title for the plots, this is created automatically
  if (missing(main)) {
        main= paste(var.type, "at station", sname, sep=" ") }
        
  
  # column index of the station identified by 'sname' within 'x'
  col.index <- which( colnames(x) == sname )
  
  # If the station name exists within 'x'
  if ( length(col.index) > 0 ) {

    # Slecting the time series within 'x' corresponding to the 'sname' station
    x <- x[ ,col.index]
    
    # Transform the vector of time series ('x') and the vector with dates ('dates')
    # into a zoo variable, using the format psecified by 'date.fmt'
    x <- vector2zoo(x, dates, date.fmt)
    
    # If the user didn't provided a value for 'tick.tstep' and 
    # the lenght of the daily ts is less than 1 year, the ticks in 
    # the 'x' axis are placed by day; if larger than a year, they are placed by month
    if ( missing(tick.tstep) ) {
      if (length(x) <= 366) {
        tick.tstep <- "days"
      } else tick.tstep <- "months"
    } # IF end
        
    # If the user didn't provided a value for 'lab.tstep' and 
    # the lenght of the daily ts is less than 1 year, the labels in 
    # the 'x' axis are placed by month; if larger than a year, they are placed by year
    if ( missing(lab.tstep) ) {
      if (length(x) <=366) {
        lab.tstep <- "months"
      } else lab.tstep <- "years"
    } # IF end
     
    # 9 plots:
    # 1: Line plot with Daily time series, with 2 moving averages, specified by 'win.len1' and 'win.len2'
    # 2: Line plot with Monthly time series, with 2 moving averages, specified by 'win.len1' and 'win.len2'
    # 3: Line plot with Annual time series
    # 4: Boxplot with daily time series
    # 5: Boxplot with monthly time series
    # 6: Boxplot with annual time series
    # 7: Histogram of the daily time series
    # 8: Histogram of the monthly time series
    # 9: Histogram of the annual time series
    hydroplot(x, sname= sname, var.type=var.type, 
              FUN=FUN, na.rm=na.rm, var.unit=var.unit, 
              main=main, win.len1=win.len1, win.len2=win.len2, 
              ptype=ptype, pfreq = pfreq, tick.tstep= tick.tstep, lab.tstep= lab.tstep)
                 
  } else stop( paste("The station name", sname, "does not exist in 'x'", sep=" ") )
 
}  # 'sname2plot' END


#####################################################
# sname2ts: Station name -> time series             #
#####################################################
#                 January 13th, 2009                #
#####################################################
# This function takes a data.frame whose columns contains the time series 
# (without missing dates) of several gauging stations, it takes the name 
# of one gauging station and extracts a time 
# series with daily, monthly or annual time step
# 'x'        : data.frame containing the complete (without missing dates) 
#              times series of all the stations.
#              It can also contain 1 column with the dates of the measurements, 
#              or they can be provided in a separated way
# 'sname'    : string representing the name of the station, which have to correspond 
#              with one column name in 'x'
# 'tstep.out': character that defines the time step of the desired output time series
#              it must be one of { "daily", "monthly", "annual" }
# 'dates'    : "numeric", "factor", "Date" indicating how to obtain the 
#              dates for correponding to the 'sname' station
#              If 'dates' is a number, it indicates the index of the column in 
#                'x' that stores the dates
#              If 'dates' is a factor, it have to be converted into 'Date' class, 
#                using the date format  specified by 'date.fmt'
#              If 'dates' is already of Date class, the following line verifies that
#                the number of days in 'dates' be equal to the number of element in the 
#                time series corresponding to the 'st.name' station
# 'date.fmt' : format in which the dates are stored in 'dates'.
#              ONLY required when class(dates)=="factor" or "numeric"
# 'var.type' : character representing the type of variable being plotted 
#              Used for determining the function used for computing the 
#              Monthly and Annual values when 'FUN' is missing
#              Valid values are:
#              -) "Precipitation" => FUN = sum
#              -) "Temperature"   => FUN = mean
#              -) "Flow"          => FUN = mean 
# 'FUN'      : ONLY required when 'var.type' is missing
#              Function that have to be applied for transforming from daily to monthly or annual time step
#              For precipitation FUN MUST be "sum"
#              For temperature and flow time series, FUN MUST be "mean"#             
# 'na.rm'    : Logical. Should missing values be removed?
#              TRUE : the monthly and annual values  are computed considering only those values different from NA
#              FALSE: if there is AT LEAST one NA within a year, the monthly and annual values are NA 
sname2ts <- function(x, sname, dates, date.fmt="%Y-%m-%d", var.type, 
                     tstep.out="daily", FUN, na.rm=TRUE) {

  # Checking that the user provied a valid argument for 'x'
  if ( is.na( match( class(x), c("data.frame") ) ) ) 
      stop("Invalid argument: 'x' must be of class 'data.frame'")
      
  # Checking the user provides 'sname'
  if (missing(sname)) { stop("Missing argument: 'sname' must be provided")  
  } else
    # Checking the the station provided for the user exist within 'x'
    if ( !(sname %in% colnames(x) ) )
      stop(paste("Invalid argument: ' The station '", sname, "' is not a column name in 'x'", sep="") )
      
  # If monthly or annual values are required, 'FUN' or 'var.type' must be provided
  if (tstep.out != "daily") {
    # Checking that the user provied a valid argument for 'var.type'       
    if (missing(FUN)) {
       # If the user did not provide a title for the plots, this is created automatically
       if (missing(var.type)) { 
         stop("Missing argument: 'var.type' OR 'FUN' must be provided")  
       } else # If 'var.type' is provided
        # Checking that the user provied a valid argument for 'var.type'       
        if (is.na(match(var.type, c("Precipitation", "Temperature", "Flow") ) ) ) {
              stop("Invalid argument: 'var.type' must be in c('Precipitation', 'Temperature', 'Flow')") 
          } else {
                 if (var.type=="Precipitation") { FUN= sum  } 
                 else if (var.type=="Temperature") { FUN= mean  } 
                 else if (var.type=="Flow") { FUN= mean  }
                 } #ELSE end 
    } 
  } # IF end
  
  # Checking that the user provied a valid argument for 'tstep.out'  
  if (is.na(match( tstep.out, c("daily", "monthly", "annual") ) ) ) 
      stop("Invalid argument: 'tstep.out' must be in c('daily', 'monthly', 'annual'")
      
  # Checking the user provides the dates
  if (missing(dates)) { stop("Missing argument: 'dates' must be provided")  
  } else
    # Checking that the user provided a valid argument for 'dates'  
    if (is.na(match(class(dates), c("numeric", "factor", "Date")))) 
        stop("Invalid argument: 'dates' must be of class 'numeric', 'factor', 'Date'")
        
  # If 'dates' is a number, it indicates the index of the column of 'x' that stores the dates
  if ( class(dates) == "numeric" ) dates <- as.Date(x[, dates], format= date.fmt)
  
  # If 'dates' is a factor, it have to be converted into 'Date' class, 
  # using the date format  specified by 'date.fmt'
  if ( class(dates) == "factor" ) dates <- as.Date(dates, format= date.fmt)
  
  # If 'dates' is already of Date class, the following line verifies that
  # the number of days in 'dates' be equal to the number of element in the 
  # time series corresponding to the 'st.name' station
  if ( ( class(dates) == "Date") & (length(dates) != nrow(x) ) ) 
     stop("Invalid argument: 'length(dates)' must be equal to 'nrow(x)'")                    
        
  # column index of the station identified by 'sname' within 'x'
  col.index <- which( colnames(x) == sname )
      
  # If the station name exists within 'x'
  if ( length(col.index) > 0 ) {

  # Selecting the time series within 'x' corresponding to the 'sname' station
  x <- x[ ,col.index]
    
  # Transform the vector of time series ('x') and the vector with dates ('dates')
  # into a zoo variable, using the format psecified by 'date.fmt'
  x.daily   <- vector2zoo(x, dates, date.fmt="%Y-%m-%d")
            
  if (tstep.out == "daily") { return (x.daily)
      
    } else if (tstep.out =="monthly") {
    
      # Transformation from daily to monthly
      x.monthly <- daily2monthly(x.daily, FUN, na.rm )
      
      return (x.monthly)
      
      } else if (tstep.out =="annual") {
      
        # Transformation from daily to annual
        x.annual  <- daily2annual(x.daily, FUN, na.rm )
        
        return (x.annual)
        
        } # IF/ELSE/ELSE END

  } # IF end

} # 'sname2ts' END



#########################################################################
#  stdx : standarizes a vector or matrix, i.e., scales all the values   #
#         in in a way that the transformed values will be within the    #
#         range [0,1]. z = scale(x) = [ x - xmin ] / [ xmax - xmin ]    #  
#########################################################################
#               February 19th, 2009                 #
#####################################################
# This function 

# 'x'     : vector or matrix to be scaled
# 'result': standarized 'x', where all the values of each column of 'x' 
#           are within the range [0,1]

# If you are (very likely) interested in Back transforming this standarized 
# values into the original ranges, you sould keep the following values:
#        xmin   <-  apply(x,2,min, na.rm=TRUE)
#        xrange <-  apply(x,2,range, na.rm=TRUE)
#        xrange <-  apply(xrange,2, diff, na.rm=TRUE) 

stdx <-function(x,...) UseMethod("stdx")
 
stdx.default <- function (x,...) {
  
     if (is.na(match(class(x), c("ts", "zoo") ) ) ) x <- as.numeric(x)

     # range of 'x', i.e., r = xmax - xmin
     r <- diff( range(x, na.rm=TRUE) )
     
     if ( r ==0 ) { std <- x
     
     } else { std <- scale(x, center=min(x, na.rm=TRUE), scale= r) }
     
     names(std) <- names(x)
     
     return( as.numeric(std) )
     
  } # 'stdx.default' end
  
  
stdx.matrix <- function(x,...) { 
          
  std <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
  
  std <- sapply(1:ncol(x), function(i) { 
		 std[,i] <- stdx.default( x[,i] )
		 })
   
  colnames(std) <- colnames(x)
  rownames(std) <- rownames(x)                             
     
  return(std)  
 
}  # 'stdx.matrix' end  


stdx.data.frame <- function(x,...) { 
          
  x <- as.matrix(x)
  NextMethod("stdx.matrix")
 
}  # 'stdx.data.frame' end   




#########################################################################
#  istdx : This function transforms back a standarized vector/matrix    #
#           into their original values, i.e., re-scales all the values  #
#           in the [0,1] interval to the original range of values       #
#           x = re-scale(z) = z*[ zmax - zmin ] + xmin                  #
#########################################################################
#               February 19th, 2009                 #
#####################################################
# 

# 'z'     : standarized vector or matriz to be re-scaled, all the values
#           have to be in the range [0,1]
# 'xmin'  : numeric with the minimum value(s) in the original 'x'
#           -) if 'z' is a vector, 'xmin' has to be a real
#           -) if 'z' is a matrix/data.frame, 'xmin' has to be a vector, with
#              the minimum values for each column of the original 'x' 
#              in this case, the vector of minimums can be otained as: 
#               xmin <-  apply(x,2,min, na.rm=TRUE)
# 'xrange'  : numeric with the range of value(s) in the original 'x'
#           -) if 'z' is a vector, 'xrange' has to be a real
#           -) if 'z' is a matrix/data.frame, 'xrange' has to be a vector, with
#              the range of values for each column of the original 'x' 
#              in this case, the vector of ranges can be otained as: 
#               xrange <-  apply(x, 2, range, na.rm=TRUE)
#               xrange <-  apply(xrange,2, diff, na.rm=TRUE) 
# 'result': re-scaled 'x', where all the values of each column of 'x' 
#           are within the original range of x values

istdx <-function(x,...) UseMethod("istdx")
 
istdx.default <- function(x, xmin, xrange,...) { 

   z <- x*xrange + xmin
     
   names(z) <- names(x)
     
   return(as.numeric(z))
     
} # 'istdx.default' end


istdx.matrix <- function(x, xmin, xrange,...) {  
          
  z <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
  
  z <- sapply(1:ncol(z), function(i) { 
		 z[,i] <- istdx.default( x[,i], xmin[i], xrange[i] )
		 })
   
  colnames(z) <- colnames(x)
  rownames(z) <- rownames(x) 
     
  return(z)  
 
}  # 'istdx.matrix' end  


istdx.data.frame <- function(x, xmin, xrange,...) {  
          
  x <- as.matrix(x)
  NextMethod("istdx.matrix")
 
}  # 'istdx.data.frame' end  



#############################################################################
# time2season : This function transforms a vector of dates into a vector of #
#               seasons (summer, winter, autumm, spring), considering that: #
#               winter = DJF: December, January, February                   #
#               spring = MAM: March, April, May                             # 
#               summer = JJA: June, July, August                            #
#               autumm = SON: September, October, November                  #
#############################################################################
#                       March 18th, 2009                                    #
#############################################################################

# 'x'       : vector with the dates that have to be transformed. class(x) must be "Date"
# 'out.fmt' : format of the output seasons. Possible values are:
#             -) 'seasons' =>  "winter", "spring",  "summer", autumm"
#             -) 'months'  =>  "DJF", "MAM",  "JJA", SON" 

# 'result': vector with the wheater season to which each date in 'x' belongs

time2season <- function(x, out.fmt="months") {

 # Checking that 'class(x)==Date'
  if (is.na(match(class(x), c("Date") ) ) ) 
     stop("Invalid argument: 'x' must be of class 'Date'")
     
 # Checking the class of out.fmt
  if (is.na(match(out.fmt, c("seasons", "months") ) ) ) 
     stop("Invalid argument: 'out.fmt' must be in c('seasons', 'months')")
     
 months <- format(x, "%m")
 
 winter <- which( months %in% c("12", "01", "02") )
 spring <- which( months %in% c("03", "04", "05") )
 summer <- which( months %in% c("06", "07", "08") )
 autumm <- which( months %in% c("09", "10", "11") )
 
 # Creation of the output, with the same length of the 'x' input
 seasons <- rep(NA, length(x))
 
 if (out.fmt == "seasons") {
    seasons[winter] <- "winter"
    seasons[spring] <- "spring"
    seasons[summer] <- "summer"
    seasons[autumm] <- "autumm"
 } else { # out.fmt == "months"
    seasons[winter] <- "DJF"
    seasons[spring] <- "MAM"
    seasons[summer] <- "JJA"
    seasons[autumm] <- "SON"
 } # IF end
 
 return(seasons)

} # 'time2season' END


#########################################################################
# extractzoo: Extracts from a zoo object, all the values belonging to a #
#             given month, year or weather season                       #
#########################################################################
#           April 16th, 2009;  May 15th, 2009; Ago 30th, 2009           #
#########################################################################

# 'x'    : variable of type 'zoo'
# 'trgt' : numeric or character indicating elements to extract from 'x'
#          Valid values are:
#          1) integer from 1 to 12: 'trgt' is considered as a month, and
#             all the vaues in 'x' belonging to the month specified by 'trgt' 
#             will be extracted  (1=JAN, 2=FEB,...., 12=DEC)
#          2) integer > 12: 'trgt' is considered as a year, and all the
#             values in 'x' belonging to the year specified by 'trgt' 
#             will be extracted
#          3) character: 'trgt' is considered as a weather season, and 
#             all the values in 'x' belonging to the season specified by 
#             'trgt' will be extracted. Valid values are:
#              -) "DJF": December, January, February                    
#              -) "MAM": March, April, May                              
#              -) "JJA": June, July, August                            
#              -) "SON": September, October, November    

extractzoo <- function(x, trgt, ...) {

  # Checking that the user provied a valid argument for 'trgt'
  if ( is.na( match( class(trgt), c("integer", "numeric", "character") ) ) ) 
      stop("Invalid argument: class('trgt') must be in c('integer', 'numeric','character')")	  
	  
  # If 'trgt' is a month or a year
  if ( (class(trgt)=="numeric") | (class(trgt)=="integer")) {
  
    # Checking that 'trgt' is integer
    if ( trgt - trunc(trgt) != 0 ) 
        stop("Invalid argument: 'trgt' must be integer") 
		
	if (trgt %in% 1:12) {
	   index <- which( as.numeric(format(time(x), "%m")) == trgt )          
	} else index <- which( as.numeric(format(time(x), "%Y")) == trgt )          


  } # if END
  
    # if 'trgt' is a weather season 
    else if (class(trgt)=="character") {
    
	  # Checking a valid value for 'trgt'
      if (is.na(match(trgt, c("DJF", "MAM", "JJA", "SON") ) ) ) 
         stop("Invalid argument: 'trgt' must be in c('DJF', 'MAM', 'JJA', 'SON')")      
	  
	  # Gets the season each element of 'x' belongs to
      seasons <- time2season(time(x), out.fmt="months") 
  
      # Selects only those elements of 'x' belonging to the desired season
      index <- which(seasons == trgt)        
	} # ELSE end  
	  
  return(x[index])
  
} # 'extractzoo' END




#########################################################################
# dm2seasonal: Generic function for computing seasonal values for every #
#              year of a daily/monthly zoo object                       #
#########################################################################
#                  May 15th, 2009; Aug 31th 2009                        #
#########################################################################

# 'x   '    : variable of type 'zoo' or 'data.frame'
# 'season'  : character, indicating the weather season to be used for selecting the data
#             Valid values are:
#             -) "DJF": December, January, February                    
#             -) "MAM": March, April, May                              
#             -) "JJA": June, July, August                            
#             -) "SON": September, October, November  
# 'FUN'      : Function that will be applied to ALL the values of 'x' belonging to the given weather season          
# 'na.rm'    : Logical. Should missing values be removed?
#              TRUE : the seasonal values  are computed considering only those values different from NA
#              FALSE: if there is AT LEAST one NA within a season, the corresponding values are NA      
dm2seasonal <-function(x, ...) UseMethod("dm2seasonal")


dm2seasonal.default <- function(x, season, FUN, na.rm=TRUE, ...) {
      
  # Checking that the user provied a valid argument for 'x'
  if ( is.na( match( class(x), c("zoo") ) ) ) 
      stop("Invalid argument: 'class(x)' must be 'zoo'")
      
  # Checking that the user provied a valid argument for 'season'
  if ( missing(season) ) {
      stop("Missing argument: 'season' must be provided")
  } else # If 'season' is provided
      # Checking a valid value for 'season'
      if (is.na(match(season, c("DJF", "MAM", "JJA", "SON") ) ) ) 
         stop("Invalid argument: 'season' must be in c('DJF', 'MAM', 'JJA', 'SON')")    
  
  # Checking that the user provied a valid argument for 'FUN'
  if (missing(FUN)) stop("Missing argument: 'FUN' must be provided")
			   
  # Checking that 'x' is a Daily or Monthly ts
  if (is.na(match(sfreq(x), c("daily", "monthly") ) ) )
      stop(paste("Invalid argument: 'x' must be a daily or monthly ts, but 'sfreq(x)' is", sfreq(x), sep=" ")  )
  
  dates <- time(x)
  
  # Computing the Starting and Ending Year of the analysis
  Starting.Year <- as.numeric(format(range(dates)[1], "%Y"))
  Ending.Year   <- as.numeric(format(range(dates)[2], "%Y"))
  
  # Amount of Years belonging to the desired season
  nyears <- Ending.Year - Starting.Year + 1
  
  # Requiring the Zoo Library
  require(zoo)
  
  # If 'x' is a daily ts, the Monthly values are computed
  #if ( sfreq(x) == "Daily")  {
  #   x <- daily2monthly(x, FUN, na.rm ) }
  
  # Getting the Monthly values beloonging ONLY to the desired weather season
  s <- extractzoo(x= x, trgt=season) 
  
  # Moving forward all the December values, in order that 
  # December of 1991 be used together with Jan/92 and Feb/92, 
  # instead of with Jan/91 and Feb/91
  if (season == "DJF") {
			syears            <- as.numeric(format( time(s), "%Y" ))
			dec.index         <- which(format(time(s), "%m") == 12)
			dec.years         <- syears[dec.index]
			dec.years         <- dec.years + 1
			syears[dec.index] <- dec.years               
			
			s.a <- aggregate( s, by= syears, FUN=FUN, na.rm= na.rm) 
			
			# Removing the last value of december, because it is outside of the analysis period
			s.a <- s.a[1:(length(s.a)-1)]
  } else  {
  
			s.a <- aggregate( s, by= format( time(s), "%Y" ), FUN=FUN, na.rm= na.rm )  
	}  
		  
  # Getting the position of all the years in which there were no values
  # mean(NA:NA, na.rm=TRUE) == NaN
  nan.index <- which(is.nan(s.a))
  
  # Changing all the NaN's by NA's
  if ( length(nan.index) > 0 ) { s.a[nan.index] <- NA }
  
  # Getting the position of all the years in which there were no values
  # min(NA:NA, na.rm=TRUE) == Inf  ; max(NA:NA, na.rm=TRUE) == -Inf
  inf.index <- which(is.infinite(s.a))
  
  # Changing all the Inf and -Inf by NA's
  if ( length(inf.index) > 0 ) { s.a[inf.index] <- NA }
 
  return(s.a)
  
} # 'dm2seasonal.default' END



# 'dates'   : "numeric", "factor", "Date" indicating how to obtain the 
#             dates for correponding to the 'sname' station
#             If 'dates' is a number, it indicates the index of the column in 
#                'x' that stores the dates
#             If 'dates' is a factor, it have to be converted into 'Date' class, 
#                using the date format  specified by 'date.fmt'
#             If 'dates' is already of Date class, the following line verifies that
#                the number of days in 'dates' be equal to the number of element in the 
#                time series corresponding to the 'st.name' station
# 'date.fmt': character indicating the format in which the dates are stored in 'dates'.
#             ONLY required when class(dates)=="factor" or "numeric"
# 'out.type': string that define the desired type of output. Possible values are
#             -) "data.frame": a data.frame, with as many columns as stations 
#                              are included in 'x'
#             -) "db"        : a data.frame, with 4 colums will be produced.
#                              The first column will store the ID of the station, 
#                              The second column will store the Year, 
#                              The third column will store the season
#                              The fouth column will contain the seasonal value, corresponding to the year specified in the second column
dm2seasonal.data.frame <- function(x, season, FUN, na.rm=TRUE,
                                   dates, date.fmt="%Y-%m-%d", 
								   out.type="data.frame",... ) {
  
  # Checking that the user provied a valid argument for 'out.type'  
  if (is.na(match( out.type, c("data.frame", "db") ) ) ) 
      stop("Invalid argument: 'out.type' must be in c('data.frame', 'db'")
      
  # Checking that the user provied a valid argument for 'season'
  if ( missing(season) ) {
      stop("Missing argument: 'season' must be provided")
  } else # If 'season' is provided
      # Checking a valid value for 'season'
      if (is.na(match(season, c("DJF", "MAM", "JJA", "SON") ) ) ) 
         stop("Invalid argument: 'season' must be in c('DJF', 'MAM', 'JJA', 'SON')") 
      
  # Checking that the user provied a valid argument for 'FUN'
  if (missing(FUN)) stop("Missing argument: 'FUN' must be provided")
	  
  # Checking that the user provied a valid argument for 'dates'       
  if (missing(dates)) {
      stop("Missing argument: 'dates' must be provided") 
  } else      
     # Checking that the user provied a valid argument for 'dates'  
     if (is.na(match(class(dates), c("numeric", "factor", "Date")))) 
         stop("Invalid argument: 'dates' must be of class 'numeric', 'factor', 'Date'")
        
  # If 'dates' is a number, it indicates the index of the column of 'x' that stores the dates
  # The column with dates is then substracted form 'x' for easening the further computations
  if ( class(dates) == "numeric" ) {    
    tmp   <- dates
    dates <- as.Date(x[, dates], format= date.fmt)
    x     <- x[-tmp]
  }  # IF end 
  
  # If 'dates' is a factor, it have to be converted into 'Date' class, 
  # using the date format  specified by 'date.fmt'
  if ( class(dates) == "factor" ) dates <- as.Date(dates, format= date.fmt)
  
  # If 'dates' is already of Date class, the following line verifies that
  # the number of days in 'dates' be equal to the number of element in the 
  # time series corresponding to the 'st.name' station
  if ( ( class(dates) == "Date") & (length(dates) != nrow(x) ) ) 
     stop("Invalid argument: 'length(dates)' must be equal to 'nrow(x)'")
     
  # Checking a valid value for 'season'
  if (is.na(match(season, c("DJF", "MAM", "JJA", "SON") ) ) ) 
     stop("Invalid argument: 'season' must be in c('DJF', 'MAM', 'JJA', 'SON')")  
     
  # Amount of stations in 'x'
  nstations <- ncol(x)

  # ID of all the stations in 'x'
  snames <- colnames(x)  
  
  # Computing the Starting and Ending Year of the analysis
  Starting.Year <- as.numeric(format(range(dates)[1], "%Y"))
  Ending.Year   <- as.numeric(format(range(dates)[2], "%Y"))
  
  # Amount of Years belonging to the desired season
  nyears <- Ending.Year - Starting.Year + 1
  
  # Requiring the Zoo Library (Z’s ordered observations)
  require(zoo)
  
  print("Starting the computations...", quote=FALSE )
  
  if (out.type == "data.frame") {
  
	# Vector with the names of the field that will be used for storing the results
	field.names <- c(snames )  

	# Creating the data.frame that will store the computed averages for each station
	z <- as.data.frame(matrix(data = NA, nrow = nyears, ncol = nstations, 
						byrow = TRUE, dimnames = NULL) )
	colnames(z) <- field.names
	
	rownames(z) <- Starting.Year:Ending.Year         
	  
	z[1:nstations] <- sapply(1:nstations, function(j,y) {

		print( paste("Station: ", format(snames[j], width=10, justify="left"),
					 " : ",format(j, width=3, justify="left"), "/", 
					 nstations, " => ", 
					 format(round(100*j/nstations,2), width=6, justify="left"), 
					 "%", sep=""), quote=FALSE )
						
		# Transforming the column of 'x' into a zoo object, 
		# using the dates provided by the user
		tmp <- vector2zoo(x=x[,j], dates=dates, date.fmt=date.fmt)
		
		z[,j] <- dm2seasonal.default(x= tmp, season=season, FUN=FUN, na.rm=na.rm) 
						 
	}, y = x) # sapply END
  
  } else if (out.type == "db") {
  
        # Creating a vector with the names of the field that will be used for storing the results
        field.names <- c("StationID", "Year", "Season", "Value" )  

        # Creating the data.frame that will store the computed averages for each subcatchment
        z <- as.data.frame(matrix(data = NA, nrow = nyears*nstations, ncol = 4, 
                            byrow = TRUE, dimnames = NULL) )
        colnames(z) <- field.names      
        
        y = x
        
        for (j in 1:nstations) {
            
            print( paste("Station: ", format(snames[j], width=10, justify="left"),
                         " : ",format(j, width=3, justify="left"), "/", 
                         nstations, " => ", 
                         format(round(100*j/nstations,2), width=6, justify="left"), 
                         "%", sep=""), quote=FALSE )
                            
            # Transforming the column of 'x' into a zoo object, 
		    # using the dates provided by the user
		    tmp <- vector2zoo(x=x[,j], dates=dates, date.fmt=date.fmt)
		
		    # Computing the seasonal values
		    s.a <- dm2seasonal.default(x= tmp, season=season, FUN=FUN, na.rm=na.rm)
     
            # Putting the annual seasonal values in the output data.frame
            # The first column of 'x' corresponds to the Year
            row.ini <- (j-1)*nyears + 1
            row.fin <- j*nyears
            
            z[row.ini:row.fin, 1] <- snames[j] # it is automatically repeted 'nyears' times
            z[row.ini:row.fin, 2] <- Starting.Year:Ending.Year            
            z[row.ini:row.fin, 3] <- season    # it is automatically repeted 'nyears' times
            z[row.ini:row.fin, 4] <- s.a            
                              
        } # FOR end
  
    } # IF end
   
  return( z )
  
 } #'dm2seasonal.data.frame' END
 
 
 
########################################################################
# monthlyfunction: Generic function for applying any R function to     #
#                  ALL the values in 'x' belonging to a given month    #
########################################################################
#                  May 15th, 2009; Aug 31th 2009                       #
########################################################################
# 'x   '    : variable of type 'zoo' or 'data.frame', with daily or monthly frequency
# 'FUN'      : Function that will be applied to ALL the values in 'x' belonging to each one of the 12 months of the year
# 'na.rm'    : Logical. Should missing values be removed?
#              TRUE : the monthly values  are computed considering only those values in 'x' different from NA
#              FALSE: if there is AT LEAST one NA within a month, the FUN and monthly values are NA    
monthlyfunction <- function(x, ...) UseMethod("monthlyfunction")

monthlyfunction.default <- function(x, FUN, na.rm=TRUE,...) { 

	 # Checking that 'x' is a zoo object
	 if (is.na(match(class(x), c("zoo")))) 
			stop("Invalid argument: 'x' must be of class 'zoo'") 
            
     # Checking that the user provied a valid argument for 'FUN'       
     if (missing(FUN))
         stop("Missing argument: 'FUN' must be provided")
	 
	 # Checking the user provide a valid value for 'x'
	 if (is.na(match(sfreq(x), c("daily", "monthly")))) {
		 stop(paste("Invalid argument: 'x' is not a daily or mothly ts, it is a ", sfreq(x), " ts", sep="") ) }			 
         
     # Requiring the Zoo Library (Z’s ordered observations)
	 require(zoo) 	 
	 
	 # 'as.numeric' is necessary for being able to change the names to the output
	 totals <- as.numeric( aggregate( x, by= format( time(x), "%m" ), FUN=FUN, na.rm= na.rm ) )
     
     # Replacing the NaNs by 'NA.
     # NaN's are obtained when using theFUN=mean with complete NA values
     nan.index          <- which(is.nan(totals))                              
     totals[ nan.index] <- NA
     
     # Getting the position of all the years in which there were no values
     # min(NA:NA, na.rm=TRUE) == Inf  ; max(NA:NA, na.rm=TRUE) == -Inf
     inf.index <- which(is.infinite(totals))
  
     # Changing all the Inf and -Inf by NA's
     if ( length(inf.index) > 0 ) { totals[inf.index] <- NA }
	 
	 names(totals) <- month.abb
	 
	 return(totals)

} # 'monthlyfunction.default' end



# 'dates'   : "numeric", "factor", "Date" indicating how to obtain the 
#             dates for correponding to the 'sname' station
#             If 'dates' is a number, it indicates the index of the column in 
#                'x' that stores the dates
#             If 'dates' is a factor, it have to be converted into 'Date' class, 
#                using the date format  specified by 'date.fmt'
#             If 'dates' is already of Date class, the following line verifies that
#                the number of days in 'dates' be equal to the number of element in the 
#                time series corresponding to the 'st.name' station
# 'date.fmt': format in which the dates are stored in 'dates'.
#             ONLY required when class(dates)=="factor" or "numeric"
# 'out.type': string that define the desired type of output. Possible values are
#             -) "data.frame": a data.frame, with 12 columns representing the months, 
#                              and as many rows as stations are included in 'x'
#             -) "db"        : a data.frame, with 4 colums will be produced.
#                              The first column stores the ID of the station
#                              The second column stores the Year, 
#                              The third column stores the ID of the station,
#                              The fourth column contains the monthly value corresponding to the year specified in the second column
# 'verbose'      : logical; if TRUE, progress messages are printed 
monthlyfunction.data.frame <- function(x, FUN, na.rm=TRUE,
                                    dates, date.fmt="%Y-%m-%d", 
								    out.type="data.frame", 
									verbose=TRUE,...) {
  
  # Checking that the user provied a valid argument for 'out.type'  
  if (is.na(match( out.type, c("data.frame", "db") ) ) ) 
      stop("Invalid argument: 'out.type' must be in c('data.frame', 'db'")
      
   # Checking that the user provied a valid argument for 'FUN'       
   if (missing(FUN))
         stop("Missing argument: 'FUN' must be provided")
	  
  # Checking that the user provied a valid argument for 'dates'       
  if (missing(dates)) {
      stop("Missing argument: 'dates' must be provided") 
  } else  
    {    
       # Checking that the user provied a valid argument for 'dates'  
       if (is.na(match(class(dates), c("numeric", "factor", "Date")))) 
           stop("Invalid argument: 'dates' must be of class 'numeric', 'factor', 'Date'")
           
        # Verification that the number of days in 'dates' be equal to the number 
        # of elements in 'x'
        if ( ( class(dates) == "Date") & (length(dates) != nrow(x) ) ) 
        stop("Invalid argument: 'length(dates)' must be equal to 'nrow(x)'") 
    } # ELSE end
        
  # If 'dates' is a number, it indicates the index of the column of 'x' that stores the dates
  # The column with dates is then substracted form 'x' for easening the further computations
  if ( class(dates) == "numeric" ) {    
    tmp   <- dates
    dates <- as.Date(x[, dates], format= date.fmt)
    x     <- x[-tmp]
  }  # IF end 
  
  # If 'dates' is a factor, it have to be converted into 'Date' class, 
  # using the date format  specified by 'date.fmt'
  if ( class(dates) == "factor" ) dates <- as.Date(dates, format= date.fmt)
     
  # Amount of stations in 'x'
  nstations <- ncol(x)

  # ID of all the stations in 'x'
  snames <- colnames(x)  
  
  # Computing the Starting and Ending Year of the analysis
  Starting.Year <- as.numeric(format(range(dates)[1], "%Y"))
  Ending.Year   <- as.numeric(format(range(dates)[2], "%Y"))
  
  # Amount of Years belonging to the desired period
  nyears <- Ending.Year - Starting.Year + 1
  
  # Amount of months belonging to the desired period
  nmonths <- 12*nyears
  
  # Requiring the Zoo Library (Z’s ordered observations)
  require(zoo)
  
  if (verbose) print("Starting the computations...", quote=FALSE )
  
  if (out.type == "data.frame") {
  
	# Vector with the names of the field that will be used for storing the results
	field.names <- month.abb

	# Creating the data.frame that will store the computed averages for each station
	z <- as.data.frame(matrix(data = NA, nrow = nstations, ncol = 12, 
						byrow = TRUE, dimnames = NULL) )
	  
	z <- sapply(1:nstations, function(j,y) {

		
		if (verbose) print( paste("Station: ", format(snames[j], width=10, justify="left"),
					              " : ",format(j, width=3, justify="left"), "/", 
					              nstations, " => ", 
					              format(round(100*j/nstations,2), width=6, justify="left"), 
					              "%", sep=""), quote=FALSE )
						
		# Transforming the column of 'x' into a zoo object, 
		# using the dates provided by the user
		tmp <- vector2zoo(x=x[,j], dates=dates, date.fmt=date.fmt)
		
		# Computing the annual values
		z[,1:12] <- monthlyfunction.default(x= tmp, FUN=FUN, na.rm=na.rm)
						 
	}, y = x) # sapply END
    
    z <- t(z) # I don't know WHY !!
    rownames(z) <- snames 
  
  } else if (out.type == "db") {
  
        # Creating a vector with the names of the field that will be used for storing the results
        field.names <- c("StationID", "Year", "Month", "Value" )  

        # Creating the data.frame that will store the computed averages for each subcatchment
        z <- as.data.frame(matrix(data = NA, nrow = nmonths*nstations, ncol = 4, 
                           byrow = TRUE, dimnames = NULL) )        
        
        y = x
        
        for (j in 1:nstations) {
            
            if (verbose) print( paste("Station: ", format(snames[j], width=10, justify="left"),
                                      " : ", format(j, width=3, justify="left"), "/", 
                                      nstations, " => ", 
                                      format(round(100*j/nstations,2), width=6, justify="left"), 
                                      "%", sep=""), quote=FALSE )
                            
            # Transforming the column of 'x' into a zoo object, 
		    # using the dates provided by the user
		    tmp <- vector2zoo(x=x[,j], dates=dates, date.fmt=date.fmt)
		
		    # Computing the annual values
		    tmp <- monthlyfunction.default(x= tmp, FUN=FUN, na.rm=na.rm)
			
			# Putting the annual/monthly values in the output data.frame
            # The first column of 'x' corresponds to the Year
            row.ini <- (j-1)*nmonths + 1
            row.fin <-  j*nmonths
            
            z[row.ini:row.fin, 1] <- snames[j] # it is automatically repeted 'nmonths' times
            z[row.ini:row.fin, 2] <- rep(Starting.Year:Ending.Year, each=12)
            z[row.ini:row.fin, 3] <- month.abb        
            z[row.ini:row.fin, 4] <- tmp
                              
        } # FOR end
        
        colnames(z) <- field.names 
  
    } # IF end  
   
  return( z )
  
 } #'monthlyfunction.data.frame' END
 
 
 
 #######################################################################
# seasonalfunction: Generic function for applying any R function to    #
#                   summarize the seasonal values of 'x'               #
########################################################################
#                        Sep 11th, 2009                                #
########################################################################
# 'x   '    : variable of type 'zoo' or 'data.frame'
# 'FUN'      : Function that will be applied to ALL the values in 'x' belonging to each one of the 4 weather seasons
#              (e.g., Fun can be some of c('mean', 'max', 'min', 'sd'))          
# 'na.rm'    : Logical. Should missing values be removed?
#              TRUE : the monthly values  are computed considering only those values in 'x' different from NA
#              FALSE: if there is AT LEAST one NA within a month, the FUN and monthly values are NA    
seasonalfunction <- function(x, ...) UseMethod("seasonalfunction")

seasonalfunction.default <- function(x, FUN, na.rm=TRUE,...) { 

	 # Checking that 'x' is a zoo object
	 if (is.na(match(class(x), c("zoo")))) 
			stop("Invalid argument: 'x' must be of class 'zoo'") 
        
     # Checking that the user provied a valid argument for 'FUN'       
     if (missing(FUN))
         stop("Missing argument: 'FUN' must be provided") 
	 
	 # Checking the user provide a valid value for 'x'
	 if (is.na(match(sfreq(x), c("daily", "monthly")))) {
		 stop(paste("Invalid argument: 'x' is not a daily or monthly ts, it is a ", sfreq(x), " ts", sep="") ) }			 
         
     # Requiring the Zoo Library (Z’s ordered observations)
	 require(zoo) 
     
     seasons <- c("DJF", "MAM", "JJA", "SON")
     
     # Creating the output variable
     z <- NA*numeric(4)       
     
     z <- sapply(1:4, function(j) {
		
		s <- dm2seasonal(x, season=seasons[j], FUN=FUN, na.rm=na.rm)
        
        # 'as.numeric' is necessary for being able to change the names to the output
	    z[j] <- as.numeric( aggregate( s, by= rep(seasons[j], length(s)), FUN=FUN, na.rm= na.rm ) )
						 
	 }) # sapply END 	 
     
     ## Replacing the NaNs by 'NA.
     ## NaN's are obtained when using theFUN=mean with complete NA values
     nan.index <- which(is.nan(z))  
     
     if (length(nan.index) > 0  )                           
       z[nan.index] <- NA
       
     # Getting the position of all the years in which there were no values
     # min(NA:NA, na.rm=TRUE) == Inf  ; max(NA:NA, na.rm=TRUE) == -Inf
     inf.index <- which(is.infinite(z))
  
     # Changing all the Inf and -Inf by NA's
     if ( length(inf.index) > 0 ) { z[inf.index] <- NA } 
	 
	 names(z) <- seasons
	 
	 return(z)

} # 'seasonalfunction.default' end



# 'dates'   : "numeric", "factor", "Date" indicating how to obtain the 
#             dates for correponding to the 'sname' station
#             If 'dates' is a number, it indicates the index of the column in 
#                'x' that stores the dates
#             If 'dates' is a factor, it have to be converted into 'Date' class, 
#                using the date format  specified by 'date.fmt'
#             If 'dates' is already of Date class, the following line verifies that
#                the number of days in 'dates' be equal to the number of element in the 
#                time series corresponding to the 'st.name' station
# 'date.fmt': format in which the dates are stored in 'dates'.
#             ONLY required when class(dates)=="factor" or "numeric"
# 'out.type': string that define the desired type of output. Possible values are
#             -) "data.frame": a data.frame, with as many columns as stations 
#                              are included in 'x', and an additional column indicating the Year
#             -) "db"        : a data.frame, with 3 colums will be produced.
#                              The first column will store the Year, 
#                              The second column will store the ID of the station,
#                              The third column will contain the seasonal 
#                                value corresponding to that year and that station.
# 'verbose'      : logical; if TRUE, progress messages are printed 
seasonalfunction.data.frame <- function(x, FUN, na.rm=TRUE,
                                        dates, date.fmt="%Y-%m-%d", 
								        out.type="data.frame", 
									    verbose=TRUE,...) {
  
  # Checking that the user provied a valid argument for 'out.type'  
  if (is.na(match( out.type, c("data.frame", "db") ) ) ) 
      stop("Invalid argument: 'out.type' must be in c('data.frame', 'db'")
	  
  # Checking that the user provied a valid argument for 'FUN'       
   if (missing(FUN))
         stop("Missing argument: 'FUN' must be provided")
	  
  # Checking that the user provied a valid argument for 'dates'       
  if (missing(dates)) {
      stop("Missing argument: 'dates' must be provided") 
  } else  
    {    
       # Checking that the user provied a valid argument for 'dates'  
       if (is.na(match(class(dates), c("numeric", "factor", "Date")))) 
           stop("Invalid argument: 'dates' must be of class 'numeric', 'factor', 'Date'")
           
        # Verification that the number of days in 'dates' be equal to the number 
        # of elements in 'x'
        if ( ( class(dates) == "Date") & (length(dates) != nrow(x) ) ) 
        stop("Invalid argument: 'length(dates)' must be equal to 'nrow(x)'") 
    }  # ELSE end
        
  # If 'dates' is a number, it indicates the index of the column of 'x' that stores the dates
  # The column with dates is then substracted form 'x' for easening the further computations
  if ( class(dates) == "numeric" ) {    
    tmp   <- dates
    dates <- as.Date(x[, dates], format= date.fmt)
    x     <- x[-tmp]
  }  # IF end 
  
  # If 'dates' is a factor, it have to be converted into 'Date' class, 
  # using the date format  specified by 'date.fmt'
  if ( class(dates) == "factor" ) dates <- as.Date(dates, format= date.fmt)
     
  # Amount of stations in 'x'
  nstations <- ncol(x)

  # ID of all the stations in 'x'
  snames <- colnames(x)  
  
  # Computing the Starting and Ending Year of the analysis
  Starting.Year <- as.numeric(format(range(dates)[1], "%Y"))
  Ending.Year   <- as.numeric(format(range(dates)[2], "%Y"))
  
  # Amount of Years belonging to the desired period
  nyears <- Ending.Year - Starting.Year + 1
  
  # Amount of months belonging to the desired period
  nmonths <- 12*nyears
  
  # Requiring the Zoo Library (Z’s ordered observations)
  require(zoo)
  
  if (verbose) print("Starting the computations...", quote=FALSE )
  
  seasons <- c("DJF", "MAM", "JJA", "SON")
  
  if (out.type == "data.frame") {    
  
	# Vector with the names of the field that will be used for storing the results
	field.names <- seasons

	# Creating the data.frame that will store the computed averages for each station
	z <- as.data.frame(matrix(data = NA, nrow = nstations, ncol = 4, 
						byrow = TRUE, dimnames = NULL) )
	  
	z <- sapply(1:nstations, function(j,y) {
		
		if (verbose) print( paste("Station: ", format(snames[j], width=10, justify="left"),
					              " : ",format(j, width=3, justify="left"), "/", 
					              nstations, " => ", 
					              format(round(100*j/nstations,2), width=6, justify="left"), 
					              "%", sep=""), quote=FALSE )
						
		# Transforming the column of 'x' into a zoo object, 
		# using the dates provided by the user
		tmp <- vector2zoo(x=y[,j], dates=dates, date.fmt=date.fmt)
		
		# Computing the annual values
		z[j,] <- seasonalfunction.default(x= tmp, FUN=FUN, na.rm=na.rm)
						 
	}, y = x) # sapply END
    
    z <- t(z) # I don't know WHY !!
    rownames(z) <- snames 
  
  } else if (out.type == "db") {
  
        # Creating a vector with the names of the field that will be used for storing the results
        field.names <- c("StationID", "Season", "Value" )  

        # Creating the data.frame that will store the computed averages for each subcatchment
        z <- as.data.frame(matrix(data = NA, nrow = 4*nstations, ncol = 3, 
                           byrow = TRUE, dimnames = NULL) )
        
        for (j in 1:nstations) {
            
            if (verbose) print( paste("Station: ", format(snames[j], width=10, justify="left"),
                                      " : ", format(j, width=3, justify="left"), "/", 
                                      nstations, " => ", 
                                      format(round(100*j/nstations,2), width=6, justify="left"), 
                                      "%", sep=""), quote=FALSE )
                            
            # Transforming the column of 'x' into a zoo object, 
		    # using the dates provided by the user
		    tmp <- vector2zoo(x=x[,j], dates=dates, date.fmt=date.fmt)
		
		    # Computing the annual values
		    tmp <- seasonalfunction.default(x= tmp, FUN=FUN, na.rm=na.rm)
			
			# Putting the annual/monthly values in the output data.frame
            # The first column of 'x' corresponds to the Year
            row.ini <- (j-1)*4 + 1
            row.fin <-  j*4
            
            z[row.ini:row.fin, 1] <- snames[j] # it is automatically repeted 4 times
            z[row.ini:row.fin, 2] <- seasons        
            z[row.ini:row.fin, 3] <- tmp
                              
        } # FOR end
        
        colnames(z) <- field.names 
  
    } # IF end  
   
  return( z )
  
 } #'seasonalfunction.data.frame' END
 
 
 
 ###########################################################################
# annualfunction: Generic function for computing monthly totals/mean values  #
#               for a zoo object or data.frame                             #
############################################################################
#                  May 15th, 2009; Sep 01st 2009                           #
############################################################################
# 'x   '    :  daily, monthly or annual 'zoo' or 'data.frame' object
# 'FUN'      :Function that will be applied to ALL the values in 'x' belonging to each weather season of the year
#             (e.g., Fun can be some of c('mean', 'max', 'min', 'sd'))              
# 'na.rm'    : Logical. Should missing values be removed?
#              TRUE : the annual values are computed considering only those values different from NA
#              FALSE: if there is AT LEAST one NA within a year, the annual values are NA   
annualfunction <- function(x, FUN, na.rm=TRUE,...) UseMethod("annualfunction")

annualfunction.default <- function(x, FUN, na.rm=TRUE,...) { 

	 # Checking that 'x' is a zoo object
	 if (is.na(match(class(x), c("zoo")))) 
			stop("Invalid argument: 'x' must be of class 'zoo'") 
			
	 # If the user did not provide a title for the plots, this is created automatically
     if (missing(FUN)) stop("Missing argument: 'FUN' must be provided")     
	   
     # Requiring the Zoo Library
	 require(zoo) 
	 
	 # 'FUN' is first applied to all the values of 'x' belonging to the same year
	 totals <- aggregate( x, by= format( time(x), "%Y" ), FUN=FUN, na.rm= na.rm )
     
     #  'FUN' is applied to all the previously computed annual values to get the final result.
	 totals <- aggregate( totals, by= rep("value", length(totals)), FUN=FUN, na.rm= na.rm )
	 
	 return(totals)

} # 'annualfunction.default' end


############################################################################
# annualmean: Generic function for computing monthly totals/mean values  #
#               for a zoo object or data.frame                             #
############################################################################
#                  May 15th, 2009; Sep 01st 2009                           #
############################################################################
# 'dates'   : "numeric", "factor", "Date" indicating how to obtain the 
#             dates correponding to the 'sname' station
#             If 'dates' is a number, it indicates the index of the column in 
#                'x' that stores the dates
#             If 'dates' is a factor, it have to be converted into 'Date' class, 
#                using the date format  specified by 'date.fmt'
#             If 'dates' is already of Date class, the following line verifies that
#                the number of days in 'dates' be equal to the number of element in the 
#                time series corresponding to the 'st.name' station
# 'date.fmt': format in which the dates are stored in 'dates'.
#             ONLY required when class(dates)=="factor" or "numeric"
# 'out.type': string that define the desired type of output. Possible values are
#             -) "data.frame": a data.frame, with as many columns as stations 
#                              are included in 'x', and an additional column indicating the Year
#             -) "db"        : a data.frame, with 3 colums will be produced.
#                              The first column will store the Year, 
#                              The second column will store the ID of the station,
#                              The third column will contain the seasonal 
#                                value corresponding to that year and that station.
# 'verbose'      : logical; if TRUE, progress messages are printed 
annualfunction.data.frame <- function(x, FUN, na.rm=TRUE,
                                    dates, date.fmt="%Y-%m-%d",
									verbose=TRUE,...) {
	  
  # If the user did not provide a title for the plots, this is created automatically
  if (missing(FUN)) stop("Missing argument: 'FUN' must be provided") 
	  
  # Checking that the user provied a valid argument for 'dates'       
  if (missing(dates)) {
      stop("Missing argument: 'dates' must be provided") 
  } else  
    {    
       # Checking that the user provied a valid argument for 'dates'  
       if (is.na(match(class(dates), c("numeric", "factor", "Date")))) 
           stop("Invalid argument: 'dates' must be of class 'numeric', 'factor', 'Date'")
           
        # Verification that the number of days in 'dates' be equal to the number 
        # of elements in 'x'
        if ( ( class(dates) == "Date") & (length(dates) != nrow(x) ) ) 
        stop("Invalid argument: 'length(dates)' must be equal to 'nrow(x)'") 
    } # ELSE end
        
  # If 'dates' is a number, it indicates the index of the column of 'x' that stores the dates
  # The column with dates is then substracted form 'x' for easening the further computations
  if ( class(dates) == "numeric" ) {    
    tmp   <- dates
    dates <- as.Date(x[, dates], format= date.fmt)
    x     <- x[-tmp]
  }  # IF end 
     
  # Amount of stations in 'x'
  nstations <- ncol(x)

  # ID of all the stations in 'x'
  snames <- colnames(x)  
  
  # Computing the Starting and Ending Year of the analysis
  Starting.Year <- as.numeric(format(range(dates)[1], "%Y"))
  Ending.Year   <- as.numeric(format(range(dates)[2], "%Y"))
  
  # Amount of Years belonging to the desired period
  nyears <- Ending.Year - Starting.Year + 1
  
  # Amount of months belonging to the desired period
  nmonths <- 12*nyears
  
  # Requiring the Zoo Library
  require(zoo)
  
  if (verbose) print("Starting the computations...", quote=FALSE )

  # Creating the data.frame that will store the computed averages for each station
  z <- NA*numeric(nstations)
                        
  names(z) <- snames
      
  z[1:nstations] <- sapply(1:nstations, function(j, y) {
        
      if (verbose) print( paste("Station: ", format(snames[j], width=10, justify="left"),
                                " : ",format(j, width=3, justify="left"), "/", 
                                nstations, " => ", 
                                format(round(100*j/nstations,2), width=6, justify="left"), 
                                "%", sep=""), quote=FALSE )
                            
      # Transforming the column of 'x' into a zoo object, 
      # using the dates provided by the user
      tmp <- vector2zoo(x=y[,j], dates=dates, date.fmt=date.fmt)
            
      # Computing the annual values
      z[j] <- annualfunction.default(x= tmp, FUN=FUN, na.rm=na.rm)
                         
  }, y = x) # sapply END
  
  return(z)
  
} #'annualfunction.data.frame' END



#################################################################
# matrixplot: Plots a color matrix representing the amount of days #
#          with information in a set of gauging stations        #
#################################################################
#        May 21th, 2009; Sep 22th, 2009, Sep 24th, 2010         #
#################################################################
# Adapted (and thank you very much) from: 
# http://www2.warwick.ac.uk/fac/sci/moac/currentstudents/peter_cock/r/matrix_contour/

# 'x'         : variable of type 'matrix', with the amount of days with information in each station
#               -) The rows represent the gauging stations
#               -) The columns represetn the years, and they stores the amount of 
#                  days with information in each station
# 'ColorRamp' :  Character or function defining a personalized color ramp for ploting the maps. \cr
#                Valid character values are in c("Days", "Precipitation", "Temperature", 
#                "PCPAnomaly", "PCPAnomaly2" "TEMPAnomaly", "TEMPAnomaly2", "TEMPAnomaly3").  
# 'ncolors'   :  Number of color intervals that will be used for differentiating
#                from 0 to 366 days with information
# 'main'      :  Main title for the plot
                   
matrixplot <- function(x, ColorRamp="Days", ncolors=70, main="", ...) {

  # Checking that 'class(x)==Date'
  if (is.na(match(class(x), c("matrix", "data.frame") ) ) ) 
     stop("Invalid argument: 'x' must be of class 'matrix' or 'data.frame'")
     
  # If 'x' is a data.frame, it trys to coherce into a matrix
  if (class(x) == "data.frame") x <- as.matrix(x)
  
  # Generate the traspose of the matrix, in order to get the years in the 'x' axis
  # and the stations in the 'y' axis
  #x <- t(x)
           
  Days.cols          <- colorRampPalette(c("red3", "orange", "darkgoldenrod1", "yellow", "darkolivegreen2", "green"))  
  
  #Precipitation.cols <- colorRampPalette(c("aquamarine", "blue", "darkblue"))
  Precipitation.cols <- colorRampPalette(c("red3", "orange", "yellow", "darkolivegreen3", "lightskyblue", "royalblue3"))
  Temperature.cols   <- colorRampPalette(c("yellow4", "yellow", "orange", "red3", "darkred"))
  
  PCPAnomaly.cols    <- colorRampPalette(c("sienna4", "peachpuff", "royalblue", "blue"))
  PCPAnomaly2.cols   <- colorRampPalette(c("darkred", "red3", "orange", "yellow", "lightskyblue", "royalblue3", "darkblue"))

  TEMPAnomaly.cols   <- colorRampPalette(c("lightyellow", "yellow", "red", "darkred"))	
  TEMPAnomaly2.cols  <- colorRampPalette(c("yellow4", "yellow", "orange", "red3", "darkred"))
  TEMPAnomaly3.cols  <- colorRampPalette(c("darkblue", "royalblue3", "lightskyblue", "yellow", "orange", "red3", "darkred"))
 
     
  # Another way for temperature colors, using the reverse order (from white to red):
  # Temperature.cols <- rev(heat.colors(100)) 
  
  # # Generating palettes of colors
  if (class(ColorRamp) != "function"  ) {
     # Checking that the user provided a valid argument for 'ColorRamp'       
    if (is.na(match(ColorRamp, c("Days", "Precipitation", "Temperature", "PCPAnomaly", "PCPAnomaly2", "TEMPAnomaly", "TEMPAnomaly2", "TEMPAnomaly3") ) ) ) {
      stop("Invalid argument: 'ColorRamp' must be in c('Days', 'Precipitation', 'Temperature', 'PCPAnomaly', 'PCPAnomaly2', 'TEMPAnomaly', 'TEMPAnomaly2', 'TEMPAnomaly3')") 
    } else {
      # Assgning the color ramp, when 'ColorRamp' was given as a character
      if (ColorRamp == "Days") {   
      ColorRamp <- Days.cols   
      } else if (ColorRamp == "Precipitation") {   
        ColorRamp <- Precipitation.cols   
        } else if (ColorRamp == "Temperature") {
            ColorRamp <- Temperature.cols  
          } else if (ColorRamp == "PCPAnomaly") {
               ColorRamp <- PCPAnomaly.cols  
            } else if (ColorRamp == "PCPAnomaly2") {
               ColorRamp <- PCPAnomaly2.cols  
              }  else if (ColorRamp == "TEMPAnomaly") {
                  ColorRamp <- TEMPAnomaly.cols  
                } else if (ColorRamp == "TEMPAnomaly2") {
                    ColorRamp <- TEMPAnomaly2.cols 
                  } else if (ColorRamp == "TEMPAnomaly3") {
                      ColorRamp <- TEMPAnomaly3.cols 
                    }# ELSE end 
      } # ELSE end
  } # IF end
     
  #par(fig=c(0,0.8,0,0.8), new=TRUE)      
  require(lattice) # for levelplot()
  y <- levelplot(x, scales=list(tck=0, x=list(rot=90)), 
                 col.regions=ColorRamp(ncolors), 
                 #at= seq(0, 366, length.out= ncolors),
                 main= main, 
                 xlab=NULL, ylab=NULL,...)
                 
  #par(fig=c(0,0.8,0.55,1), new=TRUE)
  #up <- colMeans(x)
  #up.g <- barplot(up)
  
  #par(fig=c(0.65,1,0,0.8),new=TRUE)
  #right <- rowMeans(x)
  #r.g <- barplot(right)
  
  #print(y, split = c(1, 1, 10, 10))
  #print(up.g, split = c(1, 1, 10, 10), newpage = FALSE)
  
  return(y)
} # 'matrixplot' END




######################################################
# fdc: Flow Duration Curve, computation and plot     #
######################################################
#                   June 04, 2009                    #
######################################################

# Plot the flow Duration Curve in the original time units of 'x' and 
# also gives the probability of exceedence of each element

# plot   : logical. Indicates if the flow duration curve should be plotted or not
# thr.shw: logical, indicating if the stremflow values corresponding to the user-defined thresholds 'lQ.thr' and 'hQ.thr'  have to be plotted
# new    : logical, indicates if a new graphics device has to be started
# log    : character, indicates which axis has to be plotted with a logarithmic scale. By default is 'y' 

fdc <-function(x, ...) UseMethod("fdc")

fdc.default <- function (x, 
                         lQ.thr=0.7,
                         hQ.thr=0.2,
                         plot=TRUE, 
                         log="y",
                         main="Flow Duration Curve", 
                         xlab="% Time flow equalled or exceeded", 
                         ylab="Q, [m3/s]", 
                         ylim,
                         col="black", 
                         pch=1, 
                         lwd=1,
                         lty=1, 
                         cex=0.4,
                         leg.txt=NULL,
				         verbose= TRUE, 
                         thr.shw=TRUE,
                         new=TRUE,                         
                         ...) {
                         
     # Returns the position in the vector 'x' where the scalar 'Q' is located
     Qposition <- function(x, Q) {     
       Q.dist  <- abs(x - Q)
       Q.index <- which.min( Q.dist ) 
       return(Q.index)       
     } # end
     
     if (log == "y") {
       x.zero.index <- which(x==0)
       if (length(x.zero.index) > 0 ) { 
        x <- x[-x.zero.index] 
        if (verbose) print("[Warning: all 'x' equal to zero will not be plotted]", quote=FALSE)
       } # IF end      
     } # IF end

	 # If 'x' is of class 'ts' or 'zoo'
	 #if ( !is.na( match( class(x), c("ts", "zoo") ) ) ) 
	 x <- as.numeric(x)
	 
	 # Storing the original values
	 x.old <- x
	 
	 # 1) Sort 'x' in drecreasing order. This is just for avoiding misleading 
	 #lines when using 'type="o"' for plotting
	 x <- sort(x)
	 
	 # Index with the position of the original values
	 ind <- match(x.old, x)
	 
	 # 2) Compute the length of 'x'
	 n <- length(x)
	 
	 # 3) Creation of the output vector
	 dc <- rep(NA, n)
	 
	 # 4) Exceedence Probability
	 dc[1:n] <- sapply(1:n, function(j,y) {									  
						  dc[j] <- length( which(y >= y[j]) ) 
					  }, y = x)
					  
	 # Computing the probabilitites
	 dc <- dc/n	 
	 
	 # Another way
	 # Fn <- ecdf(x)
	 # dc <- 1 - Fn(x) + 1/n
	 
	 if (plot) {
     
          if ( missing(ylim) ) {
            ylim <- range(x)
          } # IF end
     
          # If a new plot has to be created
          if (new) {          
               
            # Creating the plot, but without anything on it, for allowing the call to polygon
             if (log=="y") {
               plot(dc, x,  xaxt = "n", yaxt="n", type="o", pch=pch, col=col, lty=lty, 
                    cex=cex, main=main, xlab=xlab, ylab=ylab, ylim=ylim, log=log,...)
             } else {
               plot(dc, x,  xaxt = "n", type="o", pch=pch, col=col, lty=lty, 
                    cex=cex, main=main, xlab=xlab, ylab=ylab, ylim=ylim, log=log,...)
               } # ELSE end 
                  
          } else {
             lines(dc, x,  xaxt = "n", type="o", pch=pch, col=col, lty=lty, 
                   cex=cex)
            } # ELSE end
          
          # Draws the labels corresponding to Annual ticks in the X axis
          Axis(side = 1, at = seq(0.0, 1, by=0.05), labels = FALSE)
          Axis(side = 1, at = seq(0.0, 1, by=0.1), 
               labels = paste(100*seq(0.0, 1, by=0.1),"%", sep="") )
               
          if (log=="y") {
            # Draws the labels corresponding to Annual ticks in the Y axis
            ylabels <- union( c(0.01, 0.1, 1,10), pretty(ylim) )
            Axis( side = 2, at =ylabels, labels = ylabels )
          } # IF end
               
          # If the user provided a value for 'lQ.thr', a vertical line is drawn
          if ( !is.na(lQ.thr) ) {
            abline(v=lQ.thr, col="grey", lty=3, lwd=2)
          } # IF end
      
          # If the user provided a value for 'hQ.thr', a vertical line is drawn
          if ( !is.na(hQ.thr) ) {
            abline(v=hQ.thr, col="grey", lty=3, lwd=2)
          } # IF end	
          
          # Drawing a legend. bty="n" => no border
          if ( !missing(leg.txt) ) {
           legend("topright", legend=leg.txt, cex=cex*1.5, col=col, lty=lty, pch=pch, bty="n")
          } # IF end   
          
          if (thr.shw) {
              # Finding the flow values corresponding to the 'lQ.thr' and 'hQ.thr' pbb of excedence
              x.lQ <- x[Qposition(dc, lQ.thr)]
              x.hQ <- x[Qposition(dc, hQ.thr)]
     
              legend("bottomleft", c(paste("Qhigh.thr=", round(x.hQ, 2), sep=""), 
                                     paste("Qlow.thr=", round(x.lQ, 2), sep="") ), 
                     cex=0.7, bty="n" ) #bty="n" => no box around the legend    
          } # IF end
	 } # IF end
	 
	 # Restoring the original positions
	 dc <- dc[ind]
					  
	 return(dc)
  
 
} # 'fdc' END


######################################################################
# fdc.matrix: (ONLY) Plot of Multiple Flow Duration Curves,          #
#                  for comparison                                    #
######################################################################
#                   June 04, 2009                                    #
######################################################################

fdc.matrix <- function (x, 
                        lQ.thr=0.7,
                        hQ.thr=0.2,
                        plot=TRUE,  
                        log="y",
                        main= "Flow Duration Curve",                         
                        xlab="% Time flow equalled or exceeded", 
						ylab="Q, [m3/s]",
                        ylim,
                        col=palette("default")[1:ncol(x)],
                        pch=1:ncol(x),
                        lwd=rep(1, ncol(x)),
                        lty=1:ncol(x), 
                        cex=0.4,
                        leg.txt= colnames(x),
                        verbose=TRUE, 
                        thr.shw=TRUE,
                        new=TRUE,
                        ...) {
      
  n <- ncol(x)
  
  if (missing(ylim)) { ylim <- range(x) }
  
  j <- 1 # starting column for the analysis
  
  if (verbose) print( paste("Column: ", format(j, width=10, justify="left"),
                            " : ", format(j, width=3, justify="left"), "/", 
                            n, " => ", 
                            format(round(100*j/n,2), width=6, justify="left"), 
                            "%", sep=""), quote=FALSE )
   
  # Computing and plotting the Flow duration Curve for the first vector
  fdc(x=x[,1], plot=plot, pch=pch[1], col=col[1], lty=lty[1], 
      cex=cex, main=main, xlab= xlab, ylab=ylab, log=log, verbose=verbose, 
      thr.shw=FALSE, new=TRUE, ...)
   
  # Plotting the Flow Duration Curves
  sapply(2:n, function(j) {
  
        if (verbose) print( paste("Column: ", format(j, width=10, justify="left"),
                                " : ", format(j, width=3, justify="left"), "/", 
                                n, " => ", 
                                format(round(100*j/n,2), width=6, justify="left"), 
                                "%", sep=""), quote=FALSE )
			# Computing and plotting the Flow duration Curve for the first vector
            fdc(x=x[,j], plot=plot, pch=pch[j], col=col[j], lty=lty[j], 
                cex=cex, main=main, xlab= xlab, ylab=ylab, log=log, verbose=verbose, 
                thr.shw=FALSE, new=FALSE, ...)
        } )
                  
  if (plot) { 
      # Drawing a legend. bty="n" => no border
      if ( is.null(colnames(x)) ) {
        leg.txt <- paste("Q", 1:ncol(x), sep="")
      } # IF end
      legend("topright", legend=leg.txt, cex=cex*2.2, col=col, lty=lty, pch=pch, bty="n")
  } # IF end
 
} # 'fdc.matrix' END


fdc.data.frame <- function(x, 
                           lQ.thr=0.7,
                           hQ.thr=0.2,
                           plot=TRUE,  
                           log="y",
                           main= "Flow Duration Curve",                         
                           xlab="% Time flow equalled or exceeded", 
						   ylab="Q, [m3/s]",
                           ylim,
                           col=palette("default")[1:ncol(x)],
                           pch=1:ncol(x),                            
                           lwd=rep(1, ncol(x)),
                           lty=1:ncol(x), 
                           cex=0.4,
                           leg.txt= colnames(x),
                           verbose=TRUE, 
                           thr.shw=TRUE,
                           new=TRUE,
                           ...) {
                  
   x <- as.matrix(x)
   
   NextMethod("fdc", x, 
               lQ.thr=lQ.thr,
               hQ.thr=hQ.thr,
               plot=plot, 
               main=main,                         
               xlab=xlab, 
               ylab=ylab, 
               log=log,
               col=col,
               pch=pch, 
               lty=lty, 
               cex=cex,
               verbose=verbose, 
               leg.txt= leg.txt,
               thr.shw=thr.shw,
               new=new,
               ylim=ylim,
               ...)
                  
} # 'fdc.data.frame' END


######################################################
# fdcu: Flow Duration Curve with Uncertainty Bounds  #
######################################################
#  Started    :  January 29th, 2010                  #
# Last updated:  February 04th, 2010                 #
######################################################

# Plot the flow Duration Curve in the original time units of 'x' and 
# also gives the probability of exceedence of each element

# 'x'            : 'numeric', 'matrix' or 'data.frame' whose columns contains the values of the time series of observed streamflows for which the flow duration curve will be computed.
# 'lband'        : 'numeric', 'matrix' or 'data.frame' whose columns contains the values of the time series with the lower uncertainty bound of 'x'.
# 'uband'        : 'numeric', 'matrix' or 'data.frame' whose columns contains the values of the time series with the upper uncertainty bound of 'x'.
# 'sim'          : 'numeric', 'matrix' or 'data.frame' whose columns contains the values of the time series with the simulated values of 'x', for which the flow duration curve will be computed.
# 'lQ.thr'       : 'numeric', low flows separation threshold. If this value is different from 'NA', a vertical line is drawn in this value, and all the values to the left of it are deemed low flows.
# 'hQ.thr'       : 'numeric', high flows separation threshold. If this value is different from 'NA', a vertical line is drawn in this value, and all the values to the right of it are deemed high flows
# 'plot'         : 'logical'. Indicates if the flow duration curve should be plotted or not
# 'log'          : 'character', indicates which axis has to be plotted with a logarithmic scale. By default is 'y' 
# 'main'         : See '?plot'. An overall title for the plot: see 'title'.
# 'xlab'         : See '?plot'. A title for the x axis: see 'title'.
# 'ylab'         : See '?plot'. A title for the y axis: see 'title'.
# 'ylim'         : See '?plot.default'.  The y limits of the plot.
# 'col'          : See '?plot.default'. The colors for lines and points.  Multiple colors can be specified so that each point can be given its own color.  If there are fewer colors than points they are recycled in the standard fashion. Lines will all be plotted in the first colour specified.
# 'pch'          : See '?plot.default'. A vector of plotting characters or symbols: see 'points'.
# 'lwd'          : See '?plot.default'. The line width, see 'par'.
# 'lty'          : See '?plot.default'. The line type, see 'par'. 
# 'cex'          : See '?plot.default'. A numerical vector giving the amount by which plotting characters and symbols should be scaled relative to the default.  This works as a multiple of 'par("cex")'. 'NULL' and 'NA' are equivalent to '1.0'.  Note that this does not affect annotation
# 'leg.txt'      : vector with the names that have to be used for each column of 'x'.
# 'verbose'      : logical; if TRUE, progress messages are printed 
# 'thr.shw'      : logical, indicating if the stremflow values corresponding to the user-defined thresholds 'lQ.thr' and 'hQ.thr' have to be shown in the plot.
# 'border'       : See '?polygon'. The color to draw the border of the polygon with the uncertainty bounds. The default, 'NA', means to omit borders. 
# 'bands.col'    : See '?polygon'. The color for filling the polygon. The default, 'NA', is to leave polygons unfilled, unless 'density' is specified. If 'bands.density' is specified with a positive value this gives the color of the shading lines.
# 'bands.density': See '?polygon'. The density of shading lines for the polygon with the uncertainty bounds, in lines per inch.  The default value of 'NULL' means that no shading lines are drawn. A zero value of 'bands.density' means no shading nor filling whereas negative values (and 'NA') suppress shading (and so allow color filling).
# 'bands.angle'  : See '?polygon'. The slope of shading lines for the polygon with the uncertainty bounds, given as an angle in degrees (counter-clockwise). 
# 'new'          : logical, if TRUE, a new plotting window is created. 


# Example:
# data(EbroQts)
# q <- sname2ts(EbroQts, "Q071", dates=1)
# q <- window(q, end=as.Date("1965-1-10"))
# lband <- q-min(q, na.rm=T)
# uband <- q+mean(q, na.rm=T)
# fdcu(q, lband, uband)

fdcu <-function(x, lband, uband, ...) UseMethod("fdcu")
                          
fdcu.default <- function (x, 
                          lband,
                          uband,
                          sim,
                          lQ.thr=0.7,
                          hQ.thr=0.2,
                          plot=TRUE,
                          log="y",
                          main="Flow Duration Curve", 
                          xlab="% Time flow equalled or exceeded", 
                          ylab="Q, [m3/s]", 
                          ylim,
                          col=c("black", "red"), 
                          pch=c(1, 15), 
                          lwd=c(1, 0.8),
                          lty=c(1, 3),                          
                          cex=0.2,  
                          leg.txt= c("Qobs", "Qsim", "95PPU"),
                          verbose= TRUE,                             
                          thr.shw=TRUE,
                          border=NA,
                          bands.col="lightcyan",   
                          bands.density=NULL,
                          bands.angle=45,
                          new=TRUE,                                    
                          ...) {
                          
                          
     #checking that the user provided 'x'
     if (missing(x))  stop("Missing argument: 'x' must be provided")
     
     #checking that the user provided 'lband'
     if (missing(lband))  stop("Missing argument: 'lband' must be provided")
     
     #checking that the user provided 'uband'
     if (missing(uband))  stop("Missing argument: 'uband' must be provided")
     
     # Function for finding the position of 'Q' within 'x'
     Qposition <- function(x, Q) {     
       Q.dist  <- abs(x - Q)
       Q.index <- which.min( Q.dist ) 
       return(Q.index)       
     } # end

	 # In case 'x', 'laband' and 'uband' or 'sim' be of 'zoo' or 'ts' class
     x     <- as.numeric(x)
	 lband <- as.numeric(lband)
	 uband <- as.numeric(uband)     
     if (!missing(sim)) { sim <- as.numeric(sim)   }
     
     # Removing zero values when using a logarithmic scale
     if (log == "y") {
     
       x.zero.index     <- which(x==0)
       lband.zero.index <- which(lband==0)
       uband.zero.index <- which(uband==0)
       if (!missing(sim)) { sim.zero.index <- which(sim==0) }
       
       if (length(x.zero.index) > 0 ) { 
        x <- x[-x.zero.index] 
        if (verbose) print("[Warning: all 'x' equal to zero will not be plotted]", quote=FALSE)
       } # IF end 
       
       if (length(lband.zero.index) > 0 ) { 
        lband <- lband[-lband.zero.index] 
        if (verbose) print("[Warning: all 'lband' equal to zero will not be plotted]", quote=FALSE)
       } # IF end 
       
       if (length(uband.zero.index) > 0 ) { 
        uband <- uband[-uband.zero.index] 
        if (verbose) print("[Warning: all 'uband' equal to zero will not be plotted]", quote=FALSE)
       } # IF end 
       
       if (!missing(sim)) { 
            if (length(sim.zero.index) > 0 ) { 
            sim <- sim[-sim.zero.index] 
            if (verbose) print("[Warning: all 'sim' equal to zero will not be plotted]", quote=FALSE)
           } # IF end
       } # IF end
            
     } # IF end
	 
	 # 1) Sort 'x' in decreasing order. This is just for avoiding misleading 
	 #lines when using 'type="o"' for plotting
	 x.sort     <- sort(x, decreasing=TRUE)
	 lband.sort <- sort(lband, decreasing=TRUE)
	 uband.sort <- sort(uband, decreasing=TRUE)
     if (!missing(sim)) { sim.sort <- sort(sim, decreasing=TRUE) }
	 
	 fdc.x     <- fdc(x.sort, log=log, plot=FALSE,...)
	 fdc.lband <- fdc(lband.sort, log=log, plot=FALSE)
	 fdc.uband <- fdc(uband.sort, log=log, plot=FALSE)
     if (!missing(sim)) { fdc.sim <- fdc(sim.sort, log=log, plot=FALSE) }
	 
	 #na.index <- which(is.na(x))

     # Avoiding plotting the uncertainty bands for the Na's
     #uband[na.index] <- uband[na.index-1]
     #lband[na.index] <- lband[na.index+1]
   
     if (plot) {
     
         if ( missing(ylim) ) {
          ylim <- range(lband, uband, na.rm=TRUE)
         } else {
             # In case 'ylim' is 'NA'
             if ( is.na(ylim[1]) ) { ylim[1] <- range(x, na.rm=TRUE)[1] }
             if ( is.na(ylim[2]) ) { ylim[2] <- range(x, na.rm=TRUE)[2] }
           } # ELSE end
         
         # Logarithmic scale can NOT start in zero
         if (log=="y") { ylim <- c(0.01, ylim[2]) }

         # Creating the 'x' values of the polygons of the bands
         t <- c( fdc.lband, rev(fdc.uband) )
      
         # Creating the 'y' values of the polygons of the bands
         bands <- c( lband.sort, rev(uband.sort) )
         
         if (new) {
             # Creating the plot, but without anything on it, for allowing the call to polygon
             if (log=="y") {
               plot(fdc.x, x.sort, type="n", xaxt = "n", yaxt = "n", main=main, xlab=xlab, ylab=ylab, log=log, ylim=ylim,...)  
             } else {
               plot(fdc.x, x.sort, type="n", xaxt = "n", main=main, xlab=xlab, ylab=ylab, log=log, ylim=ylim,...) 
               } # ELSE end      
         } # IF end
      
         # Plotting the polygons between the lower and upper bands
         polygon(t, bands, col=bands.col, density=bands.density, angle=bands.angle, border=border)   
         
         # Ploting the OBServations over the polygon
         lines(fdc.x, x.sort, cex=cex, col=col[1], pch=pch[1], lwd=lwd[1], lty=lty[1]) 
         points(fdc.x, x.sort, cex=cex, col=col[1], pch=pch[1]) 
         
         # Ploting the SIMulated values over the polygon
         if (!missing(sim)) { 
           lines(fdc.sim, sim.sort, cex=cex, col=col[2], pch=pch[2], lwd=lwd[2], lty=lty[2]) 
           points(fdc.sim, sim.sort, cex=cex, col=col[2], pch=pch[2]) 
         } # IF end
         
         # Draws the labels corresponding to Annual ticks in the X axis
         Axis(side = 1, at = seq(0.0, 1, by=0.05), labels = FALSE)
         Axis(side = 1, at = seq(0.0, 1, by=0.1), labels = paste(100*seq(0.0, 1, by=0.1),"%", sep="") )
         
         if (log=="y") {
           # Draws the labels corresponding to Annual ticks in the Y axis
           ylabels <- union( c(0.01, 0.1, 1,10), pretty(ylim) )
           Axis( side = 2, at =ylabels, labels = ylabels )
         } # IF end
         
         # Drawing a legend. bty="n" => no border
         if (!is.null(leg.txt)) {
         
             # Drawing a legend. bty="n" => no border
             ifelse (log=="y", leg.pos <- "bottomleft", leg.pos <- "topright")
         
             if (!missing(sim)) { #Legend with the OBS + simulations + 95PPU
              legend(leg.pos, legend=leg.txt,  #inset=0.03,
                 bty="n", cex =0.9, col=c(col[1], col[2], bands.col), lwd=c(lwd[1], lwd[2], 0), lty=c(lty[1], lty[2], 0), pch=c(NA,NA,15), pt.cex=3)
             } else { #Legend only with the OBS + 95PPU
              legend(leg.pos, legend=c(leg.txt[1], leg.txt[3]),  #inset=0.03,
                 bty="n", cex =0.9, col=c(col[1], bands.col), lwd=c(lwd[1], 0), lty=c(lty[1], 0), pch=c(NA,15), pt.cex=3)
             }# IF end 
                
         } # IF end     
         
          # If the user provided a value for 'lQ.thr', a vertical line is drawn
         if ( !is.na(lQ.thr) ) {
            abline(v=lQ.thr, col="grey", lty=3, lwd=2)
         } # IF end
    
         # If the user provided a value for 'hQ.thr', a vertical line is drawn
         if ( !is.na(hQ.thr) ) {
            abline(v=hQ.thr, col="grey", lty=3, lwd=2)
         } # IF end	          
         
         # If the user want to see the Q values corresponding to 'lQ.thr' and 'hQ.thr'
         if (thr.shw) {
            # Drawing a legend. bty="n" => no border
             ifelse (log=="y", leg.pos <- "topright", leg.pos <- "bottomleft")
             
            # Finding the flow values corresponding to the 'lQ.thr' and 'hQ.thr' pbb of excedence
            x.lQ <- x.sort[Qposition(fdc.x, lQ.thr)]
            x.hQ <- x.sort[Qposition(fdc.x, hQ.thr)]
         
            legend(leg.pos, c(paste("Qhigh.thr=", round(x.hQ, 2), sep=""), 
                                   paste("Qlow.thr=", round(x.lQ, 2), sep="") ), 
                   cex=0.7, bty="n" ) #bty="n" => no box around the legend 
         } # IF end
         
     } # IF 'plot' end    
	 
	 #tmp <- cbind(fdc.lband, fdc.x, fdc.uband) 
	 #colnames(tmp) <- c("fdc.l", "fdc.x", "fdc.u")    
	 
	 #return(tmp)
 
} # 'fdcu.default' END



######################################################################
# fdcu.matrix: (ONLY) Plots of Multiple Flow Duration Curves,        #
#                  with uncertainty bands                            #
######################################################################
# Started    :  January 29th, 2010                                   #
# Last updated:  April 06th, 2010                                    #
######################################################################

fdcu.matrix <- function (x, 
                         lband,
                         uband,                         
                         sim,
                         lQ.thr=0.7,
                         hQ.thr=0.2,
                         plot=TRUE,
                         log="y",
                         main="Flow Duration Curve", 
                         xlab="% Time flow equalled or exceeded", 
                         ylab="Q, [m3/s]", 
                         ylim,
                         col=matrix(c(rep("black", ncol(x)), palette("default")[2:(ncol(x)+1)]), byrow=FALSE, ncol=2),
                         pch=matrix(rep(c(1, 15), ncol(x)), byrow=TRUE, ncol=2),
                         lwd=matrix(rep(c(1, 0.8), ncol(x)), byrow=TRUE, ncol=2),
                         lty=matrix(rep(c(1, 3), ncol(x)), byrow=TRUE, ncol=2),                        
                         cex=rep(0.1, ncol(x)),   
                         leg.txt=c("OBS", colnames(x), "95PPU"),
                         verbose= TRUE,
                         thr.shw=TRUE,
                         border=rep(NA, ncol(x)),
                         bands.col=rep("lightcyan", ncol(x)),
                         bands.density=rep(NULL, ncol(x)),
                         bands.angle=rep(45, ncol(x)),                         
                         new=TRUE,       
                         ...) {
  n <- ncol(x)
  
  j <- 1 # starting column for the analysis
  
  if (verbose) print( paste("Column: ", format(j, width=10, justify="left"),
                      " : ", format(j, width=3, justify="left"), "/", 
                      n, " => ", 
                      format(round(100*j/n,2), width=6, justify="left"), 
                      "%", sep=""), quote=FALSE )
   
  # Computing and plotting the Flow duration Curve for the first column
  if (missing(sim)) {
  
      sim.exists <- FALSE
      fdcu(x=x[,1], lband=lband[,1], uband=uband[,1], lQ.thr=lQ.thr, hQ.thr=hQ.thr, 
           plot=TRUE, main=main, xlab= xlab, ylab=ylab, log=log, col=col[1,], pch=pch[1,], 
           lwd=lwd[1,], lty=lty[1,], cex=cex[1], verbose, leg.txt=NULL, thr.shw=thr.shw, border=border[1], 
           bands.col= bands.col[1], bands.density=bands.density[1], bands.angle=bands.angle[1], new=TRUE, ...)
  } else {
  
    sim.exists <- TRUE
    fdcu(x=x[,1], lband=lband[,1], uband=uband[,1], sim=sim[,1], lQ.thr=lQ.thr, hQ.thr=hQ.thr, 
           plot=TRUE, main=main, xlab= xlab, ylab=ylab, log=log, col=col[1,], pch=pch[1,], 
           lwd=lwd[1,], lty=lty[1,], cex=cex[1], verbose, leg.txt=NULL, thr.shw=thr.shw, border=border[1], 
           bands.col= bands.col[1], bands.density=bands.density[1], bands.angle=bands.angle[1], new=TRUE, ...)
    }
   
  # Plotting the Flow Duration Curves
  sapply(2:n, function(j) {
  
         if (verbose) print( paste("Column: ", format(j, width=10, justify="left"),
                             " : ", format(j, width=3, justify="left"), "/", 
                             n, " => ", 
                             format(round(100*j/n,2), width=6, justify="left"), 
                             "%", sep=""), quote=FALSE )
                             
            if (!sim.exists) {                 
              # Computing and plotting the Flow duration Curve for the first vector
              fdcu(x=x[,j], lband=lband[,j], uband=uband[,j], lQ.thr=lQ.thr, hQ.thr=hQ.thr, 
                   plot=TRUE, main=main, xlab= xlab, ylab=ylab, log=log, col=col[j,], pch=pch[j,], 
                   lwd=lwd[j,], lty=lty[j,], cex=cex[j], verbose, leg.txt=NULL, thr.shw=FALSE, border=border[j], 
                   bands.col= bands.col[j], bands.density=bands.density[j], bands.angle=bands.angle[j], new=FALSE, ...)
            } else {
              fdcu(x=x[,j], lband=lband[,j], uband=uband[,j], sim=sim[,j], lQ.thr=lQ.thr, hQ.thr=hQ.thr, 
                   plot=TRUE, main=main, xlab= xlab, ylab=ylab, log=log, col=col[j,], pch=pch[j,], 
                   lwd=lwd[j,], lty=lty[j,], cex=cex[j], verbose, leg.txt=NULL, thr.shw=FALSE, border=border[j], 
                   bands.col= bands.col[j], bands.density=bands.density[j], bands.angle=bands.angle[j], new=FALSE, ...)
              } # ELSE end
    } )
    
    if (verbose) print("Re-plotting the 'sim' lines")         
    
    # Re-Plotting the lines that were overdrawn by the polygons
    sapply(1:n, function(j) {                
            
            # Ploting the line with observations over the uncertainty polygon            
            tmp <- sort(x[,j], decreasing=TRUE)
            if (log == "y") {
               x.zero.index <- which(tmp==0)
               if (length(x.zero.index) > 0 ) { 
                tmp <- tmp[-x.zero.index] 
               } # IF end      
            } # IF end
            
            tmp.fdc <- fdc(tmp, plot=FALSE, log=log)
            lines( tmp.fdc, tmp, cex=cex, col=col[j, 1], pch=pch[j, 1], lwd=lwd[j, 1], lty=lty[j, 1] )
            points(tmp.fdc, tmp, cex=cex, col=col[j, 1], pch=pch[j, 1]) 
            
            # Ploting the lines with simulations over the uncertainty polygon    
            if (sim.exists) {      
              tmp <- sort(sim[,j], decreasing=TRUE)   
              if (log == "y") {
               x.zero.index <- which(tmp==0)
               if (length(x.zero.index) > 0 ) { 
                tmp <- tmp[-x.zero.index] 
               } # IF end      
             } # IF end  
            tmp.fdc <- fdc(tmp, plot=FALSE, log=log)      
            lines(tmp.fdc, tmp, cex=cex, col=col[j, 2], pch=pch[j, 2], lwd=lwd[j, 2], lty=lty[j, 2]) 
            points(tmp.fdc, tmp, cex=cex, col=col[j, 2], pch=pch[j, 2]) 
            } 
    } ) # sapply end
    
    # Drawing a legend. bty="n" => no border
    if (!is.null(leg.txt)) {
     
        # Legend with bold font 
        par(font=2)
    
        # Drawing a legend. bty="n" => no border
        ifelse (log=="y", leg.pos <- "bottomleft", leg.pos <- "topright")
    
        if (!missing(sim)) { #Legend with the OBS + SIMs + 95PPU
          legend(leg.pos, legend=leg.txt,  inset=0.01,
                 bty="n", col=c(col[1,1], col[,2], bands.col), lwd=c(3*lwd[1,1], 
                 3*lwd[,2], 0), lty=c(lty[1,1], lty[,2], 0), 
                 pch=c(rep(NA, (ncol(x)+1)), 15), pt.cex=2.5, cex=0.7)
                 
        } else { #Legend only with the OBS + 95PPU
         legend(leg.pos, legend=c(leg.txt[1], leg.txt[3]),  inset=0.01,
                bty="n", col=c(col[1,1], bands.col), lwd=c(3*lwd[1,1], 0), 
                lty=c(lty[1,1],0), pch=c(NA,15), pt.cex=2.5, cex=0.7)
        }# IF end    
        
    } # IF end
  
 
} # 'fdcu.matrix' END


fdcu.data.frame <- function(x, 
                         lband,
                         uband,
                         sim,
                         lQ.thr=0.7,
                         hQ.thr=0.2,
                         plot=TRUE,
                         log="y",
                         main="Flow Duration Curve", 
                         xlab="% Time flow equalled or exceeded", 
                         ylab="Q, [m3/s]", 
                         ylim,
                         col=matrix(c(rep("black", ncol(x)), palette("default")[2:(ncol(x)+1)]), byrow=FALSE, ncol=2),
                         pch=matrix(rep(c(1, 15), ncol(x)), byrow=TRUE, ncol=2),
                         lwd=matrix(rep(c(1, 0.8), ncol(x)), byrow=TRUE, ncol=2),
                         lty=matrix(rep(c(1, 3), ncol(x)), byrow=TRUE, ncol=2),                        
                         cex=rep(0.1, ncol(x)),  
                         leg.txt=c("OBS", colnames(x), "95PPU"),
                         verbose= TRUE, 
                         thr.shw=TRUE,
                         border=rep(NA, ncol(x)),
                         bands.col=rep("lightcyan", ncol(x)),
                         bands.density=rep(NULL, ncol(x)),
                         bands.angle=rep(45, ncol(x)),                         
                         new=TRUE,                         
                         ...) {
                  
   x <- as.data.frame(x)
   
   NextMethod("fdcu", x, 
                lband,
                uband,
                sim,
                lQ.thr=lQ.thr,
                hQ.thr=hQ.thr,
                plot=plot,
                log=log,
                main=main, 
                xlab=xlab, 
                ylab=ylab, 
                ylim=ylim,
                col=col,
                pch=pch, 
                lty=lty, 
                cex=cex,  
                leg.txt=leg.txt,
                verbose= verbose,  
                border=border,
                bands.col=bands.col,   
                bands.density=bands.density,
                bands.angle=bands.angle,                   
                ...)
                  
} # 'fdcu.data.frame' END



##############################################################################
# dmc: Monthly Double-Mass Curve for daily precipitation or streamflow data  #
#      From daily time series, in a data.frame, it computes the monthly mean #
#      double-mass curves (Homogeneity test)                                 #
##############################################################################
#                   July 28th, 2009                                          #
##############################################################################
# 'x   '    : variable of type 'zoo' or 'data.frame'  
# 'target'  : character with the ID of the target station, 
#             It has to correspond to some of the column names in 'x'
#             It also can take tha value "all", in which case the double-mass curve is computed for all the stations in 'x'
# 'dates'   : "numeric", "factor", "Date" indicating how to obtain the 
#             dates for correponding to the 'sname' station
#             If 'dates' is a number, it indicates the index of the column in 
#                'x' that stores the dates
#             If 'dates' is a factor, it have to be converted into 'Date' class, 
#                using the date format  specified by 'date.fmt'
#             If 'dates' is already of Date class, the following line verifies that
#                the number of days in 'dates' be equal to the number of element in the 
#                time series corresponding to the 'st.name' station
# 'date.fmt': character indicating the format in which the dates are stored in 'dates'.
#             ONLY required when class(dates)=="factor" or "numeric"
# 'var.type' : character representing the type of variable being plotted 
#              ONLY for determining the function used for computing the 
#              Monthly values when 'FUN' is missing
#              Valid values are:
#              -) "Precipitation" => FUN = sum
#              -) "Flow"          => FUN = mean 
# 'FUN'      : ONLY required when 'var.type' is missing
#              Function that have to be applied for transforming from daily to monthly or annual time step
#              For precipitation FUN MUST be "sum"
#              For temperature and flow time series, FUN MUST be "mean"#             
# 'na.rm'    : Logical. Should missing values be removed?
#              TRUE : the annual values are computed considering only those values different from NA
#              FALSE: if there is AT LEAST one NA within a year, the monthly values are NA  
# 'method'  : See '?cor'
#             a character string indicating which correlation coefficient
#             (or covariance) is to be computed.  
#             Valid values are: '"pearson"', (default), '"kendall"', or '"spearman"', can be abbreviated
# 'use'      : See '?cor'
#              an optional character string giving a method for computing
#              covariances in the presence of missing values.  
#             This must be (an abbreviation of) one of the strings 
#             '"everything"', '"all.obs"', '"complete.obs"', '"na.or.complete"', or '"pairwise.complete.obs"'            
# 'print.out': character. 
#              Valid values are: 
#              -) "data.frame" : a data.frame with the results (monthly values in the reference and target stations, 
#                                and cumulative values in reference and target stations) are pinted out
#              -) "plot"       : only the plot with the double-mass curve is printeed out, NO data.frame
#              -) "both"       : a data.frame with the results and a plot is printed out. Equivalent to 'print.out' = "data.frame" + "both"

.dmc <-function(x, ...) UseMethod(".dmc")

.dmc.default <- function (x, trgt, dates, date.fmt="%Y-%m-%d", var.type, FUN, na.rm=TRUE, 
                 method="pearson", use="pairwise.complete.obs",  
                 main= "Monthly Double-Mass Curve", 
                 screen=c( ceiling(sqrt(ncol(x)-1)), (ncol(x)-1)-ceiling(sqrt(ncol(x)-1)) ),                 
                 xlab="Target Station", 
                 ylab="Reference Station", 
                 col="blue", 
                 print.out="both",...) {
                 
  
       # Checking that the user provied a valid argument for 'target'
  if ( is.na( match( class(trgt), c("character") ) ) ) 
      stop("Invalid argument: 'class(trgt)' must be in c('character')")
      
  # Checking that the user provied a valid argument for 'var.type'       
  if (missing(var.type)) {
     # If the user did not provide a title for the plots, this is created automatically
     if (missing(FUN)) { stop("Missing argument: 'var.type' OR 'FUN' & 'na.rm' must be provided")  }
  } else # If 'var.type' is provided
      # Checking that the user provied a valid argument for 'var.type'       
      if (is.na(match(var.type, c("Precipitation", "Flow") ) ) ) {
            stop("Invalid argument: 'var.type' must be in c('Precipitation', 'Flow')") 
        } else {
               if (var.type=="Precipitation") { FUN= sum  } 
               else if (var.type=="Flow") { FUN= mean  }
               } #ELSE end 
               
  # Checking a valid value for 'print.out'
  if (is.na(match(print.out, c("data.frame", "plot", "both") ) ) ) 
     stop("Invalid argument: 'print.out' must be in c('data.frame', 'plot', 'both')")  
     
  # Checking a valid value for 'method'
  if (is.na(match(method, c("pearson", "kendall", "spearman") ) ) ) 
     stop("Invalid argument: 'method' must be in c('pearson', 'kendall', 'spearman')") 
     
  # Checking a valid value for 'use'
  if (is.na(match(use, c("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs") ) ) ) 
     stop("Invalid argument: 'method' must be in c('everything', 'all.obs', 'complete.obs', 'na.or.complete', 'pairwise.complete.obs')") 
               
  # Checking that the user provied a valid argument for 'dates'  
  if (is.na(match(class(dates), c("numeric", "factor", "Date")))) 
      stop("Invalid argument: 'dates' must be of class 'numeric', 'factor', 'Date'")
        
  # If 'dates' is a number, it indicates the index of the column of 'x' that stores the dates
  # The column with dates is then substracted from 'x' for easening the further computations
  if ( class(dates) == "numeric" ) {    
    tmp   <- dates
    dates <- as.Date(x[, dates], format= date.fmt)
    x <- x[-tmp]
  }  # IF end 
  
  # If 'dates' is a factor, it have to be converted into 'Date' class, 
  # using the date format  specified by 'date.fmt'
  if ( class(dates) == "factor" ) dates <- as.Date(dates, format= date.fmt)
  
  # If 'dates' is already of Date class, the following line verifies that
  # the number of days in 'dates' be equal to the number of element in the 
  # time series corresponding to the 'st.name' station
  if ( ( class(dates) == "Date") & (length(dates) != nrow(x) ) ) 
     stop("Invalid argument: 'length(dates)' must be equal to 'nrow(x)'") 
  
      col.index <- which(colnames(x) == trgt) 
      
      if ( length(col.index) == 0 ) { stop(paste("Invalid argument: '", trgt, "' is not a colum name in 'x'"), sep="" )
      } else {
          # Assigning the daily values of the target station
          x.target <- x[, col.index]
      
          # The column with daily values on the target station is then substracted from 'x' for easening the further computations
          x <- x[-col.index]
        } # ELSE end        
        
      # computing the number of stations
      nstations <- ncol(x)        
      
      # Transforming the vector of time series ('x.target') and the vector with dates ('dates')
      # into a zoo variable, using the format specified by 'date.fmt'
      x.target  <- vector2zoo(x= x.target, dates= dates, date.fmt= date.fmt)
      
      # Computing the monthly values for the target station
      x.target.monthly <- daily2monthly(x.target, FUN=FUN, na.rm= na.rm)
      
      # Computing the cumulative values, for each month, on the target station
      x.target.cum <- cumsum(x.target.monthly) 
      
      # Computing the daily mean signal on the reference station, skipping the missing values
      x.ref <- rowMeans(x, na.rm=TRUE)
      
      # Transforming the vector of time series ('x.ref') and the vector with dates ('dates')
      # into a zoo variable, using the format specified by 'date.fmt'
      x.ref  <- vector2zoo(x= x.ref, dates= dates, date.fmt= date.fmt)
      
      # Computing the monthly values for the reference station
      x.ref.monthly <- daily2monthly(x.ref, FUN=FUN, na.rm= na.rm)
       
      # Computing the cumulative values, for each month, on the reference station
      x.ref.cum <- cumsum(x.ref.monthly) 
      
      # Pearson's product-moment correlation coefficient between the monthly series
      r <- cor(x.ref.cum, x.target.cum, method= method, use=use)
      
      # Coefficient of cdetermination between the 2 monthly series
      r2 <- r^2
      
      #if (!is.na(match(print.out, c("plot", "both")))) {
          ## Plotting the Double-Mass Curve      
          #plot(x.target.cum, x.ref.cum, type="o", cex=0.6, main=main, xlab=xlab, ylab=ylab, col=col)
                          
          ## Drawing a legend
          #legend("topleft", legend=paste("R2=", round(r2,3), sep=""), cex=0.9)
      #} # IF end
      
      if (!is.na(match(print.out, c("data.frame", "both")))) {
        # Creating a data.frame with the results
        res <- data.frame(Date=time(x.target.cum), Tar.Monthly=x.target.monthly, 
                          Ref.Monthly=x.ref.monthly, Tar.Cum=as.numeric(x.target.cum), 
                          Ref.Cum=as.numeric(x.ref.cum) )                      
        return(res)  
      } # IF end      
      
  } # 'dmc.default' END
  



.dmc.data.frame <- function (x, trgt, dates, date.fmt, var.type, FUN, na.rm=TRUE, 
                            method="pearson", use="pairwise.complete.obs",  
                            main= "Monthly Double-Mass Curve", 
                            screen=c( ceiling(sqrt(ncol(x)-1)), (ncol(x)-1)-ceiling(sqrt(ncol(x)-1)) ),                 
                            xlab="Target Station", 
                            ylab="Reference Station", 
                            col="blue", 
                            print.out="both",...) {
                   
      
  # Checking that the user provied a valid argument for 'x'
  if ( is.na( match( class(x), c("data.frame") ) ) ) 
      stop("Invalid argument: 'class(x)' must be in c('data.frame')")

  # Checking that the user provied a valid argument for 'trgt'
  if ( is.na( match( class(trgt), c("character") ) ) ) 
      stop("Invalid argument: 'class(trgt)' must be in c('character')")
      
  # Checking that the user provied a valid argument for 'var.type'       
  if (missing(var.type)) {
     # If the user did not provide a title for the plots, this is created automatically
     if (missing(FUN)) { stop("Missing argument: 'var.type' OR 'FUN' & 'na.rm' must be provided")  }
  } else # If 'var.type' is provided
      # Checking that the user provied a valid argument for 'var.type'       
      if (is.na(match(var.type, c("Precipitation", "Flow") ) ) ) {
            stop("Invalid argument: 'var.type' must be in c('Precipitation', 'Flow')") 
        } else {
               if (var.type=="Precipitation") { FUN= sum  } 
               else if (var.type=="Flow") { FUN= mean  }
               } #ELSE end 
               
  # Checking a valid value for 'print.out'
  if (is.na(match(print.out, c("data.frame", "plot", "both") ) ) ) 
     stop("Invalid argument: 'print.out' must be in c('data.frame', 'plot', 'both')")  
     
  # Checking a valid value for 'method'
  if (is.na(match(method, c("pearson", "kendall", "spearman") ) ) ) 
     stop("Invalid argument: 'method' must be in c('pearson', 'kendall', 'spearman')") 
     
  # Checking a valid value for 'use'
  if (is.na(match(use, c("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs") ) ) ) 
     stop("Invalid argument: 'method' must be in c('everything', 'all.obs', 'complete.obs', 'na.or.complete', 'pairwise.complete.obs')") 
               
  # Checking that the user provied a valid argument for 'dates'  
  if (is.na(match(class(dates), c("numeric", "factor", "Date")))) 
      stop("Invalid argument: 'dates' must be of class 'numeric', 'factor', 'Date'")
        
  # If 'dates' is a number, it indicates the index of the column of 'x' that stores the dates
  # The column with dates is then substracted from 'x' for easening the further computations
  if ( class(dates) == "numeric" ) {    
    tmp   <- dates
    dates <- as.Date(x[, dates], format= date.fmt)
    x <- x[-tmp]
  }  # IF end 
  
  # If 'dates' is a factor, it have to be converted into 'Date' class, 
  # using the date format  specified by 'date.fmt'
  if ( class(dates) == "factor" ) dates <- as.Date(dates, format= date.fmt)
  
  # If 'dates' is already of Date class, the following line verifies that
  # the number of days in 'dates' be equal to the number of element in the 
  # time series corresponding to the 'st.name' station
  if ( ( class(dates) == "Date") & (length(dates) != nrow(x) ) ) 
     stop("Invalid argument: 'length(dates)' must be equal to 'nrow(x)'")  
  
  
  # If the user selected 1 single 'target' station, the colum index corresponding to that station is looked for
  if ( trgt != "all" ) {
  
        # Computes and plots the double-mass curve for the target station
        .dmc.default(x, trgt=trgt, dates=dates, date.fmt=date.fmt, var.type=var.type,
                   FUN=FUN, na.rm=na.rm, method=method, use=use, main=main,
                   xlab=xlab, ylab=ylab, col=col, print.out=print.out)
        
  } else  if ( trgt == "all" ) {
       
          # ID of all the stations in 'x'
          IDs <- colnames(x)
          
          reference.cum <- x*NA
          target.cum    <- x*NA
          
           # Setting up the screen with the amount of rows and columns defined by screen
          par(mfrow=screen)
          
          for (i in 1:length(IDs) ) {
          
            # Computes and plots the double-mass curve for the target station
            tmp <- .dmc.default(x, trgt=IDs[i], dates=dates, date.fmt=date.fmt, 
                              var.type=var.type, FUN=FUN, na.rm=na.rm, 
                              method=method, use=use, main=paste(main, ", ", IDs[i], sep=""), 
                              xlab=paste(xlab, ", ", IDs[i], sep=""), ylab=ylab, 
                              col=col, print.out=print.out)    
           
            if (i==1) {
              reference.cum <- reference.cum[1:nrow(tmp), ]
              target.cum    <- target.cum[1:nrow(tmp), ]
            }
           
            reference.cum[, i] <- tmp$Ref.Cum
            target.cum[, i]    <- tmp$Tar.Cum      
          } # # FOR end
          
          if (!is.na(match(print.out, c("plot", "both")))) {
             # Adding the date to the reference and target monthly cumulative values
             reference.cum <- cbind(tmp$Date, reference.cum)
             target.cum    <- cbind(tmp$Date, target.cum)
          
             return(list(Ref.Cum=reference.cum, Tar.Cum= target.cum))
          } # IF end
          
         
    } # ELSE end    
 
} # 'dmc.data.frame' END



############################################################
# 'hydropairs' : Visualization of a Correlation Matrix    #
############################################################
#                    July 29th, 2009                       #
############################################################
# On top the (absolute) value of the correlation plus the result of the cor.test as points. #
# On botttom, the bivariate scatterplots, with a fitted line
# On diagonal, histograms (from '?pairs')

# Original idea taken from: http://addictedtor.free.fr/graphiques/graphcode.php?graph=137
# Histogram panles was taken form the R help of the original 'pairs' function

# x     : a numeric vector, matrix or data frame
# dec   : decimals places to be used for showing the correlation values

# use   : an optional character string giving a method for computing
#          covariances in the presence of missing values.  This must be
#          (an abbreviation of) one of the strings '"everything"',
#          '"all.obs"', '"complete.obs"', '"na.or.complete"', or
#          '"pairwise.complete.obs"'.

# method: a character string indicating which correlation coefficient
#          (or covariance) is to be computed.  One of '"pearson"'
#          (default), '"kendall"', or '"spearman"', can be abbreviated.

hydropairs <- function(x, dec=3, use="pairwise.complete.obs", method="pearson",... ) {

  # Checking that the user provied a valid argument for 'x'
  if ( is.na( match( class(x), c("matrix", "data.frame") ) ) ) 
      stop("Invalid argument: 'class(x)' must be in c('data.frame')")

  panel.cor <- function(x, y, digits=dec, prefix="", cex.cor) 
    {
        usr <- par("usr"); on.exit(par(usr)) 
        par(usr = c(0, 1, 0, 1)) 
        
        r <- abs(cor(x, y, method= method, use= use)) 
        
        txt <- format(c(r, 0.123456789), digits=dec)[1] 
        txt <- paste(prefix, txt, sep="") 
        if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 

        test <- cor.test(x,y) 
        # borrowed from printCoefmat
        Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                      symbols = c("***", "**", "*", ".", " ")) 

        text(0.5, 0.5, txt, cex = cex * r) 
        text(.8, .8, Signif, cex=cex, col=2) 
    } # 'panel.cor' END

  panel.hist <- function(x, ...)
         {
             usr <- par("usr"); on.exit(par(usr))
             par(usr = c(usr[1:2], 0, 1.5) )
             h <- hist(x, plot = FALSE)
             breaks <- h$breaks; nB <- length(breaks)
             y <- h$counts; y <- y/max(y)
             rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
         } # 'panel.hist' END    


  # 'font.labels' =2 : bold font for the variables
  # 'cex.labels' controls the size of the fonts
  pairs(x, lower.panel=panel.smooth, upper.panel=panel.cor, 
        diag.panel=panel.hist, ...)
  
} # 'hydropairs' END



####################################################################
# Smry: Min, 1stQ, Mean, Median, 3rdQ, Max, IQR, sd, cv, skewness, #
#       kurtosis, amount of elements and amount of NA's            #
#       For numerical variables.                                   #
#       Skewness and Kurtosis are computed with the e1071 package  #
####################################################################
#	Date: 14-Jun-2008; 11-Sep-2009                                 #
####################################################################
smry <-function(x, ...) UseMethod("smry")


smry.default <- function(x, na.rm=TRUE, digits = max(3, getOption("digits")-3), ...)  {
   
    if ( class(x) %in% c("ts", "zoo") ) { 
        x <- as.numeric(x)
    } # IF end
    
    
    # Creating the resulting object
    z <- matrix(NA, ncol=1, nrow=13)
    
    # If 'x' is not numeric but Date
    if ( class(x) == "Date" ) {
      z <- as.data.frame(z)
    } # IF end
    
    
    if ( class(x) %in% c("numeric", "integer", "Date") ) {
      
        n <- length(x)
        
        na.index <- which(is.na(x))
        nna      <- length(na.index)
      
        if (na.rm) {          
          if (nna > 0) {
            x  <- x[-na.index] }
        } # IF end
        
        s    <- summary(x, ..., digits=digits)  
        
        if ( class(x) %in% c("numeric", "integer") ) {
        
            z[1,1] <- s[1] # min
            z[2,1] <- s[2] # q1
            z[3,1] <- s[3] # median
            z[4,1] <- s[4] # mean
            z[5,1] <- s[5] # q3
            z[6,1] <- s[6] # max
            
            z[7,1] <- IQR(x, na.rm = na.rm)	                             # Interquantile Range IQR = Q(0.75) – Q(0.25)
            z[8,1] <- sd(x, na.rm = na.rm)	                             # Standard Deviation
            z[9,1] <- sd(x, na.rm = na.rm) / abs(mean(x, na.rm = na.rm)) # Coefficient of variation ( coef. of variation = sd / |mean| )
        
            require(e1071) # for the following 2 functions
            z[10,1] <- skewness(x, na.rm = na.rm)  # Skewness (using  e1071 package)   
            z[11,1] <- kurtosis(x, na.rm = na.rm)  # Kurtosis (using  e1071 package)
            
            z <- round( z, digits)
            
        } else { # if 'x' is a Date object
          
            z[1,1] <- as.character(s[1]) # min
            z[2,1] <- as.character(s[2]) # q1
            z[3,1] <- as.character(s[3]) # median
            z[4,1] <- as.character(s[4]) # mean
            z[5,1] <- as.character(s[5]) # q3
            z[6,1] <- as.character(s[6]) # max
            
            z[7,1]  <- NA # IQR
            z[8,1]  <- NA # sd
            z[9,1]  <- NA # cv
            z[10,1] <- NA # sk
            z[11,1] <- NA # kur
            
          } # ELSE end
        
        z[12,1] <- nna # Amount of NA's
        z[13,1] <- n   # Number of elements
    
    } # IF end
      
    row.names(z) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.",
                      "Max.", "IQR", "sd", "cv", "Skewness", "Kurtosis", 
                      "NA's", "n")  
      
    return(z)

} # 'smry.default' end


smry.data.frame <- function(x, na.rm=TRUE, digits = max(3, getOption("digits")-3), ...)  {
        
    # Creating a copy of the original observed values
	z <- as.data.frame( matrix(NA, nrow=13, ncol=ncol(x)) )
	  
	z[,1:ncol(z)] <- sapply(1:ncol(x), function(j,y) {
 
		# Putting the monthly values in the output data.frame
		# The first column of 'x' corresponds to the Year
		z[,j] <- smry.default(x= y[,j], na.rm=na.rm, digits=digits)
						 
	}, y = x) # sapply END
    
    
    rownames(z) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", 
                     "Max.", "IQR", "sd", "cv", "Skewness", "Kurtosis", 
                     "NA's", "n") 
                     
    colnames(z) <- colnames(x)  
    
    return(z)

} # 'smry.data.frame' end




smry.matrix <- function(x, na.rm=TRUE, digits = max(3, getOption("digits")-3), ...)  {
        
    x <- as.data.frame(x)
    smry.data.frame(x, na.rm=na.rm, digits=digits,...)
    #NextMethod("smry", x, na.rm=TRUE, digits=digits,...)

} # 'smry.data.frame' end



#####################################################
# Computes: Hypometric curve corresponding to a DEM #
#####################################################
#   Requires: 'rgdal'                               #
#####################################################
#	Date: 12-Sep-2009, Jul 2010                     #
#####################################################
# Taken from: http://osgeo-org.1803224.n2.nabble.com/hypsometric-integral-from-ecdf-curve-td2231345.html

# 'x'   : 'SpatialGridDataFrame' object with the elevations of the catchment
# Result: Plot in the 'y' axis the elevations of the catchment, and in 
#         the 'x' axis, the percentage of the catchment area that is BELOW a given elevation

# Example:
#require(rgdal)
#dem <- readGDAL(x)
#hypsometric(dem)
hypsometric <- function(x, main="Hypsometric Curve", 
                        xlab="Relative Area above Elevation, (a/A)", 
						ylab="Relative Elevation, (h/H)", col="blue",...) {
  
  if (class(x) != "SpatialGridDataFrame") 
    stop("Invalid argument: 'class(x)' must be 'SpatialGridDataFrame'")	
  
  # Minimum and Maximum elevations in 'dem'
  z.min <- min(x@data, na.rm=TRUE)
  z.max <- max(x@data, na.rm=TRUE)
  
  # Horizontal dimension of the cells of 'x'
  x.dim <- x@grid@cellsize[1]
  # Vertical dimension of the cells of 'x'
  y.dim <- x@grid@cellsize[2]
  
  # Maximum area of the 'dem', in [m2]
  max.area <- length(which(!is.na(x@data$band1))) * x.dim * y.dim
  
  # res$t: elevation values, plus a first value that I don't know what it is
  # res$y: accumulated area BELOW a given elevation value.
  res <- plot.stepfun(ecdf(as.matrix(x)), lwd=0, cex.points=0)  
  
  # Mean elevation in 'dem'
  z.mean.index <- which(round(res$y,3)==0.5)[1]
  z.mean       <- res$t[z.mean.index]
  
  #plot(1 - res$y[-1], res$t[-c(1, length(res$t))], 
  #     main=main, xlim=c(0, 1), 
  #     type="l", ylab=ylab, xlab=xlab, col=col,...)
  
  # res$t[1] represent the minimum elevation within the basin
  # res$y[1] represent the accumulated area below the minimum elevation
  
  # Relative area ABOVE a given elevation
  relative.area <- ( 1 - res$y[-1] )
  
  # Relative elevation  
  relative.elev <- ( res$t[-c(1, length(res$t))] -z.min ) / ( z.max - z.min )
  
  plot(relative.area, relative.elev, 
       xaxt="n", yaxt="n",
       main=main, xlim=c(0, 1), ylim=c(0, 1), 
       type="l", ylab=ylab, xlab=xlab, col=col,...)
       
  # Draws the tickets and labels corresponding to the X axis
  Axis(side = 1, at = seq(0.0, 1, by=0.05), labels = TRUE)
       
  # Draws the tickets and labels corresponding to the Y axis
  Axis(side = 2, at = seq(0.0, 1, by=0.05), labels = TRUE)
  
  # Obtaining a functional form of the spline approximation to the hyposometric integral. 
  # Possible methods are in c("fmm", "periodic", "natural", "monoH.FC")     
  f <- splinefun(relative.area, relative.elev, method="monoH.FC")
  
  # Computing the hypsometric integral
  hi <- integrate(f=f, lower=0, upper=1)
  
  # Drawing a legend with the values of the min, mean, and max elevations
  legend("topright", c(
         paste("Min Elev. :", round(z.min, 2), "[m.a.s.l.]", sep=" "), 
         paste("Mean Elev.:", round(z.mean, 1), "[m.a.s.l.]", sep=" "), 
		 paste("Max Elev. :", round(z.max, 1), "[m.a.s.l.]", sep=" "),
         paste("Max Area  :", round(max.area/1000000, 1), "[km2]", sep=" "),
         "",
         paste("Integral value :", round(hi$value, 3), sep=" "),
         paste("Integral error :", round(hi$abs.error, 3), sep=" ")
         ), 
		 bty="n", cex =0.9, col = c("black","black","black"), 
		 lty=c(NULL,NULL,NULL,NULL) ) #bty="n" => no box around the legend 

} # 'hypsometric' END


    
    
########################################################################
# zoo2RHtest: input file to the 'RHtest_dlyPrcp.r' script              #
########################################################################
#	                   Date: 26-Nov-2009                               #
########################################################################
# function for creating the input file to the 'RHtest_dlyPrcp.r' script,  
# that test the homogeneity of climatological time series (http://ccma.seos.uvic.ca)

# 'x'        : time series that will be written. class(x) must be 'zoo'
# 'fname'    : filename of the output precipitation time series 
# 'tstep.out': Character indicating the time step that have to be used for 
#              writting 'x' into the output file
# 'dec'      : the string to use for decimal points in numeric or complex
#              columns: must be a single character.
# 'na'       : the string to use for missing values in the data

zoo2RHtest <- function(x, fname="pcp.txt", tstep.out="daily", dec=".", na="-999.0") {

  # Checking that 'class(x)' is 'ts' or 'zoo'
  if ( is.na( match(class(x), c("zoo") ) ) )
      stop("Invalid argument: 'class(x)' must be 'zoo'")
      
  # Checking 'tstep.out'
    if ( is.na( match(tstep.out, c("daily", "monthly", "annual") ) ) ) 
      stop("Invalid argument: 'tstep.out' must be in c('daily', monthly', 'annual')")
      
  pfreq <- sfreq(x)
      
  years  <- as.numeric(format(time(x), "%Y"))
  months <- as.numeric(format(time(x), "%m"))
  days   <- as.numeric(format(time(x), "%d"))
  values <- coredata(x)
  
  # Checking the user provide a valid value for 'x'
  if (is.na(match(sfreq(x), c("daily", "monthly", "annual")))) {
		 stop(paste("Invalid argument: 'x' is not a daily, mothly or annual ts, it is a ", sfreq(x), " ts", sep="") ) }	
  
  if ( tstep.out=="daily") {
    out <- data.frame(years=years, months=months, days=days, values=values)
  } else if ( tstep.out=="monthly") {
     out <- data.frame(years=years, months=months, days=0, values=values)
    } else if ( tstep.out=="monthly") {
      out <- data.frame(years=years, months=0, days=0, values=values)
      } # ELSE end
  
  write.table(out, file = fname, sep = " ",
              eol = "\n", na = na, dec = dec, row.names = FALSE,
              col.names = FALSE)
 
 
} # 'zoo2RHtest' end


########################################################################
# dwdays:  average amount of dry/wet days per each month               #
########################################################################
#	                   Date: 24-Jan-2010                               #
########################################################################
# Given a daily time series of precipitation, this function computes the average amount 
# of dry/wet days (pcp > thr or pcp < thr for wet and dry days, respectively) on each month

# 'x'  : zoo. Daily time series of precipitation.
# 'thr': numeric. Value of daily precipitation used as threshold for classifying a day as dry/wet or not. 
#        Days with a precipitation value larger or equal to 'thr' are classified as 'wet days', 
#        whereas precipitation values lower or equal to 'thr' are classified as 'dry days'
# 'type': character, indicating if the daily values have to be calssified as dry or wet days. It works liked to the values specified in 'thr'. Valid values aer c('wet', dry')
dwdays <-function(x, ...) UseMethod("dwdays")
      
dwdays.default <- function(x, thr=0, type="wet", na.rm=TRUE, ... ) { 

  # Requiring the Zoo Library (Zoos ordered observations)
  require(zoo)

  # Checking the user provide a valid value for 'x'
  if (is.na(match(class(x), c("zoo")))) 
        stop("Invalid argument: 'x' must be of class 'zoo'") 
        
  # Checking the user provide a valid value for 'x'
  if (is.na(match(sfreq(x), c("daily")))) {
		 stop(paste("Invalid argument: 'x' is not a daily ts, it is a ", sfreq(x), " ts", sep="") ) }
         
  # Checking the user provide a valid value for 'type'
  if ( is.na(match(type, c("dry", "wet"))) ) 
        stop("Invalid argument: 'type' must be in c('dry', 'wet'") 
		 
  # getting the dates of 'x'
  dates <- time(x)
		 
  # Computing the Starting and Ending Year of the analysis
  Starting.Year <- as.numeric(format(range(dates)[1], "%Y"))
  Ending.Year   <- as.numeric(format(range(dates)[2], "%Y"))
  
  # Total amount of Years of 'x'
  nyears <- Ending.Year - Starting.Year + 1
  
  # Array with the average monthly amount of wet days
  wdays <- rep(NA, 12)
  
  for (m in 1:12) {
    #Extracts all the days of 'x' belonging to the month 'm' 
    pcp.m <- extractzoo(x, m)
    
    if (type=="wet") {    
      wdays[m] <- length(which(pcp.m > thr )) / nyears  
    } else if (type=="dry") {    
        wdays[m] <- length(which(pcp.m < thr )) / nyears  
      } # ELSE end 
  } # FOR m end

  names(wdays) <- month.abb

  return(wdays)         
 
} # 'dwdays.default' end



# 'dates'   : "numeric", "factor", "Date" indicating how to obtain the 
#             dates for correponding to the 'sname' station
#             If 'dates' is a number, it indicates the index of the column in 
#                'x' that stores the dates
#             If 'dates' is a factor, it have to be converted into 'Date' class, 
#                using the date format  specified by 'date.fmt'
#             If 'dates' is already of Date class, the following line verifies that
#                the number of days in 'dates' be equal to the number of element in the 
#                time series corresponding to the 'st.name' station
# 'date.fmt': format in which the dates are stored in 'dates'.
#             ONLY required when class(dates)=="factor" or "numeric"
# 'verbose' : logical; if TRUE, progress messages are printed 
dwdays.data.frame <- function(x, thr=0, type="wet", na.rm=TRUE,
                              dates, 
                              date.fmt="%Y-%m-%d", 
							  verbose=TRUE,...) {
  
  
  # Checking that the user provied a valid argument for 'dates'       
  if (missing(dates)) {
      stop("Missing argument: 'dates' must be provided") 
  } else      
     # Checking that the user provied a valid argument for 'dates'  
     if (is.na(match(class(dates), c("numeric", "factor", "Date")))) 
         stop("Invalid argument: 'dates' must be of class 'numeric', 'factor', 'Date'")
        
  # If 'dates' is a number, it indicates the index of the column of 'x' that stores the dates
  # The column with dates is then substracted form 'x' for easening the further computations
  if ( class(dates) == "numeric" ) {    
    tmp   <- dates
    dates <- as.Date(x[, dates], format= date.fmt)
    x     <- x[-tmp]
  }  # IF end 
  
  # If 'dates' is a factor, it have to be converted into 'Date' class, 
  # using the date format  specified by 'date.fmt'
  if ( class(dates) == "factor" ) dates <- as.Date(dates, format= date.fmt)
  
  # If 'dates' is already of Date class, the following line verifies that
  # the number of days in 'dates' be equal to the number of element in the 
  # time series corresponding to the 'st.name' station
  if ( ( class(dates) == "Date") & (length(dates) != nrow(x) ) ) 
     stop("Invalid argument: 'length(dates)' must be equal to 'nrow(x)'")  
     
  # Amount of stations in 'x'
  nstations <- ncol(x)

  # ID of all the stations in 'x'
  snames <- colnames(x)  
  
  # Computing the Starting and Ending Year of the analysis
  Starting.Year <- as.numeric(format(range(dates)[1], "%Y"))
  Ending.Year   <- as.numeric(format(range(dates)[2], "%Y"))
  
  # Amount of Years belonging to the desired period
  nyears <- Ending.Year - Starting.Year + 1
  
  
  # Requiring the Zoo Library (Zoos ordered observations)
  require(zoo)
  
  if (verbose) print("Starting the computations...", quote=FALSE ) 
  

  # Creating the data.frame that will store the computed averages for each station
  z <- as.data.frame(matrix(data = NA, ncol = 12, nrow = nstations, 
                     byrow = TRUE, dimnames = NULL) )
  
  colnames(z) <- month.abb
  rownames(z) <- snames
  
  y = x
  
  for (j in 1:nstations) {
      
      if (verbose) print( paste("Station: ", format(snames[j], width=10, justify="left"),
                                " : ",format(j, width=3, justify="left"), "/", 
                                nstations, " => ", 
                                format(round(100*j/nstations,2), width=6, justify="left"), 
                                "%", sep=""), quote=FALSE )
                      
      # Transforming the column of 'x' into a zoo object, 
	    # using the dates provided by the user
	    tmp <- vector2zoo(x=y[,j], dates=dates, date.fmt=date.fmt)
	
	    # Computing the monthly values
	    z[j, ] <- as.numeric(dwdays.default(x=tmp, thr=thr, type=type, na.rm=na.rm, ...))
	    
                        
  } # FOR end
		
   
  return( z )
  
 } #'dwdays.data.frame' END
 
 
dwdays.matrix  <- function(x, thr=0, type="wet", na.rm=TRUE,
                            dates, 
                            date.fmt="%Y-%m-%d", 
							verbose=TRUE,...) {
                  
   x <- as.data.frame(x) 
      
   NextMethod("dwdays", x, thr=thr, type=type,  na.rm=na.rm,
                         dates=dates, 
                         date.fmt=date.fmt, 
                         verbose=verbose,...)
                  
} # 'dwdays.matrix  ' END
 
 
 
 
 
