##### 
loadBINdata <- 
function(filename, view=TRUE)  {
    UseMethod("loadBINdata")
} #
### 2017.03.27.
loadBINdata.default <- 
function(filename, view=TRUE)  {
    ### Stop if not. 
    stopifnot(is.character(filename),
              is.logical(view)) 
    ###  
    TotalLumType <- c("TL", "OSL", "IRSL", "M-IR", "M-VIS", "TOL", "TRPOSL",
                      "RIR", "RBR", "USER", "POSL", "SGOSL", "RL", "XRF")
    TotalDataType <- c("Natural", "N+dose", "Bleach", "Bleach+dose", 
                       "Natural(Bleach)", "N+dose(Bleach)", "Dose", "Background")
    TotalLightSource <- c("None", "Lamp", "IRDiodes", "CalibraitionLED",
                          "BlueDiodes", "WhiteLight", "GreenLaser", "IRLaser")
    ###
    length_filename <- length(filename)
    ###
    all_records <- list()
    all_tab <- data.frame(stringsAsFactors=FALSE)
    ###
    ###
    for (i in seq(length_filename)) {
        binFile <- file(description=filename[i], open="rb")
        ###
        records <- list()
        iter <- 0L
        ###
        ### Read each record iteratively.
        repeat  {
            ###  Data format version number (Version), Byte (2); ----------- USE.
            Version <- readBin(binFile, what="integer", n=1L, size=2L, endian="little")
            if (length(Version)==0L)  break 
            ### 
            iter <- iter + 1L
            ###
            if ( Version %in% c(3L,4L) )  {             
                ### Length of this record (Length), Small Integer (2);
                ### Length of previous record (Previous), Small Integer (2);
                ### Number of data points (NPoints), Small Integer (2). ----------- USE.
                pass <- readBin(binFile, what="integer", n=3L, size=2L, endian="little") 
                NPoints <- pass[3L]
                ###

                ### Luminescence type (LType), Byte (1).  ----------- USE.
                LType <- readBin(binFile, what="integer", n=1L, size=1L, endian="little")
                ###

                ### Low (temperature, time, wavelength) (Low), Single (4); ----------- USE.
                ### High (temperature, time, wavelength) (High), Single (4); ----------- USE.
                ### Rate (heating rate, scan rate) (Rate), Single (4). ----------- USE.
                pass <- readBin(binFile, what="double", n=3L, size=4L, endian="little")
                Low <- pass[1L]
                High <- pass[2L]
                Rate <- pass[3L]
                ###

                ### Sample temperature (Temperature), Small Integer (2); ----------- USE.
                ### X position of a single grain (XCoord), Small Integer (2);
                ### Y position of a single grain (YCoord), Small Integer (2); 
                ### TOL "delay" channels (Delay), Small Integer (2); ----------- USE.
                ### TOL "on" channels (On), Small Integer (2); ----------- USE.
                ### TOL "off" channels (Off), Small Integer (2). ----------- USE.
                pass <- readBin(binFile, what="integer", n=6L, size=2L, endian="little")
                Temperature <- pass[1L]
                Delay <- pass[4L]
                On <- pass[5L]
                Off <- pass[6L]
                ###

                ### Carousel position (Position), Byte (1); ----------- USE.
                ### Run number (Run), Byte (1). ----------- USE.
                pass <- readBin(binFile, what="integer", n=2L, size=1L, endian="little")
                Position <- pass[1L]
                Run <- pass[2L]
                ###

                ### Data collection time (Time), String (7); ----------- USE.
                length_Time <- readBin(binFile, what="integer", n=1L, size=1L, endian="little")
                Time <- readChar(binFile, nchars=length_Time, useBytes=TRUE)
                if (length_Time<6L) pass <- readChar(binFile, nchars=6L-length_Time, useBytes=TRUE)

                ### Data collection date (Date), String (7). ----------- USE.
                length_Date <- readBin(binFile, what="integer", n=1L, size=1L, endian="little")
                Date <- readChar(binFile, nchars=length_Date, useBytes=TRUE)
                if (length_Date<6L) pass <- readChar(binFile, nchars=6L-length_Date, useBytes=TRUE)
                ###

                ### Sequence name (Sequence), String (9);
                ### User name (User), String (9).
                pass <- readChar(binFile, nchars=18L, useBytes=TRUE)
                ###

                ### Data type (DType), Byte (1). ----------- USE.
                DType <- readBin(binFile, what="integer", n=1L, size=1L, endian="little")
                ###

                ### Irradiation time (IRR_Time), Single (4). ----------- USE.
                IRRTime <- readBin(binFile, what="double", n=1L, size=4L, endian="little") 
                ###

                ### Irradiation type (alpha, beta or gamma) (IRR_Type), Byte (1);
                ### Irradiation unit (Gy, Rads, secs, mins, hrs) (IRR_Unit), Byte (1). 
                pass <- readBin(binFile, what="integer", n=2L, size=1L, endian="little")
                ###

                ### Bleaching time (BL_Time), Single (4).    
                pass <- readBin(binFile, what="double", n=1L, size=4L, endian="little")
                ###

                ### Bleaching unit (mJ, J, secs, mins, hrs) (BL_Unit), Byte (1).    
                pass <- readBin(binFile, what="integer", n=1L, size=1L, endian="little")
                ###

                ### Annealing temperature (AN_Temp), Single (4); ----------- USE.
                ### Annealing time (AN_Time), Single(4);
                ### Normalisation factor (1) (Norm1), Single (4);
                ### Normalisation factor (2) (Norm2), Single (4);
                ### Normalisation factor (3) (Norm3), Single (4);
                ### Background level (BG), Single (4).
                pass <- readBin(binFile, what="double", n=6L, size=4L, endian="little")
                AnTemp <- pass[1L]
                TimeSinceIrr <- -1L
                ###

                ### Number of channels to shift data (Shift), Small Integer (2).    
                pass <- readBin(binFile, what="integer", n=1L, size=2L, endian="little")
                ###

                ### Sample name (Sample), String (21).
                pass <- readChar(binFile, nchars=21L, useBytes=TRUE) 
                ###

                ### Comment (Comment), String (81).
                pass <- readChar(binFile, nchars=81L, useBytes=TRUE) 
                ###

                ### Light Source (LightSource), Byte (1); ----------- USE.
                ### Set Number (SET), Byte (1); ----------- USE.
                ### Tag (TAG), Byte(1).
                pass <- readBin(binFile, what="integer", n=3L, size=1L, endian="little") 
                LightSource <- pass[1L]
                Set <- pass[2L]
                ###

                ### Grain number (Grain), Small Integer (2). ----------- USE. 
                Grain <- readBin(binFile, what="integer", n=1L, size=2L, endian="little")      
                ###

                ### Optical stimulation power (LightPower), Single (4).    
                pass <- readBin(binFile, what="double", n=1L, size=4L, endian="little")  
                ###

                ### System ID (ystemID), Small Integer (2).
                pass <- readBin(binFile, what="integer", n=1L, size=2L, endian="little")
                ###

                ### Reserved (Reserved), Byte (54).    
                pass <- readBin(binFile, what="integer", n=54L, size=1L, endian="little")
                ###

                ### Data array of nPoints long integers (DPoints), Integer (4). ----------- USE.
                DPoints <- readBin(binFile, what="integer", n=NPoints, size=4L, endian="little")
                ###
            } else if (Version %in% c(6L, 7L, 8L)) {
                ###
                ### Length of this record (Length), Long Integer (4);
                ### Length of previous record (Previous), Long Integer (4);
                ### Number of data points (NPoints), Long Integer (4). ----------- USE.
                pass <- readBin(binFile, what="integer", n=3L, size=4L, endian="little") 
                NPoints <- pass[3L]

                ### Version 8 only, Record Type, Byte(1).
                if (Version==8L) {
                    pass <- readBin(binFile, what="integer", n=1L, size=1L, endian="little")
                } # end if.

                ###
                ### Run number (Run), Small Integer (2); ----------- USE.
                ### Set number (Set), Small Integer (2); ----------- USE.
                ### Carousel position (Position), Small Integer (2); ----------- USE.
                ### Grain number (GrainNumber), Small Integer (2); ----------- USE.
                ### Curve number (for multiple curve operations) (CurveNo), Small Integer (2);
                ###  X position of a single grain (Xcoord), Small Integer (2);
                ###  Y position of a single grain (Ycoord), Small Integer (2). 
                pass <- readBin(binFile, what="integer", n=7L, size=2L, endian="little")     
                Run <- pass[1L]
                Set <- pass[2L]
                Position <- pass[3L] 
                Grain <- pass[4L]
                ###

                ### Sample name (Sample), String (21).
                pass <- readChar(binFile, nchars=21L, useBytes=TRUE) 
                ###

                ### Comment (comment), String (81).
                pass <- readChar(binFile, nchars=81L, useBytes=TRUE) 
                ###

                ### System ID (SystemID), Small Integer (2). 
                pass <- readBin(binFile, what="integer", n=1L, size=2L, endian="little") 
                ###

                ### File name (.SEC, BINX ect) (FName), String (101).
                pass <- readChar(binFile, nchars=101L, useBytes=TRUE)
                ###

                ### User name (User), String (31).
                pass <- readChar(binFile, nchars=31L, useBytes=TRUE)
                ###

                ### Data collection time (hh-mm-ss) (Time), String (7); ----------- USE.
                length_Time <- readBin(binFile, what="integer", n=1L, size=1L, endian="little")
                Time <- readChar(binFile, nchars=length_Time, useBytes=TRUE)
                if (length_Time<6L) pass <- readChar(binFile, nchars=6L-length_Time, useBytes=TRUE)

                ### Data collection date (dd-mm-yy) (Date), String (7). ----------- USE.
                length_Date <- readBin(binFile, what="integer", n=1L, size=1L, endian="little")
                Date <- readChar(binFile, nchars=length_Date, useBytes=TRUE)
                if (length_Date<6L) pass <- readChar(binFile, nchars=6L-length_Date, useBytes=TRUE)
                ###

                ### Data type (DType), Byte (1). ----------- USE.
                DType <- readBin(binFile, what="integer", n=1L, size=1L, endian="little")
                ###

                ### Bleaching time (BL_Time), Single (4).     
                pass <- readBin(binFile, what="double", n=1L, size=4L, endian="little")
                ###

                ### Bleaching unit (mJ, J, secs, mins, hrs) (BL_Unit), Byte (1).    
                pass <- readBin(binFile, what="integer", n=1L, size=1L, endian="little")
                ###

                ### Normalisation factor (1) (Norm1), Single (4);
                ### Normalisation factor (2) (Norm2), Single (4);
                ### Normalisation factor (3) (Norm3), Single (4);
                ### Background level (BG), single (4).
                pass <- readBin(binFile, what="double", n=4L, size=4L, endian="little") 
                ### 

                ### Number of channels to shift data (Shift), Small Integer (2).    
                pass <- readBin(binFile, what="integer", n=1L, size=2L, endian="little")
                ###

                ### Tag (Tag), Byte (1).
                pass <- readBin(binFile, what="integer", n=1L, size=1L, endian="little")
                ###   

                ### Reserved for internal use (Reserved), (20).
                pass <- readBin(binFile, what="raw", n=20L, size=1L, endian="little")
                ###

                ### Luminescence type (LType), Byte (1); ----------- USE.
                ### Light Source (LightSource), Byte (1).  ----------- USE.
                pass <- readBin(binFile, what="integer", n=2L, size=1L, endian="little")
                LType <- pass[1L]
                LightSource <- pass[2L]
                ###

                ### Optical Stimulation Power (LightPower), Single (4); 
                ### Low (temperature, time, wavelength) (Low), Single (4); ----------- USE.
                ### High (temperature, time, wavelength) (High), Single (4); ----------- USE.
                ### Rate (heating rate, scan rate) (Rate), Single (4). ----------- USE.
                pass <- readBin(binFile, what="double", n=4L, size=4L, endian="little")
                Low <- pass[2L]
                High <- pass[3L]
                Rate <- pass[4L]
                ###

                ### Sample temperature (Temperature), Small Integer (2); ----------- USE.
                ### Measurement temperature (MeasTemp),Small Integer (2).
                pass <- readBin(binFile, what="integer", n=2L, size=2L, endian="little")
                Temperature <- pass[1L]
                ###

                ### Preheating temperature (An_Temp), Single (4); ----------- USE.
                ### Preheating time (An_Time), Single (4).
                pass <- readBin(binFile, what="double", n=2L, size=4L, endian="little")
                AnTemp <- pass[1L]
                ###

                ### TOL "delay" channels (Delay), Small Integer (2); ----------- USE.
                ### TOL "on" channels (On), Small Integer (2); ----------- USE.
                ### TOL "off" channels (off), Small Integer (2).----------- USE.
                pass <- readBin(binFile, what="integer", n=3L, size=2L, endian="little")
                Delay <- pass[1L]
                On <- pass[2L]
                Off <- pass[3L]
                ###

                ### Irradiation time (IRR_Time), Single (4). ----------- USE.
                IRRTime <- readBin(binFile, what="double", n=1L, size=4L, endian="little")
                ###

                ### Irradiation type (alpha, beta or gamma) (IRR_Type), Byte (1).
                pass <- readBin(binFile, what="integer", n=1L, size=1L, endian="little")
                ### 

                ### Irradiation dose rate (Gy/s) (IRR_DoseRate), Single (4);
                ### Irradiation dose rate error (Gy/s) (DoseRateErr), Single (4).
                pass <- readBin(binFile, what="double", n=2L, size=4L, endian="little")
                ###

                ### Time since last irradiation (s) (TimeSinceIrr), Long Integer (4); ----------- USE.
                TimeSinceIrr <- readBin(binFile, what="integer", n=1L, size=4L, endian="little")
                ### 
    
                ### Time unit (time tick) for pulse parameters (s) (TimeTick), Single (4).
                pass <- readBin(binFile, what="double", n=1L, size=4L, endian="little")
                ###

                ### On-time for pulsed stimulation (in time ticks) (OnTime), Long Integer (4);
                ### Stimulation period (on+off time in time ticks) (StimPeriod), Long Integer (4).
                pass <- readBin(binFile, what="integer", n=2L, size=4L, endian="little")
                ###  
                ### PMT signal gating enabled (GateEnabled), Byte (1).
                pass <- readBin(binFile, what="integer", n=1L, size=1L, endian="little")
                ###

                ### Start of gating (in time ticks from start of on pulse) (GateStart), Long Integer (4);
                ### End of gating (in time ticks from start of on pulse) (GateEnd), Long Integer (4).
                pass <- readBin(binFile, what="integer", n=2L, size=4L, endian="little")
                ###

                ### Photon Timer enabled (PTenabled), Byte (1);
                ### PMT dead time correction enabled (DTenabled), Byte (1).
                pass <- readBin(binFile, what="integer", n=2L, size=1L, endian="little")
                ### 

                ### PMT dead time (s) (DeadTime), Single (4);
                ### Stimulation power corresponding to 100% (mw/cm2) (MaxLPower), Single (4);
                ### XRF acquisition time (s) (XrfAcqTime), Single (4);
                ### XRF X-ray high voltage (V) (XrfHV), Single (4). 
                pass <- readBin(binFile, what="double", n=4L, size=4L, endian="little")
                ###

                ### XRF X-ray current (uA) (XrfCurr), Long Integer (4).
                pass <- readBin(binFile, what="integer", n=1L, size=4L, endian="little")
                ###

                ### XRF dead time fraction (XrfDeadTimeF), Single (4).
                pass <- readBin(binFile, what="double", n=1L, size=4L, endian="little")
                ###

                ### Reserved for internal use (Reserved),   Byte (24 or 33). 
                if (Version==6L)  {
                    pass <- readBin(binFile, what="raw", n=24L, size=1L, endian="little")
                } else if (Version==7L)  {
                    pass <- readBin(binFile, what="raw", n=24L, size=1L, endian="little")
                } else if (Version==8L)  {
                    pass <- readBin(binFile, what="raw", n=83L, size=1L, endian="little")
                } # end if.
                ###

                ###
                ### Data array of nPoints Long Integers (DPoints), Long Integer (4).----------- USE.
                DPoints <- readBin(binFile, what="integer", n=NPoints, size=4L, endian="little")
                ###
            } else {
                stop("Error: invalid version number!")
            } # end if.
            ###
            ###
            attr(DPoints, "NPoints") <- NPoints
            attr(DPoints, "LType") <- TotalLumType[LType+1L]
            attr(DPoints, "Low") <- Low
            attr(DPoints, "High") <- High
            attr(DPoints, "Rate") <- Rate
            attr(DPoints, "Temperature") <- Temperature
            attr(DPoints, "Delay") <- Delay
            attr(DPoints, "On") <- On
            attr(DPoints, "Off") <- Off
            attr(DPoints, "Position") <- Position 
            attr(DPoints, "Run") <- Run
            attr(DPoints, "Time") <- Time
            attr(DPoints, "Date") <- Date
            attr(DPoints, "DType") <- TotalDataType[DType+1L] 
            attr(DPoints, "IRRTime") <- round(IRRTime, digits=3L) 
            attr(DPoints, "AnTemp") <- AnTemp
            attr(DPoints, "LightSource") <- TotalLightSource[LightSource+1L]
            attr(DPoints, "Set") <- Set
            attr(DPoints, "Grain") <- Grain
            attr(DPoints, "TimeSinceIrr") <- TimeSinceIrr
            ###
            ###
            records[[iter]] <- DPoints
        } # end repeat.
        ### 
        close(binFile)
        ###
        tab <- data.frame("Position"=sapply(X=records, FUN=attr, "Position"),
                          "Grain"=sapply(X=records, FUN=attr, "Grain"),
                          "Run"=sapply(X=records, FUN=attr, "Run"),
                          "Set"=sapply(X=records, FUN=attr, "Set"),
                          "DType"=sapply(X=records, FUN=attr, "DType"),
                          "IRRTime"=sapply(X=records, FUN=attr, "IRRTime"),
                          "NPoints"=sapply(X=records, FUN=attr, "NPoints"),
                          "LType"=sapply(X=records, FUN=attr, "LType"),
                          "Low"=sapply(X=records, FUN=attr, "Low"),
                          "High"=sapply(X=records, FUN=attr, "High"),
                          "Rate"=sapply(X=records, FUN=attr, "Rate"),
                          "Temperature"=sapply(X=records, FUN=attr, "Temperature"),
                          "Delay"=sapply(X=records, FUN=attr, "Delay"),
                          "On"=sapply(X=records, FUN=attr, "On"),
                          "Off"=sapply(X=records, FUN=attr, "Off"),
                          "LightSource"=sapply(X=records, FUN=attr, "LightSource"),
                          "AnTemp"=sapply(X=records, FUN=attr, "AnTemp"), 
                          "TimeSinceIrr"=sapply(X=records, FUN=attr, "TimeSinceIrr"),              
                          "Time"=sapply(X=records, FUN=attr, "Time"),
                          "Date"=sapply(X=records, FUN=attr, "Date"),
                          stringsAsFactors=FALSE)
        ###
        cat(paste("File [",i,"]: successfully loaded ", length(records), " records!\n\n", sep=""))
        ###
        all_records <- c(all_records, records)
        ###
        all_tab <- rbind(all_tab, tab)
    } # end for.
    ###
    if (view==TRUE) View(x=all_tab, title="Loaded Table of Records")
    ###
    output <- list("records"=all_records, "tab"=all_tab)
    class(output) <- "loadBIN"
    ###
    invisible(output)
} # end function loadBINdata.
#####
