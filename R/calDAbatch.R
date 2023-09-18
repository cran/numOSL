######### 
calDAbatch <- 
function(inputfile="inputDRtable", cfType="Liritzis2013", rdcf=0, rba=0, ShallowGamma=TRUE, 
         nsim=5000, rejectNeg=TRUE, outfile=paste(inputfile,"_Results",sep=""), 
         outpdf=paste(inputfile,"_Results",sep="")) {

  UseMethod("calDAbatch")

} # end function calDAbatch.
### 2023.09.17.
calDAbatch.default <- 
function(inputfile="inputDRtable", cfType="Liritzis2013", rdcf=0, rba=0, ShallowGamma=TRUE,
         nsim=5000, rejectNeg=TRUE, outfile=paste(inputfile,"_Results",sep=""), 
         outpdf=paste(inputfile,"_Results",sep="")) {

  stopifnot(length(inputfile)==1L, is.character(inputfile),
            length(cfType)==1L, cfType %in% c("Liritzis2013","Guerin2011","AdamiecAitken1998"),
            length(rdcf)==1L, is.numeric(rdcf), rdcf>=0.0,
            length(rba)==1L, is.numeric(rba), rba>=0.0,
            length(ShallowGamma)==1L, is.logical(ShallowGamma), 
            length(nsim)==1L, is.numeric(nsim), nsim>=10,
            length(rejectNeg)==1L, is.logical(rejectNeg), 
            is.null(outfile) || (length(outfile)==1L && is.character(outfile)),
            is.null(outpdf) || (length(outpdf)==1L && is.character(outpdf)))

  ###
  name1 <- c("SampleName", "De(Gy)", "Se.De(Gy)", "minGrainSize(um)", "maxGrainSize(um)",
             "Ucontent(ppm)", "Se.Ucontent(ppm)", "Thcontent(ppm)", "Se.Thcontent(ppm)",
             "Kcontent(percent)", "Se.Kcontent(percent)", "inKcontent.percent", "Se.inKcontent.percent",
             "Wct(percent)", "Se.Wct(percent)", "Depth(m)", "Se.Depth(m)", "Longitude(0)", "Se.Longitude(0)",
             "Latitude(0)", "Se.Latitude(0)", "Altitude(m)", "Se.Altitude(m)", "AlphaValue", "Se.AlphaValue", 
             "BulkDensity(g/cm3)", "Se.BulkDensity(g/cm3)")

  ###
  name2 <- c("Alpha.DR(Gy/ka)", "Se.Alpha.DR(Gy/ka)", "Beta.DR(Gy/ka)", "Se.Beta.DR(Gy/ka)", 
             "inBeta.DR(Gy/ka)", "Se.inBeta.DR(Gy/ka)", "Gamma.DR(Gy/ka)", "Se.Gamma.DR(Gy/ka)",  
             "Cosmic.DR(Gy/ka)", "Se.Cosmic.DR(Gy/ka)", "Total.DR(Gy/ka)", "Se.Total.DR(Gy/ka)", 
             "Age(ka)", "Se.Age(ka)")
 
  ###
  mycsv <- paste(inputfile,".csv",sep="")
  if (file.exists(mycsv)==FALSE) {
    TempTable <- data.frame(rbind(
    c("Coarse1",  2.21,  0.1105, 212, 250, 2.29, 0.1145,  4.07, 0.2035, 1.02, 0.051,  0,    0,     2.5,  1, 0.01, 0.005, -180,  0.01, -90,   0.01,   -5, 10, 0,     0,     2, 0.1),
    c("Coarse2",  15.01, 0.7505, 180, 212, 3.41, 0.1705,  4.64, 0.232,  1.26, 0.063,  0,    0,     8.1,  4, 0.05, 0.01,  -150,  0.01, -75,   0.01,    0, 10, 0,     0,     2, 0.1),
    c("Coarse3",  16.99, 0.8495, 150, 180, 3.87, 0.1935,  5.72, 0.286,  1.29, 0.0645, 0,    0,    11.9,  5, 0.1,  0.02,  -120,  0.01, -60,   0.01,    5, 10, 0,     0,     2, 0.1),
    c("Coarse4",  19.36, 0.968,  125, 150, 4.04, 0.202,   7.34, 0.367,  2.28, 0.114,  0,    0,    15.1,  5, 0.2,  0.05,   -90,  0.01, -45,   0.01,   25, 10, 0,     0,     2, 0.1),
    c("Coarse5",  23.08, 1.154,   90, 125, 4.2,  0.21,    8.98, 0.449,  2.64, 0.132,  0,    0,    19.8, 10, 0.3,  0.05,   -70,  0.01, -35,   0.01,   50, 10, 0,     0,     2, 0.1),
    c("Coarse6",  29.09, 1.4545,  63,  90, 4.3,  0.215,   9.91, 0.4955, 2.65, 0.1325, 0,    0,    33,   10, 0.4,  0.05,   -50,  0.01, -25,   0.01,  150, 10, 0,     0,     2, 0.1),
    c("Coarse7",  29.37, 1.4685,  38,  63, 5.34, 0.267,   9.93, 0.4965, 2.66, 0.133,  12.5, 0.5,  36.2, 10, 0.8,  0.05,   -30,  0.01, -15,   0.01,  300, 10, 0,     0,     2, 0.1),
    c("Coarse8",  32.36, 1.618,   63, 125, 5.61, 0.2805, 10.53, 0.5265, 3.32, 0.166,  12.5, 0.5,  40.6, 10, 1.2,  0.05,   -20,  0.01, -10,   0.01,  500, 10, 0,     0,     2, 0.1),
    c("Coarse9",  33.96, 1.698,   63, 150, 5.91, 0.2955, 10.98, 0.549,  3.45, 0.1725, 12.5, 0.5,  61.2, 20, 1.6,  0.05,   -10,  0.01,  -5,   0.01,  800, 10, 0,     0,     2, 0.1),
    c("Coarse10", 54.36, 2.718,   63, 250, 6.48, 0.324,  12.2,  0.61,   3.74, 0.187,  12.5, 0.5,  63.1, 20, 2,    0.05,    -5,  0.01,  -2.5, 0.01, 1100, 10, 0,     0,     2, 0.1),
    c("Coarse11", 60.99, 3.0495,  90, 180, 6.81, 0.3405, 15.82, 0.791,  3.87, 0.1935, 12.5, 0.5,  65.9, 20, 2.4,  0.05,     0,  0.01,   0,   0.01, 1500, 10, 0,     0,     2, 0.1),
    c("Coarse12", 71.03, 3.5515, 150, 250, 6.96, 0.348,  16.3,  0.815,  4.23, 0.2115, 12.5, 0.5,  69.7, 20, 2.8,  0.05,     5,  0.01,   2.5, 0.01, 2000, 10, 0,     0,     2, 0.1),
    c("Coarse13", 73.06, 3.653,  150, 180, 7.14, 0.357,  19.43, 0.9715, 5.22, 0.261,  12.5, 0.5,  69.4, 20, 3.2,  0.05,    10,  0.01,   5,   0.01, 2500, 10, 0,     0,     2, 0.1),
    c("Coarse14", 74.11, 3.7055,  90, 150, 7.17, 0.3585, 19.52, 0.976,  5.58, 0.279,  12.5, 0.5,  73.7, 40, 3.6,  0.05,    20,  0.01,  10,   0.01, 3300, 10, 0,     0,     2, 0.1),
    c("Coarse15", 77.67, 3.8835,  90, 212, 7.44, 0.372,  22.02, 1.101,  5.71, 0.2855, 12.5, 0.5,  77.5, 40, 4,    0.05,    40,  0.01,  20,   0.01, 3600, 10, 0,     0,     2, 0.1),
    c("Coarse16", 81.04, 4.052,   90, 250, 7.59, 0.3795, 23.28, 1.164,  5.73, 0.2865, 12.5, 0.5,  78.5, 40, 4.4,  0.05,    70,  0.01,  35,   0.01, 4000, 10, 0,     0,     2, 0.1),
    c("Coarse17", 85.01, 4.2505, 125, 180, 8.21, 0.4105, 23.94, 1.197,  5.83, 0.2915, 12.5, 0.5,  83.3, 40, 4.8,  0.05,    90,  0.01,  45,   0.01, 4200, 10, 0,     0,     2, 0.1),
    c("Coarse18", 89.73, 4.4865, 125, 212, 8.24, 0.412,  26.23, 1.3115, 6.07, 0.3035, 12.5, 0.5,  86,   40, 5.2,  0.05,   130,  0.01,  65,   0.01, 4500, 10, 0,     0,     2, 0.1),
    c("Coarse19", 91.26, 4.563,  150, 212, 8.26, 0.413,  31.55, 1.5775, 6.41, 0.3205, 12.5, 0.5,  93,   40, 5.6,  0.05,   150,  0.01,  75,   0.01, 4800, 10, 0,     0,     2, 0.1),
    c("Coarse20", 98.24, 4.912,  125, 250, 8.45, 0.4225, 36.51, 1.8255, 7.81, 0.3905, 12.5, 0.5, 100,   40, 6,    0.05,   180,  0.01,  90,   0.01, 4960, 10, 0,     0,     2, 0.1),
    c("Fine1",     2.21, 0.1105,   4,  11, 2.29, 0.1145,  4.07, 0.2035, 1.02, 0.051,  12.5, 0.5,   2.5,  1, 0.01, 0.005, -180,  0.01, -90,   0.01,   -5, 10, 0.038, 0.002, 2, 0.1),
    c("Fine2",    15.01, 0.7505,   4,  11, 3.41, 0.1705,  4.64, 0.232,  1.26, 0.063,  12.5, 0.5,   8.1,  4, 0.05, 0.01,  -150,  0.01, -75,   0.01,    0, 10, 0.038, 0.002, 2, 0.1),
    c("Fine3",    16.99, 0.8495,   4,  11, 3.87, 0.1935,  5.72, 0.286,  1.29, 0.0645, 12.5, 0.5,  11.9,  5, 0.1,  0.02,  -120,  0.01, -60,   0.01,    5, 10, 0.038, 0.002, 2, 0.1),
    c("Fine4",    19.36, 0.968,    4,  11, 4.04, 0.202,   7.34, 0.367,  2.28, 0.114,  12.5, 0.5,  15.1,  5, 0.2,  0.05,   -90,  0.01, -45,   0.01,   25, 10, 0.038, 0.002, 2, 0.1),
    c("Fine5",    23.08, 1.154,    4,  11, 4.2,  0.21,    8.98, 0.449,  2.64, 0.132,  12.5, 0.5,  19.8, 10, 0.3,  0.05,   -70,  0.01, -35,   0.01,   50, 10, 0.038, 0.002, 2, 0.1),
    c("Fine6",    29.09, 1.4545,   4,  11, 4.3,  0.215,   9.91, 0.4955, 2.65, 0.1325, 12.5, 0.5,  33,   10, 0.4,  0.05,   -50,  0.01, -25,   0.01,  150, 10, 0.038, 0.002, 2, 0.1),
    c("Fine7",    29.37, 1.4685,   4,  11, 5.34, 0.267,   9.93, 0.4965, 2.66, 0.133,  12.5, 0.5,  36.2, 10, 0.8,  0.05,   -30,  0.01, -15,   0.01,  300, 10, 0.038, 0.002, 2, 0.1),
    c("Fine8",    32.36, 1.618,    4,  11, 5.61, 0.2805, 10.53, 0.5265, 3.32, 0.166,  12.5, 0.5,  40.6, 10, 1.2,  0.05,   -20,  0.01, -10,   0.01,  500, 10, 0.038, 0.002, 2, 0.1),
    c("Fine9",    33.96, 1.698,    4,  11, 5.91, 0.2955, 10.98, 0.549,  3.45, 0.1725, 12.5, 0.5,  61.2, 20, 1.6,  0.05,   -10,  0.01,  -5,   0.01,  800, 10, 0.038, 0.002, 2, 0.1),
    c("Fine10",   54.36, 2.718,    4,  11, 6.48, 0.324,  12.2,  0.61,   3.74, 0.187,  12.5, 0.5,  63.1, 20, 2,    0.05,    -5,  0.01,  -2.5, 0.01, 1100, 10, 0.038, 0.002, 2, 0.1),
    c("Fine11",   60.99, 3.0495,   4,  11, 6.81, 0.3405, 15.82, 0.791,  3.87, 0.1935, 12.5, 0.5,  65.9, 20, 2.4,  0.05,     0,  0.01,   0,   0.01, 1500, 10, 0.038, 0.002, 2, 0.1),
    c("Fine12",   71.03, 3.5515,   4,  11, 6.96, 0.348,  16.3,  0.815,  4.23, 0.2115, 12.5, 0.5,  69.7, 20, 2.8,  0.05,     5,  0.01,   2.5, 0.01, 2000, 10, 0.038, 0.002, 2, 0.1),
    c("Fine13",   73.06, 3.653,    4,  11, 7.14, 0.357,  19.43, 0.9715, 5.22, 0.261,  12.5, 0.5,  69.4, 20, 3.2,  0.05,    10,  0.01,   5,   0.01, 2500, 10, 0.038, 0.002, 2, 0.1),
    c("Fine14",   74.11, 3.7055,   4,  11, 7.17, 0.3585, 19.52, 0.976,  5.58, 0.279,  12.5, 0.5,  73.7, 40, 3.6,  0.05,    20,  0.01,  10,   0.01, 3300, 10, 0.038, 0.002, 2, 0.1),
    c("Fine15",   77.67, 3.8835,   4,  11, 7.44, 0.372,  22.02, 1.101,  5.71, 0.2855, 12.5, 0.5,  77.5, 40, 4,    0.05,    40,  0.01,  20,   0.01, 3600, 10, 0.038, 0.002, 2, 0.1),
    c("Fine16",   81.04, 4.052,    4,  11, 7.59, 0.3795, 23.28, 1.164,  5.73, 0.2865, 12.5, 0.5,  78.5, 40, 4.4,  0.05,    70,  0.01,  35,   0.01, 4000, 10, 0.038, 0.002, 2, 0.1),
    c("Fine17",   85.01, 4.2505,   4,  11, 8.21, 0.4105, 23.94, 1.197,  5.83, 0.2915, 12.5, 0.5,  83.3, 40, 4.8,  0.05,    90,  0.01,  45,   0.01, 4200, 10, 0.038, 0.002, 2, 0.1),
    c("Fine18",   89.73, 4.4865,   4,  11, 8.24, 0.412,  26.23, 1.3115, 6.07, 0.3035, 12.5, 0.5,  86,   40, 5.2,  0.05,   130,  0.01,  65,   0.01, 4500, 10, 0.038, 0.002, 2, 0.1),
    c("Fine19",   91.26, 4.563,    4,  11, 8.26, 0.413,  31.55, 1.5775, 6.41, 0.3205, 12.5, 0.5,  93,   40, 5.6,  0.05,   150,  0.01,  75,   0.01, 4800, 10, 0.038, 0.002, 2, 0.1),
    c("Fine20",   98.24, 4.912,    4,  11, 8.45, 0.4225, 36.51, 1.8255, 7.81, 0.3905, 12.5, 0.5, 100,   40, 6,    0.05,   180,  0.01,  90,   0.01, 4960, 10, 0.038, 0.002, 2, 0.1)),
    stringsAsFactors=FALSE)

    ###
    colnames(TempTable) <- name1
  
    ###
    write.csv(TempTable, file=paste(inputfile,".csv",sep=""), row.names=FALSE)
    cat(paste("An editable dose-rate input template [", inputfile, ".csv] has been created at: ", paste(getwd(),"/", inputfile,".csv\n",sep=""),sep=""))
    
  } else {
    
    xoxoxo <- as.character(read.csv(file=mycsv, header=FALSE, stringsAsFactors=FALSE)[1L,])
    TempTable <- read.csv(file=mycsv, header=FALSE, skip=1L, stringsAsFactors=FALSE)
    
    ###
    idx1 <- which(!xoxoxo %in% name1)
    if (length(idx1)>=1L)  stop("Error: file [", inputfile, ".csv], invalid elements: ", paste(xoxoxo[idx1],collapse=", "))
    
    ###
    idx2 <- which(!name1 %in% xoxoxo)
    if (length(idx2)>=1L)  stop("Error: file [", inputfile, ".csv] cannot find elements: ", paste(name1[idx2],collapse=", "))
    
    ###
    N <- nrow(TempTable)
    drMAT <- matrix(nrow=N, ncol=14L)
    
    ###
    cat("Dose-rate and age calculations use [", inputfile, ".csv] (N=", N,") is in progress, please wait, ...\n",sep="")
    if (N>1L)  pb <- txtProgressBar(min=1, max=N, initial=1, style=3)

    ###
    if (!is.null(outpdf))  pdf(paste(outpdf, ".pdf", sep=""))
    
    ###
    LIST <- vector(length=N, mode="list")
    for (i in 1:N) {
      
      dose <- as.numeric(TempTable[i,c(2L,3L)])
      minGrainSize <- as.numeric(TempTable[i,4L])
      maxGrainSize <- as.numeric(TempTable[i,5L])
      Ucontent <- as.numeric(TempTable[i,c(6L,7L)])
      Thcontent <- as.numeric(TempTable[i,c(8L,9L)])
      Kcontent <- as.numeric(TempTable[i,c(10L,11L)])
      inKcontent <- as.numeric(TempTable[i,c(12L,13L)])
      Wct <- as.numeric(TempTable[i,c(14L,15L)])
      depth <- as.numeric(TempTable[i,c(16L,17L)])
      longitude <- as.numeric(TempTable[i,c(18L,19L)])
      latitude <- as.numeric(TempTable[i,c(20L,21L)])
      altitude <- as.numeric(TempTable[i,c(22L,23L)])
      alphaValue <- as.numeric(TempTable[i,c(24L,25L)])
      bulkDensity <- as.numeric(TempTable[i,c(26L,27L)])
      sampleName <- as.character(TempTable[i,1L])
     
      ###
      if (N>1L) setTxtProgressBar(pb,i) 
      ith_calDA <- try(calDA(dose=dose, minGrainSize=minGrainSize, maxGrainSize=maxGrainSize,
                       Ucontent=Ucontent, Thcontent=Thcontent, Kcontent=Kcontent, Wct=Wct, 
                       depth=depth, longitude=longitude, latitude=latitude, altitude=altitude, 
                       alphaValue=alphaValue, inKcontent=inKcontent, bulkDensity=bulkDensity, 
                       cfType=cfType, rdcf=rdcf, rba=rba, ShallowGamma=ShallowGamma, nsim=nsim,
                       rejectNeg=rejectNeg, plot=TRUE, sampleName=sampleName), silent=TRUE)
     
      ###
      if (inherits(ith_calDA, what="try-error")==TRUE) {

        cat("Failed in dose-rate and age calculations for the ", i, "-th sample!\n", sep="")
        print(paste("<",sampleName,"> ",attr(ith_calDA,"condition"),sep=""))
        LIST[[i]] <- NA

      } else {

        drMAT[i,] <- c(t(ith_calDA[,c(1L,2L),drop=FALSE]))
        LIST[[i]] <- ith_calDA

      } # end if.

    } # end for.

    ###
    if (N>1L) close(pb)
    if (!is.null(outpdf))  dev.off()
      
    ###
    TempTable[,28L:41L] <- drMAT
    colnames(TempTable) <- c(name1, name2)
    if (!is.null(outfile))  {

      write.csv(TempTable, file=paste(outfile,".csv",sep=""))

    } # end if.

    ###
    names(LIST) <- as.character(TempTable[,1L])
    return(invisible(LIST))
     
  } # end if.

} # end function calDAbatch.default.
###
