######### 
calDAbatch <- 
function(inputfile="inputDRtable", cfType="Liritzis2013", rdcf=0, rba=0, calRbfromK=FALSE, 
         ShallowGamma=TRUE, nsim=5000, reject=TRUE, outfile=paste(inputfile,"_Results",sep=""), 
         outpdf=paste(inputfile,"_Results",sep=""), digits=4) {

  UseMethod("calDAbatch")

} # end function calDAbatch.
### 2023.12.06.
calDAbatch.default <- 
function(inputfile="inputDRtable", cfType="Liritzis2013", rdcf=0, rba=0, calRbfromK=FALSE, 
         ShallowGamma=TRUE, nsim=5000, reject=TRUE, outfile=paste(inputfile,"_Results",sep=""), 
         outpdf=paste(inputfile,"_Results",sep=""), digits=4) {

  stopifnot(length(inputfile)==1L, is.character(inputfile),
            length(cfType)==1L, cfType %in% c("Liritzis2013","Guerin2011","AdamiecAitken1998"),
            length(rdcf)==1L, is.numeric(rdcf), rdcf>=0.0,
            length(rba)==1L, is.numeric(rba), rba>=0.0,
            length(calRbfromK)==1L, is.logical(calRbfromK),
            length(ShallowGamma)==1L, is.logical(ShallowGamma), 
            length(nsim)==1L, is.numeric(nsim), nsim>=10,
            length(reject)==1L, is.logical(reject), 
            is.null(outfile) || (length(outfile)==1L && is.character(outfile)),
            is.null(outpdf) || (length(outpdf)==1L && is.character(outpdf)),
            length(digits)==1L, is.numeric(digits), digits>=0)

  ###
  name1 <- c("SampleName", "De(Gy)", "Se.De(Gy)", "minGrainSize(um)", "maxGrainSize(um)",
             "Ucontent(ppm)", "Se.Ucontent(ppm)", "Thcontent(ppm)", "Se.Thcontent(ppm)",
             "Kcontent(percent)", "Se.Kcontent(percent)", "Rbcontent(ppm)", "Se.Rbcontent(ppm)", 
             "inKcontent(percent)", "Se.inKcontent(percent)", "inRbcontent(ppm)", "Se.inRbcontent(ppm)", 
             "Wct(percent)", "Se.Wct(percent)", "Depth(m)", "Se.Depth(m)", "Longitude(degree)", "Se.Longitude(degree)",
             "Latitude(degree)", "Se.Latitude(degree)", "Altitude(m a.s.l.)", "Se.Altitude(m a.s.l.)", "AlphaValue",  
             "Se.AlphaValue", "BulkDensity(g/cm3)", "Se.BulkDensity(g/cm3)")

  ###
  name2 <- c("Alpha.DR(Gy/ka)", "Se.Alpha.DR(Gy/ka)", "Beta.DR(Gy/ka)", "Se.Beta.DR(Gy/ka)", 
             "inBeta.DR(Gy/ka)", "Se.inBeta.DR(Gy/ka)", "Gamma.DR(Gy/ka)", "Se.Gamma.DR(Gy/ka)",  
             "Cosmic.DR(Gy/ka)", "Se.Cosmic.DR(Gy/ka)", "Total.DR(Gy/ka)", "Se.Total.DR(Gy/ka)", 
             "Age(ka)", "Se.Age(ka)")
 
  ###
  mycsv <- paste(inputfile,".csv",sep="")
  if (file.exists(mycsv)==FALSE) {
    TempTable <- data.frame(rbind(
    c("CoarseQuartz1",    2.21,0.1105,212,250,2.29,0.1145, 4.07,0.2035,1.02,0.051,   0,      0,        0,   0,    0,  0,  2.5, 1,0.01,0.005,-179,0.01,-89,  0.01,  -5,10,0,    0,    2,  0.1),
    c("CoarseQuartz2",   15.01,0.7505,180,212,3.41,0.1705, 4.64,0.232, 1.26,0.063,   0,      0,        0,   0,    0,  0,  8.1, 4,0.05,0.01, -150,0.01,-75,  0.01,   0,10,0,    0,    2,  0.1),
    c("CoarseQuartz3",   16.99,0.8495,150,180,3.87,0.1935, 5.72,0.286, 1.29,0.0645,  0,      0,        0,   0,    0,  0, 11.9, 5,0.1, 0.02, -120,0.01,-60,  0.01,   5,10,0,    0,    2,  0.1),
    c("CoarseQuartz4",   19.36,0.968, 125,150,4.04,0.202,  7.34,0.367, 2.28,0.114,   0,      0,        0,   0,    0,  0, 15.1, 5,0.2, 0.05,  -89,0.01,-45,  0.01,  25,10,0,    0,    2,  0.1),
    c("CoarseQuartz5",   23.08,1.154,  90,125,4.2, 0.21,   8.98,0.449, 2.64,0.132,   0,      0,        0,   0,    0,  0, 19.8,10,0.3, 0.05,  -70,0.01,-35,  0.01,  50,10,0,    0,    2.5,0.1),
    c("CoarseQuartz6",   29.09,1.4545, 63, 90,4.3, 0.215,  9.91,0.4955,2.65,0.1325,  0,      0,        0,   0,    0,  0, 33,  10,0.4, 0.05,  -50,0.01,-25,  0.01, 150,10,0,    0,    2.5,0.1),
    c("CoarseQuartz7",   32.36,1.618,  63,125,5.61,0.2805,10.53,0.5265,3.32,0.166, 117.4216, 5.87108,  0,   0,    0,  0, 40.6,10,1.2, 0.05,  -20,0.01,-10,  0.01, 500,10,0,    0,    2.5,0.1),
    c("CoarseQuartz8",   33.96,1.698,  63,150,5.91,0.2955,10.98,0.549, 3.45,0.1725,122.3785, 6.11892,  0,   0,    0,  0, 61.2,20,1.6, 0.05,  -10,0.01, -5,  0.01, 800,10,0,    0,    2.5,0.1),
    c("CoarseQuartz9",   54.36,2.718,  63,250,6.48,0.324, 12.2, 0.61,  3.74,0.187, 133.4362, 6.67181,  0,   0,    0,  0, 63.1,20,2,   0.05,   -5,0.01, -2.5,0.01,1100,10,0,    0,    2.5,0.1),
    c("CoarseKfeldspar1",60.99,3.0495, 90,180,6.81,0.3405,15.82,0.791, 3.87,0.1935,138.3931, 6.919655,12.5, 0.5,400,100, 65.9,20,2.4, 0.05,    0,0.01,  0,  0.01,1500,10,0,    0,    2,  0.1),
    c("CoarseKfeldspar2",71.03,3.5515,150,250,6.96,0.348, 16.3, 0.815, 4.23,0.2115,152.1199, 7.605995,12.5, 0.5,400,100, 69.7,20,2.8, 0.05,    5,0.01,  2.5,0.01,2000,10,0,    0,    2,  0.1),
    c("CoarseKfeldspar3",73.06,3.653, 150,180,7.14,0.357, 19.43,0.9715,5.22,0.261, 189.8686, 9.49343, 12.5, 0.5,400,100, 69.4,20,3.2, 0.05,   10,0.01,  5,  0.01,2500,10,0,    0,    2,  0.1),
    c("CoarseKfeldspar4",74.11,3.7055, 90,150,7.17,0.3585,19.52,0.976, 5.58,0.279, 203.5954,10.17977, 12.5, 0.5,400,100, 73.7,40,3.6, 0.05,   20,0.01, 10,  0.01,3300,10,0,    0,    2,  0.1),
    c("CoarseKfeldspar5",77.67,3.8835, 90,212,7.44,0.372, 22.02,1.101, 5.71,0.2855,208.5523,10.427615,13,   1,  400,100, 77.5,40,4,   0.05,   40,0.01, 20,  0.01,3600,10,0,    0,    2.5,0.1),
    c("CoarseKfeldspar6",81.04,4.052,  90,250,7.59,0.3795,23.28,1.164, 5.73,0.2865,209.3149,10.465745,13,   1,  400,100, 78.5,40,4.4, 0.05,   70,0.01, 35,  0.01,4000,10,0,    0,    2.5,0.1),
    c("CoarseKfeldspar7",85.01,4.2505,125,180,8.21,0.4105,23.94,1.197, 5.83,0.2915,213.1279,10.656395,13,   1,  400,100, 83.3,40,4.8, 0.05,   89,0.01, 45,  0.01,4200,10,0,    0,    2.5,0.1),
    c("CoarseKfeldspar8",89.73,4.4865,125,212,8.24,0.412, 26.23,1.3115,6.07,0.3035,222.2791,11.113955,10,   2,  400,100, 86,  40,5.2, 0.05,  130,0.01, 65,  0.01,4500,10,0,    0,    2.5,0.1),
    c("CoarseKfeldspar9",91.26,4.563, 150,212,8.26,0.413, 31.55,1.5775,6.41,0.3205,235.2433,11.762165,10,   2,  400,100, 93,  40,5.6, 0.05,  150,0.01, 75,  0.01,4800,10,0,    0,    2.5,0.1),
    c("MediumQuartz1",   98.24,4.912,  38, 63,8.45,0.4225,36.51,1.8255,7.81,0.3905,288.6253,14.431265, 0,   0,    0,  0,100,  40,6,   0.05,  179,0.01, 89,  0.01,4960,10,0.035,0.03, 2,  0.1),
    c("MediumQuartz2",   29.37,1.4685, 38, 63,5.34,0.267,  9.93,0.4965,2.66,0.133,  92.2558, 4.61279,  0,   0,    0,  0, 36.2,10,0.8, 0.05,  -30,0.01,-15,  0.01, 300,10,0.035,0.03, 2,  0.1),
    c("MediumKfeldspar1",98.24,4.912,  38, 63,8.45,0.4225,36.51,1.8255,7.81,0.3905,288.6253,14.431265,12.5, 0.5,400,100,100,  40,6,   0.05,  179,0.01, 89,  0.01,4960,10,0.1,  0.02, 2,  0.1),
    c("MediumKfeldspar2", 2.21,0.1105, 38, 63,2.29,0.1145, 4.07,0.2035,1.02,0.051,  92.2558, 4.61279, 12.5, 0.5,400,100,  2.5, 1,0.01,0.005,-179,0.01,-89,  0.01,  -5,10,0.1,  0.02, 2,  0.1),
    c("FineQuartz1",     15.01,0.7505,  4, 11,3.41,0.1705, 4.64,0.232, 1.26,0.063,   0,      0,        0,   0,    0,  0,  8.1, 4,0.05,0.01, -150,0.01,-75,  0.01,   0,10,0.04, 0.01, 2,  0.1),
    c("FineQuartz2",     16.99,0.8495,  4, 11,3.87,0.1935, 5.72,0.286, 1.29,0.0645,  0,      0,        0,   0,    0,  0, 11.9, 5,0.1, 0.02, -120,0.01,-60,  0.01,   5,10,0.04, 0.01, 2,  0.1),
    c("FineQuartz3",     19.36,0.968,   4, 11,4.04,0.202,  7.34,0.367, 2.28,0.114,   0,      0,        0,   0,    0,  0, 15.1, 5,0.2, 0.05,  -89,0.01,-45,  0.01,  25,10,0.04, 0.01, 2,  0.1),
    c("FineQuartz4",     23.08,1.154,   4, 11,4.2, 0.21,   8.98,0.449, 2.64,0.132,   0,      0,        0,   0,    0,  0, 19.8,10,0.3, 0.05,  -70,0.01,-35,  0.01,  50,10,0.04, 0.01, 2,  0.1),
    c("FineQuartz5",     29.09,1.4545,  4, 11,4.3, 0.215,  9.91,0.4955,2.65,0.1325,  0,      0,        0,   0,    0,  0, 33,  10,0.4, 0.05,  -50,0.01,-25,  0.01, 150,10,0.04, 0.01, 2.5,0.1),
    c("FineQuartz6",     29.37,1.4685,  4, 11,5.34,0.267,  9.93,0.4965,2.66,0.133,   0,      0,        0,   0,    0,  0, 36.2,10,0.8, 0.05,  -30,0.01,-15,  0.01, 300,10,0.038,0.002,2.5,0.1),
    c("FineQuartz7",     32.36,1.618,   4, 11,5.61,0.2805,10.53,0.5265,3.32,0.166, 117.4216, 5.87108,  0,   0,    0,  0, 40.6,10,1.2, 0.05,  -20,0.01,-10,  0.01, 500,10,0.038,0.002,2.5,0.1),
    c("FineQuartz8",     33.96,1.698,   4, 11,5.91,0.2955,10.98,0.549, 3.45,0.1725,122.3785, 6.118925, 0,   0,    0,  0, 61.2,20,1.6, 0.05,  -10,0.01, -5,  0.01, 800,10,0.038,0.002,2.5,0.1),
    c("FineQuartz9",     54.36,2.718,   4, 11,6.48,0.324, 12.2, 0.61,  3.74,0.187, 133.4362, 6.67181,  0,   0,    0,  0, 63.1,20,2,   0.05,   -5,0.01, -2.5,0.01,1100,10,0.038,0.002,2.5,0.1),
    c("FinePolymineral1",60.99,3.0495,  4, 11,6.81,0.3405,15.82,0.791, 3.87,0.1935,138.3931, 6.919655,13,   1,  400,100, 65.9,20,2.4, 0.05,    0,0.01,  0,  0.01,1500,10,0.086,0.004,2,  0.1),
    c("FinePolymineral2",71.03,3.5515,  4, 11,6.96,0.348, 16.3, 0.815, 4.23,0.2115,152.1199, 7.605995,13,   1,  400,100, 69.7,20,2.8, 0.05,    5,0.01,  2.5,0.01,2000,10,0.086,0.004,2,  0.1),
    c("FinePolymineral3",73.06,3.653,   4, 11,7.14,0.357, 19.43,0.9715,5.22,0.261, 189.8686, 9.49343, 13,   1,  400,100, 69.4,20,3.2, 0.05,   10,0.01,  5,  0.01,2500,10,0.086,0.004,2,  0.1),
    c("FinePolymineral4",74.11,3.7055,  4, 11,7.17,0.3585,19.52,0.976, 5.58,0.279, 203.5954,10.17977, 13,   1,  400,100, 73.7,40,3.6, 0.05,   20,0.01, 10,  0.01,3300,10,0.086,0.004,2,  0.1),
    c("FinePolymineral5",77.67,3.8835,  4, 11,7.44,0.372, 22.02,1.101, 5.71,0.2855,208.5523,10.427615,12.5, 0.5,400,100, 77.5,40,4,   0.05,   40,0.01, 20,  0.01,3600,10,0.086,0.004,2.5,0.1),
    c("FinePolymineral6",81.04,4.052,   4, 11,7.59,0.3795,23.28,1.164, 5.73,0.2865,209.3149,10.465745,12.5, 0.5,400,100, 78.5,40,4.4, 0.05,   70,0.01, 35,  0.01,4000,10,0.09, 0.02, 2.5,0.1),
    c("FinePolymineral7",85.01,4.2505,  4, 11,8.21,0.4105,23.94,1.197, 5.83,0.2915,213.1279,10.656395,12.5, 0.5,400,100, 83.3,40,4.8, 0.05,   89,0.01, 45,  0.01,4200,10,0.09, 0.02, 2.5,0.1),
    c("FinePolymineral8",89.73,4.4865,  4, 11,8.24,0.412, 26.23,1.3115,6.07,0.3035,222.2791,11.113955,10,   2,  400,100, 86,  40,5.2, 0.05,  130,0.01, 65,  0.01,4500,10,0.09, 0.02, 2.5,0.1),
    c("FinePolymineral9",91.26,4.563,   4, 11,8.26,0.413, 31.55,1.5775,6.41,0.3205,235.2433,11.762165,10,   2,  400,100, 93,  40,5.6, 0.05,  150,0.01, 75,  0.01,4800,10,0.09, 0.02, 2.5,0.1)),
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
      Rbcontent <- as.numeric(TempTable[i,c(12L,13L)])
      inKcontent <- as.numeric(TempTable[i,c(14L,15L)])
      inRbcontent <- as.numeric(TempTable[i,c(16L,17L)])
      Wct <- as.numeric(TempTable[i,c(18L,19L)])
      depth <- as.numeric(TempTable[i,c(20L,21L)])
      longitude <- as.numeric(TempTable[i,c(22L,23L)])
      latitude <- as.numeric(TempTable[i,c(24L,25L)])
      altitude <- as.numeric(TempTable[i,c(26L,27L)])
      alphaValue <- as.numeric(TempTable[i,c(28L,29L)])
      bulkDensity <- as.numeric(TempTable[i,c(30L,31L)])
      sampleName <- as.character(TempTable[i,1L])
     
      ###
      if (N>1L) setTxtProgressBar(pb,i) 
      ith_calDA <- try(calDA(dose=dose, minGrainSize=minGrainSize, maxGrainSize=maxGrainSize,
                       Ucontent=Ucontent, Thcontent=Thcontent, Kcontent=Kcontent, Rbcontent=Rbcontent, 
                       Wct=Wct, depth=depth, longitude=longitude, latitude=latitude, altitude=altitude, 
                       alphaValue=alphaValue, inKcontent=inKcontent, inRbcontent=inRbcontent, 
                       calRbfromK=calRbfromK, bulkDensity=bulkDensity, cfType=cfType, rdcf=rdcf, 
                       rba=rba, ShallowGamma=ShallowGamma, nsim=nsim, reject=reject, plot=TRUE, 
                       sampleName=sampleName), silent=TRUE)
     
      ###
      if (inherits(ith_calDA, what="try-error")==TRUE) {

        cat("Failed in dose-rate and age calculations for the ", i, "-th sample!\n", sep="")
        print(paste("<",sampleName,"> ",attr(ith_calDA,"condition"),sep=""))
        LIST[[i]] <- NA

      } else {

        drMAT[i,] <- c(t(ith_calDA[,c(1L,2L),drop=FALSE]))
        LIST[[i]] <- ith_calDA

      } # end if.

      ###
      if (calRbfromK==TRUE) {

          TempTable[i,c(12L,13L)] <- c(ifelse(Kcontent[1L]<=9.17/38.13, 0, -9.17+38.13*Kcontent[1L]), 0)
          TempTable[i,c(16L,17L)] <- c(ifelse(inKcontent[1L]<=9.17/38.13, 0, -9.17+38.13*inKcontent[1L]), 0)

      } # end if.

    } # end for.

    ###
    if (N>1L) close(pb)
    if (!is.null(outpdf))  dev.off()
      
    ###
    TempTable[,-1L] <- apply(TempTable[,-1L],MARGIN=2L,function(x) round(as.numeric(x),digits=digits))
    TempTable[,32L:45L] <- round(drMAT, digits=digits)
    colnames(TempTable) <- c(name1, name2)

    ###
    TempTable0 <- data.frame(cbind(TempTable[,1L],
                                   paste(TempTable[,2L],TempTable[,3L],sep="\u00B1"),
                                   paste(TempTable[,4L],TempTable[,5L],sep="\u00B1"),
                                   paste(TempTable[,6L],TempTable[,7L],sep="\u00B1"),
                                   paste(TempTable[,8L],TempTable[,9L],sep="\u00B1"),
                                   paste(TempTable[,10L],TempTable[,11L],sep="\u00B1"),
                                   paste(TempTable[,12L],TempTable[,13L],sep="\u00B1"),
                                   paste(TempTable[,14L],TempTable[,15L],sep="\u00B1"),
                                   paste(TempTable[,16L],TempTable[,17L],sep="\u00B1"),
                                   paste(TempTable[,18L],TempTable[,19L],sep="\u00B1"),
                                   paste(TempTable[,20L],TempTable[,21L],sep="\u00B1"),
                                   paste(TempTable[,22L],TempTable[,23L],sep="\u00B1"),
                                   paste(TempTable[,24L],TempTable[,25L],sep="\u00B1"),
                                   paste(TempTable[,26L],TempTable[,27L],sep="\u00B1"),
                                   paste(TempTable[,28L],TempTable[,29L],sep="\u00B1"),
                                   paste(TempTable[,30L],TempTable[,31L],sep="\u00B1"),
                                   paste(TempTable[,32L],TempTable[,33L],sep="\u00B1"),
                                   paste(TempTable[,34L],TempTable[,35L],sep="\u00B1"),
                                   paste(TempTable[,36L],TempTable[,37L],sep="\u00B1"),
                                   paste(TempTable[,38L],TempTable[,39L],sep="\u00B1"),
                                   paste(TempTable[,40L],TempTable[,41L],sep="\u00B1"),
                                   paste(TempTable[,42L],TempTable[,43L],sep="\u00B1"),
                                   paste(TempTable[,44L],TempTable[,45L],sep="\u00B1")),
                                   stringsAsFactors=FALSE)

    ###
    colnames(TempTable0) <- c("SampleName","De(Gy)","GrainSize(um)","Ucontent(ppm)","Thcontent(ppm)","Kcontent(\u0025)","Rbcontent(ppm)",
                              "inKcontent(\u0025)","inRbcontent(ppm)","Wct(\u0025)","Depth(m)","Longitude(\u00BA)","Latitude(\u00BA)",
                              "Altitude(m a.s.l.)","AlphaValue", "BulkDensity(g/cm3)","Alpha.DR(Gy/ka)","Beta.DR(Gy/ka)","inBeta.DR(Gy/ka)",
                              "Gamma.DR(Gy/ka)","Cosmic.DR(Gy/ka)","Total.DR(Gy/ka)","Age(ka)")
    ###
    if (!is.null(outfile))  {

      write.csv(TempTable, file=paste(outfile,".csv",sep=""))
      ###write.csv(TempTable0, file=paste(outfile,"0.csv",sep=""))

      ## output html file.
      # insert <br> in header.
      modified_column_names <- sapply(names(TempTable0), function(x) gsub("\\(", "<br>(", x))
      
      ### start set html.
      html_content <- c("<html>", "<head>", "<style>",
        "table {border-collapse: collapse; width: 100%; font-family: 'Times New Roman', Times, serif; font-size: smaller;}",
        "th, td {border: none; text-align: center; padding: 8px; vertical-align: middle;}", 
        "th {border-top: 2px solid black; border-bottom: 1.5px solid black; font-weight: bold;}",  
        "tr:last-child td {border-bottom: 2px solid black;}",
        "th {white-space: pre-line;}", "</style>", "<title>", outfile, "</title>", 
        "</head>", "<body>", "<table>", "<tr><th>", paste(modified_column_names, "</th><th>", collapse = ""), "</th></tr>")
      
      ###
      for (i in 1:nrow(TempTable0)) {

        row_html <- paste("<tr><td>", paste(TempTable0[i,], "</td><td>", collapse=""), "</td></tr>")
        html_content <- c(html_content, row_html)

      } # end for.
      
      ###
      html_content <- c(html_content, "</table>", "</body>", "</html>")  
      # end set html.
      
      # write html file
      writeLines(html_content, paste(outfile,".html",sep=""))

    } # end if.

    ###
    names(LIST) <- as.character(TempTable[,1L])
    return(invisible(LIST))
     
  } # end if.

} # end function calDAbatch.default.
###
