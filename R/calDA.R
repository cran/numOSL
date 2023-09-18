######### 
calDA <-
function(dose, minGrainSize, maxGrainSize,
         Ucontent, Thcontent, Kcontent, Wct, depth, longitude, 
         latitude, altitude, alphaValue=0, inKcontent=0, bulkDensity=2.5, 
         cfType="Liritzis2013", rdcf=0, rba=0, ShallowGamma=TRUE, 
         nsim=5000, rejectNeg=TRUE, plot=TRUE, sampleName="")   {

  UseMethod("calDA") 

} # end function calDA.
### 2023.09.18.
calDA.default <-
function(dose, minGrainSize, maxGrainSize,
         Ucontent, Thcontent, Kcontent, Wct, depth, longitude,
         latitude, altitude, alphaValue=0, inKcontent=0, bulkDensity=2.5, 
         cfType="Liritzis2013", rdcf=0, rba=0, ShallowGamma=TRUE,  
         nsim=5000, rejectNeg=TRUE, plot=TRUE, sampleName="")  {

    stopifnot(length(dose)==2L, is.numeric(dose), all(is.finite(dose)),
              length(minGrainSize)==1L, is.numeric(minGrainSize), minGrainSize>0.0,
              length(maxGrainSize)==1L, is.numeric(maxGrainSize), maxGrainSize<1000.0,
              minGrainSize<maxGrainSize, 
              length(Ucontent)==2L, is.numeric(Ucontent), all(Ucontent>=0.0),
              length(Thcontent)==2L, is.numeric(Thcontent), all(Thcontent>=0.0),
              length(Kcontent)==2L, is.numeric(Kcontent), all(Kcontent>=0.0),
              length(Wct)==2L, is.numeric(Wct), all(Wct>=0.0),
              length(depth) %in% c(1L,2L), is.numeric(depth), all(depth>=0.0),
              length(longitude) %in% c(1L,2L), is.numeric(longitude), all(longitude>=-180),all(longitude<=180),
              length(latitude) %in% c(1L,2L), is.numeric(latitude), all(latitude>=-90), all(latitude<=90),
              length(altitude) %in% c(1L,2L), is.numeric(altitude), all(is.finite(altitude)),
              length(alphaValue) %in% c(1L,2L), is.numeric(alphaValue), all(alphaValue>=0.0),
              length(inKcontent) %in% c(1L,2L), is.numeric(inKcontent),  all(inKcontent>=0.0),
              length(bulkDensity) %in% c(1L,2L), is.numeric(bulkDensity), all(bulkDensity>=0.0),
              length(cfType)==1L, cfType %in% c("Liritzis2013","Guerin2011","AdamiecAitken1998"),
              length(rdcf)==1L, is.numeric(rdcf), rdcf>=0.0,
              length(rba)==1L, is.numeric(rba), rba>=0.0,
              length(ShallowGamma)==1L, is.logical(ShallowGamma), 
              length(nsim)==1L, is.numeric(nsim), nsim>=10,
              length(rejectNeg)==1L, is.logical(rejectNeg), 
              length(plot)==1L, is.logical(plot), 
              length(sampleName)==1L, is.character(sampleName))

    ###
    ### Function 1: BetaAttenuation.
    BetaAttenuation <- function(GrainSize1, GrainSize2, rba, rejectNeg)  {
      ###
      meanGrainSize <- (GrainSize1+GrainSize2)/2.0
      ###
      GSL <-  c(5,10,15,20,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,400,500,
                600,800,1000,1200,1400,1600,1800,2000,2500,3000,4000,5000,6000,8000,10000)
      
      Uphi <- c(0.012,0.0214,0.0296,0.0366,0.0475,0.0564,0.0642,0.0713,0.0779,0.084,0.0901,
                0.0957,0.106,0.117,0.127,0.137,0.146,0.169,0.189,0.23,0.263,0.295,0.351,0.4,
                0.442,0.479,0.512,0.541,0.568,0.623,0.667,0.731,0.773,0.803,0.842,0.866)
      
      THphi <- c(0.0225,0.0366,0.0484,0.0582,0.0743,0.0875,0.0988,0.1088,0.1181,0.1269,
                 0.1351,0.1427,0.158,0.171,0.184,0.195,0.206,0.229,0.251,0.288,0.32,0.348,0.399,
                 0.443,0.482,0.518,0.55,0.578,0.604,0.659,0.701,0.76,0.798,0.825,0.859,0.879)
      
      Kphi <- c(0.0018,0.0035,0.0053,0.0071,0.0106,0.0141,0.0177,0.0212,0.0248,0.0283,0.0318,
                0.0354,0.0424,0.0494,0.0563,0.0633,0.0702,0.0877,0.1052,0.1402,0.1748,0.209,0.275,
                0.337,0.394,0.447,0.493,0.535,0.571,0.643,0.696,0.765,0.807,0.837,0.873,0.896)
      
      ###
      Uabsorption0 <-  spline(x=GSL, y=Uphi, xout=meanGrainSize)$y
      Thabsorption0 <- spline(x=GSL, y=THphi, xout=meanGrainSize)$y
      Kabsorption0 <- spline(x=GSL, y=Kphi, xout=meanGrainSize)$y
      
      ###
      if (rejectNeg==FALSE)  {

        Uabsorption <- rnorm(n=1L, mean=Uabsorption0, sd=rba*Uabsorption0)
        Thabsorption <- rnorm(n=1L, mean=Thabsorption0, sd=rba*Thabsorption0)
        Kabsorption <- rnorm(n=1L, mean=Kabsorption0, sd=rba*Kabsorption0)

      } else {

        n1 <- n2 <- n3 <- 0

        ###
        repeat {
 
          Uabsorption <- rnorm(n=1L, mean=Uabsorption0, sd=rba*Uabsorption0)
          if (Uabsorption>=0.0) break
          n1 <- n1+1
          if (n1>10) stop("BetaAttenuation(), inefficient sampling of U-Beta-Absorption, try a smaller [rba]!")

        } # end repeat.

        ###
        repeat {

          Thabsorption <- rnorm(n=1L, mean=Thabsorption0, sd=rba*Thabsorption0)
          if (Thabsorption>=0.0) break
          n2 <- n2+1
          if (n2>10) stop("BetaAttenuation(), inefficient sampling of Th-Beta-Absorption, try a smaller [rba]!")

        } # end repeat.

        ###
        repeat {

          Kabsorption <- rnorm(n=1L, mean=Kabsorption0, sd=rba*Kabsorption0)
          if (Kabsorption>=0.0) break
          n3 <- n3+1
          if (n3>10) stop("BetaAttenuation(), inefficient sampling of K-Beta-Absorption, try a smaller [rba]!")

        } # end repeat. 

      } # end if.
      
      ###
      out1 <- c(Uabsorption, Thabsorption, Kabsorption)
      return(out1)

    } # end function BetaAttenuation.
    
    ### Function 2: AlphaAttenuation.
    AlphaAttenuation <- function(GrainSize1, GrainSize2, rba, rejectNeg)  {
      
      meanGrainSize <- (GrainSize1+GrainSize2) / 2.0
      
      ###
      GSL <-  c(1,   2,   3,   4,   5,   6,   7,   8,   9,   10, 
                20,  30,  40,  50,  60,  70,  80,  90,  100,
                200, 300, 400, 500, 600, 700, 800, 900, 1000)
      
      Uphi <- c(0.02,  0.04, 0.06, 0.08,  0.11,  0.125, 0.149, 0.17, 0.19, 0.21,
                0.42,  0.56, 0.65, 0.72,  0.76,  0.79,  0.82,  0.84, 0.85, 0.92, 
                0.945, 0.96, 0.97, 0.975, 0.977, 0.98,  0.98,  0.98)
      
      THphi <- c(0.02, 0.045, 0.055, 0.07,  0.09, 0.11, 0.13,  0.14, 0.16, 0.17,
                 0.34, 0.48,  0.58,  0.635, 0.71, 0.74, 0.77,  0.79, 0.82, 0.905,
                 0.94, 0.95,  0.955, 0.965, 0.97, 0.97, 0.975, 0.98)
      
      ###
      Uabsorption0 <-  spline(x=GSL, y=Uphi, xout=meanGrainSize)$y
      Thabsorption0 <- spline(x=GSL, y=THphi, xout=meanGrainSize)$y
      
      ###
      if (rejectNeg==FALSE) {

        Uabsorption <- rnorm(n=1L, mean=Uabsorption0, sd=rba*Uabsorption0)      
        Thabsorption <- rnorm(n=1L, mean=Thabsorption0, sd=rba*Thabsorption0)

      } else {

          n1 <- n2 <- 0

          ###
          repeat {

            Uabsorption <- rnorm(n=1L, mean=Uabsorption0, sd=rba*Uabsorption0)
            if (Uabsorption>=0.0) break
            n1 <- n1+1
            if (n1>10) stop("AlphaAttenuation(), inefficient sampling of U-Alpha-Absorption, try a smaller [rba]!")

          } # end repeat. 

          ###
          repeat {

            Thabsorption <- rnorm(n=1L, mean=Thabsorption0, sd=rba*Thabsorption0)
            if (Thabsorption>=0.0) break
            n2 <- n2+1
            if (n2>10) stop("AlphaAttenuation(), inefficient sampling of Th-Alpha-Absorption, try a smaller [rba]!")

          } # end repeat. 

      } # end if.
      
      out1 <- c(Uabsorption, Thabsorption)
      return(out1)
    
    } # end function  AlphaAttenuation.
    
    ### Function 3: AlphaBetaGammaDoseRate.
    AlphaBetaGammaDoseRate <- function(minGrainSize, maxGrainSize,
                                       Ucontent, Thcontent, Kcontent, depth, inKcontent, 
                                       bulkDensity, cfType, rdcf, rba, ShallowGamma, rejectNeg)  {
      
      # change unit from "m" to "cm"
      depth <- depth*100
      
      ###
      if (cfType=="Liritzis2013") {
          
        # Liritzis et al., 2013.
        cua <- 2.7930
        cua_err <- ifelse(rdcf>0, 0.0110, 0)
          
        ctha <- 0.7375
        ctha_err <- ifelse(rdcf>0, 0.0026, 0)
          
        cub <- 0.1459
        cub_err <- ifelse(rdcf>0, 0.0004, 0)
          
        cug <- 0.1118
        cug_err <- ifelse(rdcf>0, 0.0002, 0)
          
        cthb <- 0.0275
        cthb_err <- ifelse(rdcf>0, 0.0009, 0)
          
        cthg <- 0.0481
        cthg_err <- ifelse(rdcf>0, 0.0002, 0)
          
        ckb <- 0.8011
        ckb_err <- ifelse(rdcf>0, 0.0073, 0)
          
        ckg <- 0.2498
        ckg_err <- ifelse(rdcf>0, 0.0048, 0)
        
      } else if (cfType=="Guerin2011") {
        
        # Guerin2011.
        cua <- 2.795
        cua_err <- cua*rdcf
        
        ctha <- 0.7375
        ctha_err <- ctha*rdcf
        
        cub <- 0.1457
        cub_err <- cub*rdcf
        
        cug <- 0.1116
        cug_err <- cug*rdcf
        
        cthb <- 0.0277
        cthb_err <- cthb*rdcf
        
        cthg <- 0.0479
        cthg_err <- cthg*rdcf
        
        ckb <- 0.7982
        ckb_err <- ckb*rdcf
        
        ckg <- 0.2491
        ckg_err <- ckg*rdcf
        
      } else if (cfType=="AdamiecAitken1998"){
        
        # Adamiec and Aitken1998.
        cua <- 2.7800
        cua_err <- cua*rdcf
        
        ctha <- 0.7320
        ctha_err <- ctha*rdcf
        
        cub <- 0.1460
        cub_err <- cub*rdcf
        
        cug <- 0.1130
        cug_err <- cug*rdcf
        
        cthb <- 0.0273
        cthb_err <- cthb*rdcf
        
        cthg <- 0.0476
        cthg_err <- cthg*rdcf
        
        ckb <- 0.7820
        ckb_err <- ckb*rdcf
        
        ckg <- 0.2430
        ckg_err <- ckg*rdcf
         
      } # end if.
      
      ###    
      if(depth/bulkDensity*2 > 41.0){

        U_factor <- Th_factor <- K_factor <- 1.0

      } else {

        if (ShallowGamma==TRUE) {
        
          shallow_gamma_dose_factor <-
          cbind(c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,
                  11,12,13,14,15,16,17,18,20,22,24,25,26,28,30,31,32,33,34,35,35.5,
                  36,36.5,37,37.5,38,38.5,39,39.5,40,40.5,41),                                  #depth(cm)
                c(0.5,0.5607,0.6022,0.6366,0.6663,0.6924,0.7156,0.7364,0.7552,0.7721,
                  0.7876,0.8016,0.8145,0.8263,0.8373,0.8473,0.8567,0.8653,0.8734,0.8809,
                  0.8879,0.9006,0.9118,0.9216,0.9303,0.938,0.9449,0.951,0.9564,0.9655,0.9727,
                  0.9784,0.9808,0.9829,0.9865,0.9893,0.9905,0.9916,0.9925,0.9934,0.9941,0.9944,
                  0.9948,0.995,0.9953,0.9956,0.9958,0.9961,0.9963,0.9965,0.9967,0.9969,0.9971), #U_factor
                c(0.5,0.5577,0.5974,0.6306,0.6594,0.6849,0.7076,0.7281,0.7466,0.7634,0.7787,
                  0.7927,0.8055,0.8174,0.8283,0.8384,0.8478,0.8564,0.8646,0.8722,0.8793,0.8921,
                  0.9035,0.9135,0.9224,0.9304,0.9375,0.9438,0.9495,0.9592,0.9669,0.9732,0.9759,
                  0.9783,0.9824,0.9857,0.9871,0.9884,0.9895,0.9905,0.9915,0.9919,0.9923,0.9927,
                  0.993,0.9934,0.9937,0.994,0.9943,0.9946,0.9949,0.9951,0.9953),                #Th_factor
                c(0.5,0.5555,0.5938,0.6258,0.6536,0.6782,0.7003,0.7202,0.7382,0.7547,0.7697,
                  0.7836,0.7964,0.8083,0.8193,0.8295,0.8391,0.848,0.8564,0.8643,0.8716,0.8851,
                  0.8971,0.9078,0.9174,0.9259,0.9335,0.9404,0.9465,0.957,0.9655,0.9722,0.9752,
                  0.9778,0.9822,0.9858,0.9873,0.9886,0.9898,0.9909,0.9919,0.9923,0.9928,0.9932,
                  0.9935,0.9939,0.9942,0.9946,0.9949,0.9951,0.9954,0.9957,0.9959))              #K_factor
        
          U_factor <- spline(x=shallow_gamma_dose_factor[,1], y=shallow_gamma_dose_factor[,2], xout=depth/bulkDensity*2)$y
          Th_factor <- spline(x=shallow_gamma_dose_factor[,1], y=shallow_gamma_dose_factor[,3], xout=depth/bulkDensity*2)$y
          K_factor <- spline(x=shallow_gamma_dose_factor[,1], y=shallow_gamma_dose_factor[,4], xout=depth/bulkDensity*2)$y

        } else {

          U_factor <- Th_factor <- K_factor <- 1.0

        } # end if.

      } # end if.
      
      ###
      if (rejectNeg==FALSE) {

        UAlpha <- rnorm(n=1L, mean=cua, sd=cua_err)*Ucontent
        UBeta <- rnorm(n=1L, mean=cub, sd=cub_err)*Ucontent
        UGamma <- rnorm(n=1L, mean=cug, sd=cug_err)*Ucontent*U_factor
        
        ###
        Thalpha <- rnorm(n=1L, mean=ctha, sd=ctha_err)*Thcontent
        ThBeta <- rnorm(n=1L, mean=cthb, sd=cthb_err)*Thcontent
        ThGamma <- rnorm(n=1L, mean=cthg, sd=cthg_err)*Thcontent*Th_factor
        
        ###
        fxfxfx <- rnorm(n=1L, mean=ckb, sd=ckb_err)
        KBeta <- fxfxfx*Kcontent
        KBeta_internal <- fxfxfx*inKcontent
        KGamma <- rnorm(n=1L, mean=ckg, sd=ckg_err)*Kcontent*K_factor

      } else {

        n1 <- n2 <- n3 <- n4 <- n5 <- n6 <- n7 <- n8 <- 0

        ###
        repeat {

          UAlpha <- rnorm(n=1L, mean=cua, sd=cua_err)*Ucontent
          if (UAlpha>=0.0) break
          n1 <- n1+1
          if (n1>10) stop("AlphaBetaGammaDoseRate(), inefficient sampling of U-Alpha CF, try a smaller [rdcf]!")

        } # end repeat.

        ###
        repeat {

          UBeta <- rnorm(n=1L, mean=cub, sd=cub_err)*Ucontent
          if (UBeta>=0.0) break
          n2 <- n2+1
          if (n2>10) stop("AlphaBetaGammaDoseRate(), inefficient sampling of U-Beta CF, try a smaller [rdcf]!")

        } # end repeat.

        ###
        repeat {

          UGamma <- rnorm(n=1L, mean=cug, sd=cug_err)*Ucontent*U_factor
          if (UGamma>=0.0) break
          n3 <- n3+1
          if (n3>10) stop("AlphaBetaGammaDoseRate(), inefficient sampling of U-Gamma CF, try a smaller [rdcf]!")

        } # end repeat.

        ###
        repeat {

          Thalpha <- rnorm(n=1L, mean=ctha, sd=ctha_err)*Thcontent
          if (Thalpha>=0.0) break
          n4 <- n4+1
          if (n4>10) stop("AlphaBetaGammaDoseRate(), inefficient sampling of Th-Alpha CF, try a smaller [rdcf]!")

        } # end repeat.

        ###
        repeat {

          ThBeta <- rnorm(n=1L, mean=cthb, sd=cthb_err)*Thcontent
          if (ThBeta>=0.0) break
          n5 <- n5+1
          if (n5>10) stop("AlphaBetaGammaDoseRate(), inefficient sampling of Th-Beta CF, try a smaller [rdcf]!")

        } # end repeat.

        ###
        repeat {

          ThGamma <- rnorm(n=1L, mean=cthg, sd=cthg_err)*Thcontent*Th_factor
          if (ThGamma>=0.0) break
          n6 <- n6+1
          if (n6>10) stop("AlphaBetaGammaDoseRate(), inefficient sampling of Th-Gamma CF, try a smaller [rdcf]!")

        } # end repeat.

        ###
        repeat {

          fxfxfx <- rnorm(n=1L, mean=ckb, sd=ckb_err)
          if (fxfxfx>=0.0) break
          n7 <- n7+1
          if (n7>10) stop("AlphaBetaGammaDoseRate(), inefficient sampling of K-Beta CF, try a smaller [rdcf]!")

        } # end repeat.
        KBeta <- fxfxfx*Kcontent
        KBeta_internal <- fxfxfx*inKcontent

        ###
        repeat {

          KGamma <- rnorm(n=1L, mean=ckg, sd=ckg_err)*Kcontent*K_factor
          if (KGamma>=0.0) break
          n8 <- n8+1
          if (n8>10) stop("AlphaBetaGammaDoseRate(), inefficient sampling of K-Gamma CF, try a smaller [rdcf]!")

        } # end repeat.
    
      } # end if.

      ###
      AttenuationFactors1 <- 1.0 - AlphaAttenuation(minGrainSize, maxGrainSize, rba, rejectNeg)  
      alphaDoseRate <- UAlpha*AttenuationFactors1[1] + Thalpha*AttenuationFactors1[2]
      
      ###
      AttenuationFactors2 <- 1.0 - BetaAttenuation(minGrainSize, maxGrainSize, rba, rejectNeg)   
      betaDoseRate <- UBeta*AttenuationFactors2[1] + ThBeta*AttenuationFactors2[2] + KBeta*AttenuationFactors2[3]
      internalBetaDoseRate <- KBeta_internal*(1.0-AttenuationFactors2[3])
      
      ###
      gammaDoseRate <- UGamma + ThGamma + KGamma

      ###
      out2 <- c(alphaDoseRate, betaDoseRate, gammaDoseRate, internalBetaDoseRate)
      return(out2)
      
    } # end function AlphaBetaGammaDoseRate.
    
    ### Function 4: CosmicDoseRate.
    CosmicDoseRate <- function(depth, longitude, latitude, altitude, bulkDensity)  {
      
      lambda <- asin(0.203*cos(latitude*pi/180.0)*
                    cos((longitude-291.0)*pi/180.0)+
                    0.979*sin(latitude*pi/180.0))*180.0/pi
      
      ###
      if (depth*bulkDensity>1.67) {
        
        D0 <- 6072.0/(((depth*bulkDensity+11.6)^1.68+75.0)*
                        (depth*bulkDensity+212.0))*exp(-0.00055*(depth*bulkDensity))
        
      } else {
        
        D0 <- 0.295 - 0.207*depth*bulkDensity + 0.221*(depth*bulkDensity)^2 - 
              0.135*(depth*bulkDensity)^3 + 0.0321*(depth*bulkDensity)^4
        
      } # end if.
      
      ###
      func_J <- function(lambda) {
        
        J_df <- cbind(c(0.1888087,5.1724030, 10.0451199, 15.0282434, 20.0476982, 25.0667214, 27.9078645, 
                        30.7488899, 34.7706058, 39.9761117, 45.1448546, 49.9073220, 55.0021075, 59.8015732, 
                        65.0071183, 69.8066625, 74.7168486),
                      c(0.5213598, 0.5324123, 0.5456987, 0.5656778, 0.5968145, 0.6361339, 0.6680608, 0.7022193, 
                        0.7504867, 0.7540955, 0.7547295, 0.7583477, 0.7612150, 0.7633447, 0.7662096, 0.7668515, 
                        0.7697227))

        ###
        if (lambda<=34.7706058) {

          Jv <- spline(x=J_df[1:9,1], y=J_df[1:9,2], xout=lambda)$y

        } else {
         
          pab <- coef(lm(y~x, data=data.frame(x=J_df[9:17,1],y=J_df[9:17,2])))
          Jv <- pab[1] + pab[2]*lambda
          
        } # end if.
        
        return(Jv)

      } # end function func_J.
      
      ###
      func_F <- function(lambda) {
        
        F_df <- cbind(c(0.2689652,  5.2163851, 10.0902791, 15.1122842, 20.0976048, 25.1203945, 27.5579301, 
                        30.0325816, 32.9505856, 35.1665234, 37.7513853, 45.0616770, 59.9405033, 54.9562813, 
                        49.9351002),
                      c(0.4015919, 0.3985112, 0.3894809, 0.3722647, 0.3505860, 0.3184920, 0.3028185, 0.2834248, 
                        0.2580705, 0.2438895, 0.2349079, 0.2317770, 0.2292291, 0.2300789, 0.2316734))

        ###
        if (lambda<=37.7513853) {

          Fv <- spline(x=F_df[1:11,1], y=F_df[1:11,2], xout=lambda)$y

        } else {

          pab <- coef(lm(y~x, data=data.frame(x=F_df[11:15,1],y=F_df[11:15,2])))
          Fv <- pab[1] + pab[2]*lambda

        } # end if.
          
        return(Fv)
 
      } # end function func_F.
      
      ###
      func_H <- function(lambda) {
        
        H_df <- cbind(c(0.2480531,  5.1971993, 10.0359390, 15.0968649, 20.0477767, 25.1110174, 28.6972671, 
                        30.0285797, 31.7669525, 33.7273153, 35.1687990),
                      c(4.399043, 4.381137, 4.359885, 4.332307, 4.297663, 4.248141, 4.200494, 4.179651, 
                        4.149877, 4.115637, 4.100372))

        ###
        if (lambda<=35.1687990) {

          Hv <- spline(x=H_df[1:11,1], y=H_df[1:11,2], xout=lambda)$y

        } else {

          Hv <- 4.100372

        } # end if.
        
        return(Hv)

      } # end function func_H.
      
      ###
      FF <- func_F(abs(lambda))
      JJ <- func_J(abs(lambda))
      HH <- func_H(abs(lambda))
      
      ###
      cosmicDoseRate <- D0*(FF+JJ*exp(altitude/HH/1000.0))  
      return(cosmicDoseRate)
      
    } # end function CosmicDoseRate.
    
    ### Function 5: calDoseRate.
    calDoseRate <- function(minGrainSize, maxGrainSize,
                            Ucontent, Thcontent, Kcontent, Wct, depth, longitude,  
                            latitude, altitude, alphaValue, inKcontent, bulkDensity,  
                            rdcf, rba, cfType, ShallowGamma)  {
      ### 
      abgDoseRate <- AlphaBetaGammaDoseRate(minGrainSize, maxGrainSize,
                                            Ucontent, Thcontent, Kcontent, depth, inKcontent, 
                                            bulkDensity, cfType, rdcf, rba, ShallowGamma, rejectNeg)
      
      ###
      alphaDoseRate <- alphaValue[1] * abgDoseRate[1] / (1.0+1.5*Wct/100)
      betaDoseRate <- abgDoseRate[2] / (1.0+1.25*Wct/100)
      gammaDoseRate <- abgDoseRate[3] /(1.0+1.14*Wct/100)
      internalBetaDoseRate <- abgDoseRate[4]
      
      ###
      cosmicDoseRate <- CosmicDoseRate(depth[1], longitude[1], latitude[1], altitude[1], bulkDensity[1]) 
      
      ###
      totalDoseRate <- alphaDoseRate + betaDoseRate + gammaDoseRate +  cosmicDoseRate + internalBetaDoseRate
      DR <- c(alphaDoseRate, betaDoseRate, internalBetaDoseRate, gammaDoseRate, cosmicDoseRate, totalDoseRate)
      return(DR)
     
    } # end function calDoseRate.
    
    
    ### Function 6: simDoseRateAge.
    simDoseRateAge <- function(dose, minGrainSize, maxGrainSize,
                               Ucontent, Thcontent, Kcontent, Wct, depth, longitude,
                               latitude, altitude, alphaValue, inKcontent, bulkDensity,  
                               rdcf, rba, cfType, ShallowGamma, nsim)  {
      
      alphaDoseRates <- betaDoseRates <- internalBetaDoseRates  <- gammaDoseRates <- 
      cosmicDoseRates <- DoseRates <- Ages <- vector(length=nsim)
      
      ###
      siminKcontent <- inKcontent[1L]
      simdepth <- depth[1L] 
      simlongitude <- longitude[1L]
      simlatitude <- latitude[1L]
      simaltitude <- altitude[1L]
      simbulkDensity <- bulkDensity[1L] 
      simalphaValue <- alphaValue[1L]
      
      ###
      for (i in 1:nsim)  {

        ###1.
        simdose <- rnorm(n=1L, mean=dose[1L], sd=dose[2L])
        
        if (rejectNeg==FALSE) {

          ###2.
          simUcontent <- rnorm(n=1L, mean=Ucontent[1L], sd=Ucontent[2L])

          ###3.
          simThcontent <- rnorm(n=1L, mean=Thcontent[1L], sd=Thcontent[2L])

          ###4.
          simKcontent <- rnorm(n=1L, mean=Kcontent[1L], sd=Kcontent[2L])

          ###5.
          simWct <- rnorm(n=1L, mean=Wct[1L], sd=Wct[2L])

          ###6.                            
          if (length(depth)==2L)  {

            simdepth <- rnorm(n=1L, mean=depth[1L], sd=depth[2L])

          } # end if.

          ###7.
          if (length(longitude)==2L)  {

            simlongitude <- rnorm(n=1L, mean=longitude[1L], sd=longitude[2L])

          } # end if.

          ###8.
          if (length(latitude)==2L)  {

            simlatitude <- rnorm(n=1L, mean=latitude[1L], sd=latitude[2L])

          } # end if.

          ###9.
          if (length(altitude)==2L)  {

            simaltitude <- rnorm(n=1L, mean=altitude[1L], sd=altitude[2L])

          } # end if.

          ###10.
          if (length(bulkDensity)==2L)  {
            
            simbulkDensity <- rnorm(n=1L, mean=bulkDensity[1L], sd=bulkDensity[2L])

          } # end if.

          ###11.
          if (length(alphaValue)==2L)  {

            simalphaValue <- rnorm(n=1L, mean=alphaValue[1L], sd=alphaValue[2L])

          } # end if.

          ###12.
          if (length(inKcontent)==2L)  {
            
            siminKcontent <- rnorm(n=1L, mean=inKcontent[1L], sd=inKcontent[2L])

          } # end if.

        } else {

          n1 <- n2 <- n3 <- n4 <- n5 <- n6 <- n7 <- n8 <- 0

          ###2.
          if (Ucontent[1L]>0.0) {
            repeat {

              simUcontent <- rnorm(n=1L, mean=Ucontent[1L], sd=Ucontent[2L])
              if (simUcontent>=0) break
              n1 <- n1+1
              if (n1>10) stop("simDoseRateAge(), inefficient sampling, standard error of [Ucontent] is very large!")

            } # end repeat.
          } else {
       
             simUcontent <- abs(rnorm(n=1L, mean=Ucontent[1L], sd=Ucontent[2L]))

          } # end if. 

          ###3.
          if (Thcontent[1L]>0.0) {

            repeat {

              simThcontent <- rnorm(n=1L, mean=Thcontent[1L], sd=Thcontent[2L])
              if (simThcontent>=0) break
              n2 <- n2+1
              if (n2>10) stop("simDoseRateAge(), inefficient sampling, standard error of [Thcontent] is very large!")

            } # end repeat.

          } else {

            simThcontent <- abs(rnorm(n=1L, mean=Thcontent[1L], sd=Thcontent[2L]))

          } # end if. 

          ###4.
          if (Kcontent[1L]>0.0) {

            repeat {

              simKcontent <- rnorm(n=1L, mean=Kcontent[1L], sd=Kcontent[2L])
              if (simKcontent>=0) break
              n3 <- n3+1
              if (n3>10) stop("simDoseRateAge(), inefficient sampling, standard error of [Kcontent] is very large!!")
            
            } # end repeat.

          } else {

            simKcontent <- abs(rnorm(n=1L, mean=Kcontent[1L], sd=Kcontent[2L]))

          } # end if. 

          ###5.
          if (Wct[1L]>0.0) {

            repeat {

              simWct <- rnorm(n=1L, mean=Wct[1L], sd=Wct[2L])
              if (simWct>=0) break
              n4 <- n4+1
              if (n4>10) stop("simDoseRateAge(), inefficient sampling, standard error of [Wct] is very large!")

            } # end if.

          } else {

            simWct <- abs(rnorm(n=1L, mean=Wct[1L], sd=Wct[2L]))

          } # end if. 

          ###6. 
          if (depth[1L]>0.0) {   
                       
            if (length(depth)==2L)  {

              repeat {

                simdepth <- rnorm(n=1L, mean=depth[1L], sd=depth[2L])
                if (simdepth>=0) break
                n5 <- n5+1
                if (n5>10) stop("simDoseRateAge(), inefficient sampling, standard error of [depth] is very large!")

              } # end repeat.

            } # end if.

          } else {

            simdepth <- abs(rnorm(n=1L, mean=depth[1L], sd=depth[2L]))

          } # end if. 

          ###7.
          if (length(longitude)==2L)  {

            simlongitude <- rnorm(n=1L, mean=longitude[1L], sd=longitude[2L])

          } # end if.

          ###8.
          if (length(latitude)==2L)  {

            simlatitude <- rnorm(n=1L, mean=latitude[1L], sd=latitude[2L])

          } # end if.

          ###9.
          if (length(altitude)==2L)  {

            simaltitude <- rnorm(n=1L, mean=altitude[1L], sd=altitude[2L])
            
          } # end if.

          ###10.
          if (length(bulkDensity)==2L)  {

            if (bulkDensity[1L]>0.0) {

              repeat {
            
                simbulkDensity <- rnorm(n=1L, mean=bulkDensity[1L], sd=bulkDensity[2L])
                if (simbulkDensity>=0) break
                n6 <- n6+1
                if (n6>10) stop("simDoseRateAge(), inefficient sampling, standard error of [bulkDensity] is very large!")

              } # end repeat.

            } else {

              simbulkDensity <- abs(rnorm(n=1L, mean=bulkDensity[1L], sd=bulkDensity[2L]))

            } # end if. 

          } # end if.

          ###11.
          if (length(alphaValue)==2L)  {

            if (alphaValue[1L]>0.0) {

              repeat {

                simalphaValue <- rnorm(n=1L, mean=alphaValue[1L], sd=alphaValue[2L])
                if (simalphaValue>=0) break
                n7 <- n7+1
                if (n7>10) stop("simDoseRateAge(), inefficient sampling, standard error of [alphaValue] is very large!")

              } # end repeat.

            } else {

              simalphaValue <- abs(rnorm(n=1L, mean=alphaValue[1L], sd=alphaValue[2L]))

            } # end if. 

          } # end if.

          ###12.
          if (length(inKcontent)==2L)  {

            if (inKcontent[1L]>0.0) {

              repeat {

                siminKcontent <- rnorm(n=1L, mean=inKcontent[1L], sd=inKcontent[2L])
                if (siminKcontent>=0) break
                n8 <- n8+1
                if (n8>10) stop("simDoseRateAge(), inefficient sampling, standard error of [inKcontent] is very large!")

              } # end repeat.
          
            } else {

              siminKcontent <- abs(rnorm(n=1L, mean=inKcontent[1L], sd=inKcontent[2L]))

            } # end if. 

          } # end if.

        } # end if.
        
        ###
        randomDoseRate <- calDoseRate(minGrainSize=minGrainSize, maxGrainSize=maxGrainSize,
                                      Ucontent=simUcontent, Thcontent=simThcontent, Kcontent=simKcontent, Wct=simWct,   
                                      depth=simdepth, longitude=simlongitude, latitude=simlatitude, altitude=simaltitude,
                                      alphaValue=simalphaValue, inKcontent=siminKcontent, bulkDensity=simbulkDensity,  
                                      rdcf=rdcf, rba=rba, cfType=cfType, ShallowGamma=ShallowGamma)
        
        ###
        alphaDoseRates[i] <- randomDoseRate[1]
        betaDoseRates[i] <- randomDoseRate[2]
        internalBetaDoseRates[i] <- randomDoseRate[3]
        gammaDoseRates[i] <- randomDoseRate[4]
        cosmicDoseRates[i] <- randomDoseRate[5]  
        DoseRates[i] <- randomDoseRate[6]
        Ages[i] <- simdose / randomDoseRate[6] 

      } # end for.
    
      ###
      sdalphaDoseRate <- sd(alphaDoseRates)
      sdbetaDoseRate <- sd(betaDoseRates)
      sdinternalBetaDoseRate <- sd(internalBetaDoseRates)
      sdgammaDoseRate <- sd(gammaDoseRates)
      sdcosmicDoseRate <- sd(cosmicDoseRates)
      sdDoseRate <- sd(DoseRates)
      sdAge <- sd(Ages)
      
      ###
      output <- list("SD"=c(sdalphaDoseRate, sdbetaDoseRate, 
                            sdinternalBetaDoseRate, sdgammaDoseRate,
                            sdcosmicDoseRate, sdDoseRate, sdAge),
                     "simaDR"=alphaDoseRates, "simbDR"=betaDoseRates,
                     "siminbDR"=internalBetaDoseRates, "simgDR"=gammaDoseRates, 
                     "simcDR"=cosmicDoseRates, "simDR"=DoseRates, "simAGE"=Ages)
      return(output)

    } # end function simDoseRateAge.

    ###
    if (minGrainSize<=63 && maxGrainSize<=63 && !alphaValue>0) {
   
      cat("<",sampleName, "> Note the grain size is <=63 um (fine/medium grain), but the alpha efficiency is zero!\n",sep="")
      
    } # end if.

    ###
    if (minGrainSize>=63 && maxGrainSize>=63 && alphaValue>0) {
   
      cat("<",sampleName,"> Note the grain size is >=63 um (coarse grain), but the alpha efficiency is not zero!\n",sep="")
      
    } # end if.

    ###
    actualDoseRate <- calDoseRate(minGrainSize=minGrainSize, maxGrainSize=maxGrainSize,
                                  Ucontent=Ucontent[1L], Thcontent=Thcontent[1L], Kcontent=Kcontent[1L], Wct=Wct[1L], 
                                  depth=depth[1L], longitude=longitude[1L], latitude=latitude[1L], altitude=altitude[1L],
                                  alphaValue=alphaValue[1L], inKcontent=inKcontent[1L], bulkDensity=bulkDensity[1L], 
                                  rdcf=0, rba=0, cfType=cfType, ShallowGamma=ShallowGamma)
    ###
    errorDoseRate <- simDoseRateAge(dose=dose, minGrainSize=minGrainSize, maxGrainSize=maxGrainSize,
                                    Ucontent=Ucontent, Thcontent=Thcontent, Kcontent=Kcontent, Wct=Wct,   
                                    depth=depth, longitude=longitude, latitude=latitude, altitude=altitude, 
                                    alphaValue=alphaValue, inKcontent=inKcontent, bulkDensity=bulkDensity, 
                                    rdcf=rdcf, rba=rba, cfType=cfType, ShallowGamma=ShallowGamma, nsim=nsim)

    ###
    lu1 <- quantile(errorDoseRate$simaDR, probs=c(0.025,0.975))
    lu2 <- quantile(errorDoseRate$simbDR, probs=c(0.025,0.975))
    lu3 <- quantile(errorDoseRate$siminbDR, probs=c(0.025,0.975))
    lu4 <- quantile(errorDoseRate$simgDR, probs=c(0.025,0.975))
    lu5 <- quantile(errorDoseRate$simcDR, probs=c(0.025,0.975))
    lu6 <- quantile(errorDoseRate$simDR, probs=c(0.025,0.975))
    lu7 <- quantile(errorDoseRate$simAGE, probs=c(0.025,0.975))

    ###
    if (plot==TRUE) {
 
      opar <- par("mfrow", "mar", "mgp")
      on.exit(par(opar))
 
      ###
      par("mgp"=c(2.5,1,0))
      layout(cbind(c(1,1,2,2),c(1,1,2,2)), respect=TRUE)

      ###
      par("mar"=c(4.1, 4.1, 4.1, 2.1))
      dDoseRates <- density(errorDoseRate$simDR)
      plot(dDoseRates$x, dDoseRates$y, type="n", xlab="Total dose rate (Gy/ka)",
           ylab="Density", cex.lab=1.5, cex.axis=1.5, main=sampleName, cex.main=1.5)
       
      ###
      xTicks<-axTicks(side=1L)
      maxYx<-dDoseRates$x[which.max(dDoseRates$y)]

      ###
      polygon(dDoseRates$x, dDoseRates$y, col="skyblue", border="skyblue")
      rug(errorDoseRate$simDR, quiet=TRUE, col="red")
        
      ###
      legend(ifelse(maxYx>median(xTicks),"topleft","topright"), 
             legend=c(paste("N=", nsim, sep=""),
                      paste("DR=",round(actualDoseRate[6],2L)," Gy/ka", sep=""),
                      paste("Se.DR=",round(errorDoseRate$SD[6],2L)," Gy/ka",sep=""),
                      paste("Mean.DR=",round(mean(errorDoseRate$simDR),2L)," Gy/ka",sep=""),
                      paste("L95.DR=",round(lu6[1],2L)," Gy/ka",sep=""),
                      paste("U95.DR=",round(lu6[2],2L)," Gy/ka",sep="")), 
             cex=1.1, bty="n")
      box(lwd=1.5)
     
      ###
      par("mar"=c(6.1, 4.1, 1.1, 2.1))
      dAges <- density(errorDoseRate$simAGE)
      plot(dAges$x, dAges$y, xlab="Age (ka)", ylab="Density",
           cex.lab=1.5, cex.axis=1.5, type="n")
      
      ###
      xTicks<-axTicks(side=1L)
      maxYx<-dAges$x[which.max(dAges$y)]
        
      ###
      polygon(dAges$x, dAges$y, col="purple", border="purple")
      rug(errorDoseRate$simAGE, quiet=TRUE, col="red")
      
      ###
      legend(ifelse(maxYx>median(xTicks),"topleft","topright"), 
             legend=c(paste("N=", nsim, sep=""),
                      paste("Age=",round(dose[1]/actualDoseRate[6],2L)," ka", sep=""),
                      paste("Se.Age=",round(errorDoseRate$SD[7],2L)," ka",sep=""), 
                      paste("Mean.Age=",round(mean(errorDoseRate$simAGE),2L)," ka", sep=""),
                      paste("L95.Age=",round(lu7[1],2L)," ka",sep=""),
                      paste("U95.Age=",round(lu7[2],2L)," ka",sep="")),
             cex=1.1, bty="n")
      box(lwd=1.5)

    } # end if.

    ###
    output <- cbind(cbind(c(actualDoseRate,dose[1]/actualDoseRate[6]), errorDoseRate$SD),
                    rbind(lu1,lu2,lu3,lu4,lu5,lu6,lu7))
    rownames(output) <- c("Alpha.DR", "Beta.DR", "inBeta.DR", "Gamma.DR", "Cosmic.DR", "Total.DR", "Age")
    colnames(output) <- c("Pars", "Se.Pars", "Lower95", "Upper95")
    return(output)

} # end function calDA.
###
