subroutine calED(dose,ltx,sltx,ndat,ninltx,n2,inltx,&
                 outDose,mcED,pars,stdp,upb,model,origin,&
                 method,nstart,nsim,fvec1,fvalue,message)
!---------------------------------------------------------------------------
! Subroutine calED() is used for calculating  
! ED values and assess their standard errors.
!---------------------------------------------------------------------------
!        dose(ndat):: input, real values, the dose values.
!         ltx(ndat):: input, real values, Lx/Tx values
!        sltx(ndat):: input, real values, the errors of Lx/Txs.
!              ndat:: input, integer, the number of data points.
!            ninltx:: input, integer, number of ED values.
!                n2:: input, integer, dimension of the problem.
!   inltx(ninltx,2):: input, real values, the natural Lx/Txs and errors.
! outDose(ninltx,2):: output, real values, the estimated EDs and errors.
! mcED(nsim,ninltx):: output, real values, the simulated ED values.
!          pars(n2):: output, real values, pars of growth curve.
!          stdp(n2):: output, real values, the std errors of pars.
!               upb:: input, real value, upper bound of b value.
!             model:: input, integer, 0=linear, 1=exp, 2=lexp, 3=dexp.
!            origin:: input, integer, 0=origin, 1=non-origin.
!            method:: input, integer, 0=simple transformation, 1=Monte Carlo.
!            nstart:: input, integer, number of trails.
!              nsim:: input, integer, number of simulations.
!       fvec1(ndat):: output, real values, the fitted Lx/Tx values.
!            fvalue:: output, real value, the minimized chi-square.
!           message:: output, integer, 0=success, 1=fail.
!---------------------------------------------------------------------------
! Author:: Peng Jun, 2014.09.12.
!---------------------------------------------------------------------------
! Dependence:: subroutine linearFit; subroutine inipars;--------------------
!              subroutine lmgrwcve; subroutine interpolate;-----------------
!              subroutine r8vec_normal.-------------------------------------
!---------------------------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: ndat, ninltx, n2, model,&
                                  origin, method, nstart, nsim
    real   (kind=8), intent(in):: dose(ndat), ltx(ndat), sltx(ndat),&
                                  upb, inltx(ninltx,2)
    real   (kind=8), intent(out):: outDose(ninltx,2), mcED(nsim,ninltx),&
                                   pars(n2), stdp(n2), fvec1(ndat), fvalue
    integer(kind=4), intent(out):: message
    ! Local variables.
    integer(kind=4):: i, j, info, seed, mcCount, niter
    real   (kind=8):: minCond, ran(2), outp(3),&
                      locp(5), cpars(n2), cstdp(n2),&
                      cfvec1(ndat), cfvalue, cond, grad,&
                      maxDose, ivalue, avg, low(ninltx), up(ninltx),&
                      low1(ninltx), up1(ninltx), sumdose, sumdose2,& 
                      simltx(ndat), mcSig(1), mcDose
    !
    outDose = -99.0
    mcED = -99.0
    pars = -99.0
    stdp = -99.0
    fvec1 = -99.0
    fvalue = -99.0
    message = 1
    !
    if (model==0) then
        call linearFit(dose,ltx,sltx,ndat,pars,&
                       stdp,n2,fvec1,fvalue,info)
        if (info==0) message = 0
        if (message/=0) return
    else if (model==1 .or. model==2 .or. model==3) then
        minCond = 1.0D+20
        maxDose = maxval(dose)
        ! call random_seed()
        !
        loopA: do i=1, nstart
            call random_number(ran)
            ran = ran*upb
            call inipars(ran(1),ran(2),model,origin,&
                         dose,ltx,sltx,ndat,outp,info)
            if (info/=0) cycle loopA
            !
            if (info==0) then
                locp = 0.0
                if (model==1 .or. model==2) then
                    locp(1) = outp(1)
                    locp(2) = ran(1)
                    locp(3) = outp(2)
                    locp(4) = outp(3)
                else if (model==3) then
                    locp(1) = outp(1)
                    locp(2) = ran(1)
                    locp(3) = outp(2)
                    locp(4) = ran(2)
                    locp(5) = outp(3)
                end if 
                cpars = locp(1:n2)
            end if
            !
            call lmgrwcve(dose,ltx,sltx,ndat,cpars,cstdp,n2,&
                          model,origin,cfvec1,cfvalue,cond,info)
            !
            if (info/=0) cycle loopA
            !
            if (model==1 .or. model==3) then
                locp = 0.0
                locp(1:n2) = cpars
                grad = locp(1)*locp(2)*exp(-locp(2)*maxDose)+&
                       locp(3)*locp(4)*exp(-locp(4)*maxDose)
                if (grad<1.0D-13) cycle loopA
            end if
            !
            if (cond<minCond) then
                pars = cpars
                stdp = cstdp
                fvec1 = cfvec1
                fvalue = cfvalue
                minCond = cond
                message = 0
            end if
        end do loopA
        !
        if (message/=0) return
    end if
    !
    !
    ! Calculate quivalent dose values. 
    if (model==0) then
        locp = 0.0
        locp(1:n2) = pars
        outDose(:,1) = (inltx(:,1)-locp(2))/locp(1)
    else if (model==1 .or. model==2 .or. model==3) then 
        do i=1, ninltx
            call interpolate(inltx(i,1),pars,n2,model,&
                             -5.0D+00,1.3*maxDose,outDose(i,1),ivalue)
            if (ivalue>1.0D-7) then
                message = 1
                return
            end if
        end do
    end if
    !
    ! Assess the standard errors of 
    ! equivalent dose values.
    if (method==0) then
        ! Simple transformation.
        avg = sum((ltx-fvec1)**2)/real(ndat)
        low = inltx(:,1) - sqrt((inltx(:,2))**2+avg)
        up = inltx(:,1) + sqrt((inltx(:,2))**2+avg)
        !
        if (model==0) then 
            locp = 0.0
            locp(1:n2) = pars
            low1 = (low-locp(2))/locp(1)
            up1 = (up-locp(2))/locp(1)
        else if (model==1 .or. model==2 .or. model==3) then
            do i=1, ninltx
                call interpolate(low(i),pars,n2,model,&
                                 -5.0D+00,1.3*maxDose,low1(i),ivalue)
                call interpolate(up(i),pars,n2,model,&
                                 -5.0D+00,1.3*maxDose,up1(i),ivalue)
            end do
        end if
        !
        outDose(:,2) = (up1-low1)/2.0
    else if (method==1) then
        ! Monte Carlo simulation.
        seed = 123456789
        do i=1, ninltx
            mcCount = 0
            sumdose = 0.0
            sumdose2 = 0.0
            niter = 0
            Innerloop: do
                niter = niter + 1
                if (niter>100*nsim) then
                    message = 1
                    return
                end if
                !
                do j=1, ndat
                    call r8vec_normal(1,ltx(j),sltx(j),seed,simltx(j))
                end do
                !
                if (model==0) then
                    call linearFit(dose,simltx,sltx,ndat,cpars,&
                                   cstdp,n2,cfvec1,cfvalue,info)
                    if (info/=0) cycle Innerloop
                else if (model==1 .or. model==2 .or. model==3) then
                    call random_number(ran)
                    ran = ran*upb
                    call inipars(ran(1),ran(2),model,origin,dose,&
                                 simltx,sltx,ndat,outp,info)
                    if (info/=0) cycle Innerloop
                    if (info==0) then
                        locp = 0.0
                        if (model==1 .or. model==2) then
                            locp(1) = outp(1)
                            locp(2) = ran(1)
                            locp(3) = outp(2)
                            locp(4) = outp(3)
                        else if (model==3) then
                            locp(1) = outp(1)
                            locp(2) = ran(1)
                            locp(3) = outp(2)
                            locp(4) = ran(2)
                            locp(5) = outp(3)
                        end if 
                        cpars = locp(1:n2)
                    end if
                    !
                    call lmgrwcve(dose,simltx,sltx,ndat,cpars,cstdp,n2,&
                                  model,origin,cfvec1,cfvalue,cond,info)
                    if (info/=0) cycle Innerloop
                    if (model==1 .or. model==3) then
                        locp = 0.0
                        locp(1:n2) = cpars
                        grad = locp(1)*locp(2)*exp(-locp(2)*maxDose)+&
                               locp(3)*locp(4)*exp(-locp(4)*maxDose)
                        if (grad<1.0D-13) cycle Innerloop
                    end if
                end if
                !
                call r8vec_normal(1,inltx(i,1),inltx(i,2),seed,mcSig(1))
                !
                if (model==0) then
                    locp = 0.0
                    locp(1:n2) = cpars
                    mcDose = (mcSig(1)-locp(2))/locp(1)
                else if (model==1 .or. model==2 .or. model==3) then 
                    call interpolate(mcSig(1),cpars,n2,model,&
                                     -5.0D+00,1.3*maxDose,mcDose,ivalue)
                    if (ivalue>1.0D-7) cycle Innerloop
                end if
                mcCount = mcCount + 1
                sumdose = sumdose +mcDose
                sumdose2 = sumdose2 +mcDose**2
                mcED(mcCount,i) = mcDose
                if (mcCount==nsim) exit Innerloop
            end do Innerloop
            !
            outDose(i,2) = sqrt((real(nsim)*sumdose2-sumdose**2)/&
                                 real(nsim)/real(nsim-1))
        end do
    end if
    !
    return
end subroutine calED 
