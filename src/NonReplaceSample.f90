subroutine NonReplaceSample(TotalSample,WantedSample,TotalSampleSize,WantedSampleSize)
!-----------------------------------------------------------------------------------------
! NonReplaceSample() is a subroutine used to do sampling without replace employing 
! the mothod of Bert and Green (1977).
! ========================================================================================
! TotalSample::       Input, the total samples which a subsample generate from.
!
! WantedSample::      Output, the wanted subsample which generates from the TotalSample.
!
! TotalSampleSize::   Input, the size of the TotalSample.
!
! WantedSampleSize::  Input, the size of the WantedSample.
! ========================================================================================
! Author:: Peng Jun, 2012.12.26.
!
! Reference:: Bert F and Green JR, 1977. Fortran subroutines for random sampling without 
!             replacement. Behavior Research Methodt and Instrumentation, pp 559.
!
! Dependence:: No
!-----------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::TotalSampleSize
  integer(kind=4),intent(in)::WantedSampleSize
  integer(kind=4),dimension(TotalSampleSize),intent(in)::TotalSample
  integer(kind=4),dimension(WantedSampleSize),intent(out)::WantedSample
  !local variables
  integer(kind=4)::M
  integer(kind=4)::i
  integer(kind=4)::L
  real(kind=8),dimension(1)::RanDom
  !
  M=0
  !
  do i=1,TotalSampleSize
    !
    call random_number(RanDom)
    !
    L=int((float(TotalSampleSize-i+1))*RanDom(1))+1
    !
    if(L.GT.(WantedSampleSize-M)) cycle
    !
    M=M+1
    WantedSample(M)=TotalSample(i)
    !
    if(M.GE.WantedSampleSize) then
      return
    end if 
    !
  end do
  return
end subroutine NonReplaceSample
