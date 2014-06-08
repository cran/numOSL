subroutine kmeans(ndim,npoi,nclu,maxiter,nstart, &
                  points,belo,clust,clusp,energ,ifault)
!-----------------------------------------------------------------------------------
! kmeans() is a wrapper subroutine calling subroutine kmeans_02 to perform 
! K-Means clustering using the algorithm of Hartigan and Wong (1979)
! subroutine kmeans allows choosing some random initial sets to perform
! the K-Means algorithm by specifying nstart and return with the cluster
! that gives the smallest cluster energy.
! ==================================================================================
! 
! ndim               :: input, integer, the dimension of the data matrix.
!
! npoi               :: input, integer, the number of the data points.
!
! nclu               :: input, integer, the specified number of clusters.
!
! maxiter            :: input, integer, the maximum iterative number allowed.
!
! nstart             :: input, integer, the number of random sets used.
!
! points(ndim,npoi)  :: input, real values, the data points.
!
! belo(npoi)         :: output, integer values, which cluster each point belongs to.
!
! clust(ndim,nclu)   :: output, real values, the estimated cluster centers.
!
! clusp(nclu)        :: output, integer values, how many points in each cluster.
!
! energ(nclu)        :: output, real values, the calculated energy of each cluster.
!
! ifault             :: output, integer, 0 for a successful work, otherwise 2.
! ==================================================================================
! Author:: Peng Jun, 2013.03.14; revised in 2014.06.08.
! 
! References:: Hartigan JA. and Wong MA, (1979).  A K-means clustering
!              algorithm. Applied Statistics V28, pp.100-108.
!
! Dependence:: subroutine kmeans_02, subroutine NonReplaceSample.
!
!-----------------------------------------------------------------------------------
  implicit none
  integer(kind=4),                     intent(in)::ndim
  integer(kind=4),                     intent(in)::npoi
  integer(kind=4),                     intent(in)::nclu
  integer(kind=4),                     intent(in)::maxiter
  integer(kind=4),                     intent(in)::nstart
  integer(kind=4),                     intent(out)::ifault
  integer(kind=4),dimension(npoi),     intent(out)::belo
  integer(kind=4),dimension(nclu),     intent(out)::clusp
  real   (kind=8),dimension(nclu),     intent(out)::energ
  real   (kind=8),dimension(ndim,npoi),intent(in)::points
  real   (kind=8),dimension(ndim,nclu),intent(out)::clust
  ! local variables
  integer(kind=4)::i,it
  real   (kind=8)::maxNum
  integer(kind=4),dimension(npoi)::TotalSample
  integer(kind=4),dimension(nclu)::WantedSample
  integer(kind=4),dimension(npoi)::sbelo
  integer(kind=4),dimension(nclu)::sclusp
  real   (kind=8),dimension(nclu)::senerg
  real   (kind=8),dimension(ndim,nclu)::sclust
  !
  maxNum=1.0D+30
  !
  do i=1, npoi
    TotalSample(i)=i
  end do
  !
  call random_seed()
  !
  do i=1, nstart
    !
    call NonReplaceSample(TotalSample,WantedSample,npoi,nclu)
    !
    sclust=points(:,WantedSample)
    !
    call kmeans_02 ( ndim, npoi, nclu, maxiter, it, points, &
                     sbelo, sclust, sclusp, senerg,ifault )
    !
    if (sum(senerg)<maxNum)  then
      maxNum=sum(senerg)
      belo=sbelo
      clusp=sclusp
      clust=sclust
      energ=senerg
    end if
    !
  end do
  !
  return
end subroutine kmeans
