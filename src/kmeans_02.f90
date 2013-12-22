subroutine kmeans_02 ( dim_num, point_num, cluster_num, it_max, it, point, &
  cluster, cluster_center, cluster_population, cluster_energy,ifault )
!
!*******************************************************************************
!
!! KMEANS_02 applies the K-Means algorithm to a partition problem.
!
!
!  Discussion:
!
!    The routine attempts to divide POINT_NUM points in 
!    DIM_NUM-dimensional space into CLUSTER_NUM clusters so that the within 
!    cluster sum of squares is minimized.
!
!  Reference:
!
!    J A Hartigan, M A Wong,
!    Algorithm AS 136: A K-Means Clustering Algorithm,
!    Applied Statistics,
!    Volume 28, Number 1, 1979, pages 100-108.
!
!  Modified:
!
!    10 March 2002
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, integer CLUSTER_NUM, the number of clusters.
!
!    Input, integer IT_MAX, the maximum number of iterations allowed.
!
!    Output, integer IT, the number of iterations taken.
!
!    Input, real POINT(DIM_NUM,POINT_NUM), the coordinates of the points.
!
!    Output, integer CLUSTER(POINT_NUM), the cluster each point belongs to.
!
!    Input/output, real CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the cluster
!    centers.
!
!    Output, integer CLUSTER_POPULATION(CLUSTER_NUM), the number of points 
!    in each cluster.
!
!    Output, real CLUSTER_ENERGY(CLUSTER_NUM), the within-cluster sum 
!    of squares.
!
  implicit none
!
  integer(kind=4) cluster_num
  integer(kind=4) dim_num
  integer(kind=4) point_num
!
  real(kind=8) an1(cluster_num)
  real(kind=8) an2(cluster_num)
  integer(kind=4) cluster(point_num)
  integer(kind=4) cluster2(point_num)
  real(kind=8) cluster_center(dim_num,cluster_num)
  real(kind=8) cluster_energy(cluster_num)
  integer(kind=4) cluster_population(cluster_num)
  real(kind=8) d(point_num)
  real(kind=8) dc
  real(kind=8) db
  real(kind=8) dt(2)
  integer(kind=4) i 
  integer(kind=4) ifault
  integer(kind=4) il
  integer(kind=4) indx
  integer(kind=4) it
  integer(kind=4) it2
  integer(kind=4) it_max
  integer(kind=4) itran(cluster_num)
  integer(kind=4) j
  integer(kind=4) l
  integer(kind=4) live(cluster_num)
  integer(kind=4) ncp(cluster_num)
  real(kind=8) point(dim_num,point_num)
  real(kind=8) temp
!
  it = 0
!
!  For each point I, find its two closest centers, CLUSTER(I) and CLUSTER2(I).
!  Assign it to CLUSTER(I).
!
  do i = 1, point_num

    cluster(i) = 1
    cluster2(i) = 2

    do il = 1, 2
      dt(il) = sum ( &
        ( point(1:dim_num,i) - cluster_center(1:dim_num,il) )**2 )
    end do

    if ( dt(1) > dt(2) ) then
      cluster(i) = 2
      cluster2(i) = 1
      temp = dt(1)
      dt(1) = dt(2)
      dt(2) = temp
    end if

    do l = 3, cluster_num

      db = sum ( ( point(1:dim_num,i) - cluster_center(1:dim_num,l) )**2 )

      if ( db < dt(1) ) then
        dt(2) = dt(1)
        cluster2(i) = cluster(i)
        dt(1) = db
        cluster(i) = l
      else if ( db < dt(2) ) then
        dt(2) = db
        cluster2(i) = l
      end if

    end do

  end do
!
!  Update cluster centers to be the average of points contained
!  within them.
!
  cluster_population(1:cluster_num) = 0
  cluster_center(1:dim_num,1:cluster_num) = 0.0D+00

  do i = 1, point_num
    l = cluster(i)
    cluster_population(l) = cluster_population(l) + 1
    cluster_center(1:dim_num,l) = cluster_center(1:dim_num,l) &
      + point(1:dim_num,i)
  end do
!
!  Check to see if there is any empty cluster.
!
  do l = 1, cluster_num
    cluster_center(1:dim_num,l) = cluster_center(1:dim_num,l) &
      / real ( cluster_population(l) )
!
!  Initialize AN1, AN2, ITRAN and NCP
!  AN1(L) = CLUSTER_POPULATION(L) / (CLUSTER_POPULATION(L) - 1)
!  AN2(L) = CLUSTER_POPULATION(L) / (CLUSTER_POPULATION(L) + 1)
!  ITRAN(L) = 1 if cluster L is updated in the quick-transfer stage,
!           = 0 otherwise
!  In the optimal-transfer stage, NCP(L) stores the step at which
!  cluster L is last updated.
!  In the quick-transfer stage, NCP(L) stores the step at which
!  cluster L is last updated plus POINT_NUM.
!
    an2(l) = real ( cluster_population(l) ) / real ( cluster_population(l) + 1 )

    if ( cluster_population(l) > 1 ) then
      an1(l) = real ( cluster_population(l) ) &
        / real ( cluster_population(l) - 1 )
    else
      an1(l) = huge ( an1(l) )
    end if

    itran(l) = 1
    ncp(l) = -1

  end do

  indx = 0
  ifault = 2

  do it2 = 1, it_max

    it = it2
!
!  In this stage, there is only one pass through the data.   Each
!  point is re-allocated, if necessary, to the cluster that will
!  induce the maximum reduction in within-cluster sum of squares.
!
    call kmeans_02_optra ( dim_num, point_num, cluster_num, point, &
      cluster_center, cluster, cluster2, cluster_population, an1, an2, &
      ncp, d, itran, live, indx )
!
!  Stop if no transfer took place in the last POINT_NUM optimal transfer steps.
!
    if ( indx == point_num ) then
      ifault = 0
      exit
    end if
!
!  Each point is tested in turn to see if it should be re-allocated
!  to the cluster to which it is most likely to be transferred,
!  CLUSTER2(I), from its present cluster, CLUSTER(I).   Loop through the
!  data until no further change is to take place.
!
    call kmeans_02_qtran ( dim_num, point_num, cluster_num, point, &
      cluster_center, cluster, cluster2, cluster_population, an1, an2, &
      ncp, d, itran, indx )
!
!  If there are only two clusters, there is no need to re-enter the
!  optimal transfer stage.
!
    if ( cluster_num == 2 ) then
      ifault = 0
      exit
    end if
!
!  NCP has to be set to 0 before entering OPTRA.
!
    ncp(1:cluster_num) = 0

  end do
!
!  Compute the within-cluster sum of squares for each cluster.
!
  cluster_center(1:dim_num,1:cluster_num) = 0.0D+00

  do i = 1, point_num
    cluster_center(1:dim_num,cluster(i)) = &
    cluster_center(1:dim_num,cluster(i)) + point(1:dim_num,i)
  end do

  do j = 1, dim_num
    cluster_center(j,1:cluster_num) = cluster_center(j,1:cluster_num) &
      / real ( cluster_population(1:cluster_num) )
  end do

  call cluster_energy_compute ( dim_num, point_num, cluster_num, point, &
    cluster, cluster_center, cluster_energy )

  return
end
subroutine kmeans_02_optra ( dim_num, point_num, cluster_num, point, &
  cluster_center, cluster, cluster2, cluster_population, an1, an2, ncp, &
  d, itran, live, indx )
!
!*******************************************************************************
!
!! KMEANS_02_OPTRA carries out the optimal transfer stage.
!
!
!  Discussion:
!
!    Each point is re-allocated, if necessary, to the cluster that
!    will induce a maximum reduction in the within-cluster sum of
!    squares.
!
!  Reference:
!
!    J A Hartigan, M A Wong,
!    A K-Means Clustering Algorithm,
!    Applied Statistics,
!    Volume 28, Number 1, 1979, pages 100-108.
!
!  Modified:
!
!    08 March 2002
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, integer CLUSTER_NUM, the number of clusters.
!
!    Input, real POINT(DIM_NUM,POINT_NUM), the coordinates of the points.
!
!    Input/output, real CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the cluster
!    centers.
!
!    Input/output, integer CLUSTER(POINT_NUM), the cluster each point 
!    belongs to.
!
!    Input/output, integer CLUSTER2(POINT_NUM), the cluster to which
!    each point is most likely to be transferred to.
!
!    Input/output, integer CLUSTER_POPULATION(CLUSTER_NUM), the number of 
!    points in each cluster.
!
!    ?, real AN1(CLUSTER_NUM), ?
!
!    ?, real AN2(CLUSTER_NUM), ?
!
!    ?, integer NCP(CLUSTER_NUM), ?
!
!    ?, real D(POINT_NUM), ?
!
!    Input/output, integer ITRAN(CLUSTER_NUM), ?
!
!    Input/output, integer LIVE(CLUSTER_NUM), ?
!
!    Input/output, integer INDX, ?
!
  implicit none
!
  integer(kind=4) cluster_num
  integer(kind=4) dim_num
  integer(kind=4) point_num
!
  real(kind=8) al1
  real(kind=8) al2
  real(kind=8) alt
  real(kind=8) alw
  real(kind=8) an1(cluster_num)
  real(kind=8) an2(cluster_num)
  integer(kind=4) cluster(point_num)
  integer(kind=4) cluster2(point_num)
  real(kind=8) cluster_center(dim_num,cluster_num)
  integer(kind=4) cluster_population(cluster_num)
  real(kind=8) d(point_num)
  real(kind=8) dc
  real(kind=8) dd
  integer(kind=4) i
  integer(kind=4) indx
  integer(kind=4) itran(cluster_num)
  integer(kind=4) j
  integer(kind=4) l
  integer(kind=4) l1
  integer(kind=4) l2
  integer(kind=4) live(cluster_num)
  integer(kind=4) ll
  integer(kind=4) ncp(cluster_num)
  real(kind=8) point(dim_num,point_num)
  real(kind=8) r2
  real(kind=8) rr
!
!  If cluster L is updated in the last quick-transfer stage, it
!  belongs to the live set throughout this stage.   Otherwise, at
!  each step, it is not in the live set if it has not been updated
!  in the last POINT_NUM optimal transfer steps.
!
  do l = 1, cluster_num
    if ( itran(l) == 1 ) then
      live(l) = point_num + 1
    end if
  end do

  do i = 1, point_num

    indx = indx + 1
    l1 = cluster(i)
    l2 = cluster2(i)
    ll = l2
!
!  If point I is the only member of cluster L1, no transfer.
!
    if ( cluster_population(l1) > 1 ) then
!
!  If L1 has been updated in this stage, re-compute D(I).
!
      if ( ncp(l1) /= 0 ) then
        d(i) = an1(l1) * sum ( &
          ( point(1:dim_num,i) - cluster_center(1:dim_num,l1) )**2 )
      end if
!
!  Find the cluster with minimum R2.
!
      r2 = an2(l2) * sum ( &
        ( point(1:dim_num,i) - cluster_center(1:dim_num,l2) )**2 )

      do l = 1, cluster_num
!
!  If I >= LIVE(L1), then L1 is not in the live set.   If this is
!  true, we only need to consider clusters that are in the live set
!  for possible transfer of point I.   
!
!  Otherwise, we need to consider all possible clusters.
!
!  BADLY FORMED LOGICAL HERE.
!
        if ( i >= live(l1) .and. i >= live(l) .or. l == l1 .or. l == ll ) then
          go to 60
        end if

        rr = r2 / an2(l)

        dc = sum ( ( point(1:dim_num,i) - cluster_center(1:dim_num,l) )**2 )

        if ( dc < rr ) then
          r2 = dc * an2(l)
          l2 = l
        end if

60      continue

      end do
!
!  If no transfer is necessary, L2 is the new CLUSTER2(I).
! 
      if ( r2 >= d(i) ) then

        cluster2(i) = l2

      else
!
!  Update cluster centers, LIVE, NCP, AN1 and AN2 for clusters L1 and
!  L2, and update CLUSTER(I) and CLUSTER2(I).
!
        indx = 0
        live(l1) = point_num + i
        live(l2) = point_num + i
        ncp(l1) = i
        ncp(l2) = i
        al1 = cluster_population(l1)
        alw = al1 - 1.0D+00
        al2 = cluster_population(l2)
        alt = al2 + 1.0D+00

        cluster_center(1:dim_num,l1) = ( cluster_center(1:dim_num,l1) * al1 &
          - point(1:dim_num,i) ) / alw

        cluster_center(1:dim_num,l2) = ( cluster_center(1:dim_num,l2) * al2 &
          + point(1:dim_num,i) ) / alt

        cluster_population(l1) = cluster_population(l1) - 1
        cluster_population(l2) = cluster_population(l2) + 1
        an2(l1) = alw / al1

        if ( alw > 1.0D+00 ) then
          an1(l1) = alw / ( alw - 1.0D+00 )
        else
          an1(l1) = huge ( an1(l1) )
        end if

        an1(l2) = alt / al2
        an2(l2) = alt / ( alt + 1.0D+00 )
        cluster(i) = l2
        cluster2(i) = l1

      end if

    end if

    if ( indx == point_num ) then
      return
    end if

  end do
!
!  ITRAN(L) = 0 before entering QTRAN.
!
  itran(1:cluster_num) = 0
!
!  LIVE(L) has to be decreased by POINT_NUM before re-entering OPTRA.
!
  live(1:cluster_num) = live(1:cluster_num) - point_num

  return
end
subroutine kmeans_02_qtran ( dim_num, point_num, cluster_num, point, &
  cluster_center, cluster, cluster2, cluster_population, an1, an2, ncp, &
  d, itran, indx )
!
!*******************************************************************************
!
!! KMEANS_02_QTRAN carries out the quick transfer stage.
!
!
!  Discussion:
!
!    For each point I, CLUSTER(I) and CLUSTER2(I) are switched, if necessary, 
!    to reduce within-cluster sum of squares.  The cluster centers are
!    updated after each step.
!
!  Reference:
!
!    J A Hartigan, M A Wong,
!    A K-Means Clustering Algorithm,
!    Applied Statistics,
!    Volume 28, Number 1, 1979, pages 100-108.
!
!  Modified:
!
!    10 March 2002
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, integer CLUSTER_NUM, the number of clusters.
!
!    Input, real POINT(DIM_NUM,POINT_NUM), the coordinates of the points.
!
!    Input/output, real CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the cluster
!    centers.
!
!    Input/output, integer CLUSTER(POINT_NUM), the cluster each point 
!    belongs to.
!
!    Input/output, integer CLUSTER2(POINT_NUM), the cluster to which
!    each point is most likely to be transferred to.
!
!    Input/output, integer CLUSTER_POPULATION(CLUSTER_NUM), the number of 
!    points in each cluster.
!
!    ?, real AN1(CLUSTER_NUM), ?
!
!    ?, real AN2(CLUSTER_NUM), ?
!
!    ?, integer NCP(CLUSTER_NUM), ?
!
!    ?, real D(POINT_NUM), ?
!
!    ?, integer ITRAN(CLUSTER_NUM), ?
!
!    Input/output, integer INDX, is set to 0 if any updating occurs.
!
  implicit none
!
  integer(kind=4) cluster_num
  integer(kind=4) dim_num
  integer(kind=4) point_num
!
  real(kind=8) al1
  real(kind=8) al2
  real(kind=8) alt
  real(kind=8) alw
  real(kind=8) an1(cluster_num)
  real(kind=8) an2(cluster_num)
  integer(kind=4) cluster(point_num)
  integer(kind=4) cluster2(point_num)
  real(kind=8) cluster_center(dim_num,cluster_num)
  integer(kind=4) cluster_population(cluster_num)
  integer(kind=4) count
  real(kind=8) d(point_num)
  real(kind=8) dd
  integer(kind=4) i
  integer(kind=4) indx
  integer(kind=4) itran(cluster_num)
  integer(kind=4) j
  integer(kind=4) l1
  integer(kind=4) l2
  integer(kind=4) ncp(cluster_num)
  real(kind=8) point(dim_num,point_num)
  real(kind=8) r2
  integer(kind=4) step
!
!  In the optimal transfer stage, NCP(L) indicates the step at which
!  cluster L is last updated.   In the quick transfer stage, NCP(L)
!  is equal to the step at which cluster L is last updated plus POINT_NUM.
!
  count = 0
  step = 0

  do

    do i = 1, point_num

      count = count + 1
      step = step + 1
      l1 = cluster(i)
      l2 = cluster2(i)
!
!  If point I is the only member of cluster L1, no transfer.
!
      if ( cluster_population(l1) > 1 ) then
!
!  If STEP > NCP(L1), no need to re-compute distance from point I to
!  cluster L1.   Note that if cluster L1 is last updated exactly POINT_NUM
!  steps ago, we still need to compute the distance from point I to
!  cluster L1.
!
        if ( step <= ncp(l1) ) then
          d(i) = an1(l1) * sum ( &
            ( point(1:dim_num,i) - cluster_center(1:dim_num,l1) )**2 )
        end if
!
!  If STEP >= both NCP(L1) and NCP(L2) there will be no transfer of
!  point I at this step.
!
        if ( step < ncp(l1) .or. step < ncp(l2) ) then

          r2 = d(i) / an2(l2)

          dd = sum ( ( point(1:dim_num,i) - cluster_center(1:dim_num,l2) )**2 )
!
!  Update cluster centers, NCP, CLUSTER_POPULATION, ITRAN, AN1 and AN2 
!  for clusters L1 and L2.   Also update CLUSTER(I) and CLUSTER2(I).   
!
!  Note that if any updating occurs in this stage, INDX is set back to 0.
!
          if ( dd < r2 ) then

            count = 0
            indx = 0
            itran(l1) = 1
            itran(l2) = 1
            ncp(l1) = step + point_num
            ncp(l2) = step + point_num
            al1 = cluster_population(l1)
            alw = al1 - 1.0D+00
            al2 = cluster_population(l2)
            alt = al2 + 1.0D+00

            cluster_center(1:dim_num,l1) = &
              ( cluster_center(1:dim_num,l1) * al1 &
              - point(1:dim_num,i) ) / alw

            cluster_center(1:dim_num,l2) = &
              ( cluster_center(1:dim_num,l2) * al2 &
              + point(1:dim_num,i) ) / alt

            cluster_population(l1) = cluster_population(l1) - 1
            cluster_population(l2) = cluster_population(l2) + 1
            an2(l1) = alw / al1

            if ( alw > 1.0D+00 ) then
              an1(l1) = alw / ( alw - 1.0D+00 )
            else
              an1(l1) = huge ( an1(l1) )
            end if

            an1(l2) = alt / al2
            an2(l2) = alt / ( alt + 1.0D+00 )
            cluster(i) = l2
            cluster2(i) = l1

          end if

        end if

      end if
!
!  If no re-allocation took place in the last POINT_NUM steps, return.
!
      if ( count == point_num ) then
        return
      end if

    end do

  end do

  return
end
subroutine cluster_energy_compute ( dim_num, point_num, cluster_num, point, &
  cluster, cluster_center, cluster_energy )
!
!*******************************************************************************
!
!! CLUSTER_ENERGY_COMPUTE computes the energy of the clusters.
!
!
!  Modified:
!
!    18 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer POINT_NUM, the number of data points.
!
!    Input, integer CLUSTER_NUM, the number of clusters.
!
!    Input, real POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Output, integer CLUSTER(POINT_NUM), the cluster to which each
!    data point belongs.
!
!    Input/output, real CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the 
!    centers associated with the minimal energy clustering.
!
!    Output, real CLUSTER_ENERGY(CLUSTER_NUM), the energy associated with 
!    each cluster.
!
  implicit none
!
  integer(kind=4) cluster_num
  integer(kind=4) dim_num
  integer(kind=4) point_num
!
  integer(kind=4), dimension ( point_num ) :: cluster
  real(kind=8), dimension ( dim_num, cluster_num ) :: cluster_center
  real(kind=8), dimension ( cluster_num ) :: cluster_energy
  integer(kind=4) i
  integer(kind=4) j
  real(kind=8), dimension ( dim_num, point_num ) :: point
  real(kind=8) point_energy
!
  cluster_energy(1:cluster_num) = 0.0D+00

  do i = 1, point_num

    j = cluster(i)

    point_energy = sum ( &
      ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )

    cluster_energy(j) = cluster_energy(j) + point_energy

  end do

  return
end
