  SUBROUTINE interpolate_index(ndim1, ndim2, g1, g2, k2first, k2last, k2index)
    
    ! This subroutine initializes the module variables which
    ! are used by the interpolation subroutines.
    
    ! Input:
    ! ------
    ! ndim1   : dimension of grid 1
    ! g1      : coordinates of grid 1
    ! ndim2   : dimension of grid 2
    ! g2      : coordinates of grid 2
    
    ! Output:
    ! -------
    ! k2first : index of first coordinate of grid 2
    !           in the range of grid 1
    ! k2last  : index of last coordinate of grid 2
    !           in the range of grid 1
    ! k2index : contains for each coordinate of grid2
    !           in the range of grid1 the index of the
    !           nearest bigger coordinate of grid1
    
    ! Data on grid 1 may be interpolated to the index
    ! range [k2first,k2last] of grid 2.
    
    ! Authors: Marco Giorgetta, MPI for Meteorology, Hamburg
    !          October 1999
    !          Martin Schraner, ETH Zurich, February 2009
    
    IMPLICIT NONE
    
    ! Subroutine arguments:
    INTEGER, INTENT(in)                :: ndim1, ndim2
    REAL, INTENT(in)                   :: g1(ndim1), g2(ndim2)
    INTEGER, INTENT(out)               :: k2first, k2last
    INTEGER, INTENT(out)               :: k2index(ndim2)
   
    ! Local variables:
    INTEGER  :: k1, k2

    ! Executable statements:
    
    ! Find index of first element of grid 2
    ! which is in the range of grid 1
    ! -------------------------------------
    k2first=0
    k2=1

    DO
       IF (g2(k2) .GE. g1(1)) THEN
          k2first=k2
       ELSE
          k2=k2+1
       ENDIF
       IF (k2first .NE. 0) EXIT
    ENDDO

    ! Find index of last element of grid 2
    ! which is in the range of grid 1
    ! ------------------------------------
    k2last=0
    k2=ndim2
    DO
       IF (g2(k2) .LE. g1(ndim1)) THEN
          k2last=k2
       ELSE
          k2=k2-1
       ENDIF
       IF (k2last .NE. 0) EXIT
    ENDDO
    
    ! Find indices k2index for elements k2first
    ! to k2last of grid 2
    ! -----------------------------------------
    k2index(:)=0
    k2=k2first
    k1=2
    DO
       IF (g2(k2) .LE. g1(k1)) THEN
          k2index(k2)=k1
          k2=k2+1
       ELSE
          k1=k1+1
       ENDIF
       IF (k2 .GT. k2last) EXIT
    ENDDO
    
  END SUBROUTINE interpolate_index
