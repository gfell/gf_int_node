MODULE mod_int_node
  IMPLICIT NONE

  INTERFACE sub_int_node
     MODULE PROCEDURE sub_int_node_sc, sub_int_node_ar
  END INTERFACE

  ! Set kind of real arrays
  INTEGER, PARAMETER :: ki=SELECTED_REAL_KIND(p = 15, r = 30)

CONTAINS

  SUBROUTINE sub_int_node_sc(x,xp,nodes)
    ! Returns index of right interpolating node (from grid x) for point xp

    REAL(ki), INTENT(in), DIMENSION(:) 	:: x
    REAL(ki), INTENT(in)		:: xp
    INTEGER, INTENT(out)	        :: nodes
    INTEGER, DIMENSION(1)               :: x_min
    INTEGER                             :: length_x

    length_x = SIZE(x)
    nodes = 2

    IF (xp>x(1)) THEN 
    x_min = MINLOC(ABS(x-xp)) 
    nodes = MAX(2,x_min(1))
    IF (xp>x(nodes)) nodes = MIN(length_x,nodes + 1)
    END IF
        
  END SUBROUTINE sub_int_node_sc
  
  SUBROUTINE sub_int_node_ar(x,xp,nodes)
    ! Returns array of indeces of appropriate RIGHT interpolating nodes 
    ! (from grid x) for array xp.
    !
    ! The routine exploits the fact that x and xp are ordered
    ! to speed up the location of the appropriate nodes.

    REAL (ki), DIMENSION (:), INTENT(in)	:: x
    REAL (ki), DIMENSION (:), INTENT(in)	:: xp
    INTEGER, DIMENSION (SIZE(xp)), INTENT(out)	:: nodes

    INTEGER, DIMENSION(1) :: z, i_min, i_max
    INTEGER :: j, length_x, length_xp, step, x_min, x_max
    INTEGER, DIMENSION (SIZE(xp)) :: x_index

    x_min = 1
    nodes = 2
    i_max = 1 
    length_x = SIZE(x)
    length_xp = SIZE(xp)

    ! 1. Right bracketing nodes for extrapolation regions
    !-----------------------------------------------------
    IF (xp(length_xp) < x(1)) THEN
      ! Only extrapolation to the left: nodes(:) = 2
      ! 2. is not executed
      RETURN
    ELSEIF (xp(1).GE.x(length_x-1)) THEN
      ! Only extrapolation to the right: nodes(:) = length_x 
      ! 2. is not executed
      nodes = length_x 
      RETURN
    ELSE
      ! Index of first xp for which interpolation applies
      i_min = COUNT(xp < x(1),1) + 1 
      ! Index of last xp for which interpolation applies
      i_max = COUNT(xp<x(length_x-1),1)
      IF (i_max(1)<length_xp) nodes(i_max(1)+1:) = length_x
    END IF


    ! 2. Right bracketing nodes for remaining points
    !--------------------------------------------------------
    DO j = i_min(1),i_max(1)

      ! Increment x_max until x(x_max) brackets xp(j)
      DO step = 10,length_x+10,10
        x_max = MIN(x_min+step, length_x-1)
        IF (x(x_max)>xp(j)) EXIT
      END DO

      ! Locate position of  closest node to xp(j)
      z = MINLOC(ABS(x(x_min:x_max)-xp(j))) + x_min - 1
      x_min = z(1)
      IF (xp(j)>x(x_min)) x_min = x_min + 1
      nodes(j) = x_min

    ENDDO

  END SUBROUTINE sub_int_node_ar
END MODULE mod_int_node
