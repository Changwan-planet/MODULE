MODULE MD_STACKING
IMPLICIT NONE

CONTAINS

  SUBROUTINE static_stacking(removed_bgr_signal,ROWS,DIS,stacked_signal)
  IMPLICIT NONE

  REAL*8, DIMENSION(DIS,1:ROWS), INTENT(IN) :: removed_bgr_signal
  INTEGER, INTENT(IN) :: ROWS
  INTEGER, INTENT(IN) :: DIS
  REAL*8, DIMENSION(1,1:ROWS) :: stacked_signal
  
  INTEGER :: I

  !Start from the 1 index. 
  !This is becasue I did not fix the number in the header.
  
  !CALCULATE THE AVERGE OF EACH ROW.

  !PRINT *, "DIS=",DIS

  DO I = 1,ROWS
    stacked_signal(1,I) = SUM(removed_bgr_signal(:,I)) / DIS
  END DO

  END SUBROUTINE static_stacking

  SUBROUTINE lateral_stacking(mv_bgr_signal,DIS3,TRA,ROWS,stacked_signal)
  IMPLICIT NONE
  
  INTEGER :: X, Y, Z 

  REAL*8, DIMENSION(DIS3,TRA,ROWS),INTENT(IN) :: mv_bgr_signal
  INTEGER, INTENT(IN) :: ROWS
  INTEGER, INTENT(IN) :: TRA
  INTEGER, INTENT(IN) :: DIS3

  REAL*8, DIMENSION(DIS3,1,ROWS) :: stacked_signal
  

  !Start from the 1 index. 
  !This is becasue I did not fix the number in the header.
  
  !CALCULATE THE AVERGE OF EACH ROW.


!  PRINT *, "DIS3=",DIS3, 'TRA=',TRA, 'ROWS=',ROWS
  
  DO X = 1, DIS3
  DO Z = 1, ROWS
 
     stacked_signal(X,Y,Z) = SUM(mv_bgr_signal(X,:,Z)) / TRA
  
  END DO
  END DO 

  END SUBROUTINE lateral_stacking

END MODULE MD_STACKING


