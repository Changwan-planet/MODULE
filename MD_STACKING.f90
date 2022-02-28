MODULE MD_STACKING
IMPLICIT NONE

CONTAINS

  SUBROUTINE stacking(removed_bgr_signal,ROWS,DIS,stacked_signal)
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

  END SUBROUTINE stacking

END MODULE MD_STACKING


