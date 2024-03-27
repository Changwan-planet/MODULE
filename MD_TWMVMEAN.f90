!PRINT_DATE: Thu 14 Apr 2022 04:26:54 PM KST
MODULE MD_TWMVMEAN
IMPLICIT NONE

CONTAINS
  SUBROUTINE twmvmean(B_SCAN_IMAGE4, DIS, TRA, ROWS, WW, B_SCAN_IMAGE5) 
  IMPLICIT NONE

  INTEGER :: X,Y,Z 
  INTEGER, PARAMETER :: MV_DIS = (DIS - WW) + 1
  INTEGER, PARAMETER :: MV_TRA = (TRA - WW) +1
  INTEGER :: MV_X, MV_Y

  INTEGER, INTENT(IN) :: DIS
  INTEGER, INTENT(IN) :: TRA
  INTEGER, INTENT(IN) :: ROWS
  INTEGER, INTENT(IN) :: WW !the number of values for moving average
  
  REAL*8, DIMENSION(DIS,TRA,ROWS), INTENT(IN) :: B_SCAN_IMAGE4
  REAL*8, DIMENSION(:,:,:), ALLOCATABLE ::B_SCAN_IMAGE5
  
  REAL*8, DIMENSION(1,1,ROWS) :: TEMP
  REAL*8 :: SUM_2DW_amp

    DO Z = 1, ROWS
       DO MV_X = 1, MV_DIS
       DO MV_Y = 1, MV_TRA

          SUM_2DW_amp = 0.0
          DO X = MV_X, MV_X+WW-1
          DO Y = MV_Y, MV_Y+WW-1

             SUM_2DW_amp = SUM_2DW_amp + B_SCAN_IMAGE4(X,Y,Z)
                                             
          END DO
          END DO

          B_SCAN_IMAGE5(MV_X,MV_Y,Z) = ( SUM_2DW_amp ) / (W**2)

       END DO
       END DO
    END DO 

  END SUBROUTINE twmvmean

END MODULE MD_TWMVMEAN
