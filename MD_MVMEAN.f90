!PRINT_DATE: Thu 14 Apr 2022 04:26:54 PM KST
MODULE MD_MVMEAN
IMPLICIT NONE

CONTAINS
  SUBROUTINE mv_meanx(B_SCAN_IMAGE, MV, DIS, ROWS, MV_MEAN_BSCAN) 
  IMPLICIT NONE

  INTEGER :: Z !depth
  INTEGER :: X, XX
  INTEGER :: MV_X

  INTEGER, INTENT(IN) :: DIS
  INTEGER, INTENT(IN) :: ROWS
  INTEGER, INTENT(IN) :: MV !the number of values for moving average
  
  REAL*8, DIMENSION(DIS,1,ROWS), INTENT(IN) :: B_SCAN_IMAGE
  REAL*8, DIMENSION(:,:,:), ALLOCATABLE ::MV_MEAN_BSCAN
  
  REAL*8, DIMENSION(1,1,ROWS) :: TEMP

!    ALLOCATE(MV_MEAN_BSCAN(DIS3,1,ROWS)) 
  
    MV_X = 1

    DO X = 1, DIS, MV   
     
       TEMP(1,1,:) = 0.0
    
    DO Z = 1, ROWS
        TEMP(1,1,Z) = SUM(B_SCAN_IMAGE(X:X+MV-1,1,Z)) / REAL(MV)
    END DO
         MV_MEAN_BSCAN(MV_X,1,:) = TEMP(1,1,:)
         !PRINT *, "MV_X=",MV_X
         MV_X = MV_X + 1
    END DO
   
   
! DEALLOCATE(MV_MEAN_BSCAN)
 
  END SUBROUTINE mv_meanx


  SUBROUTINE mv_meanz(B_SCAN_IMAGE, MV_WIN, DIS, ROWS, ROWS2, MV_MEAN_BSCAN) 
  IMPLICIT NONE
  
  INTEGER :: Z !depth
  INTEGER :: X, XX
  INTEGER, INTENT(IN) :: MV_WIN
  
  INTEGER, INTENT(IN) :: DIS
  INTEGER, INTENT(IN) :: ROWS
  INTEGER, INTENT(IN) :: ROWS2
  
  REAL*8, DIMENSION(DIS,1,ROWS), INTENT(IN) :: B_SCAN_IMAGE
  REAL*8, DIMENSION(DIS,1,ROWS2) :: MV_MEAN_BSCAN
  REAL*8, DIMENSION(1,1,ROWS2) :: TEMP

!    PRINT*, MV_WIN, DIS, ROWS, ROWS2
    
    DO X = 1, DIS  
     
       TEMP(1,1,:) = 0.0
    
    DO Z = 1, ROWS2
         TEMP(1,1,Z) = SUM(B_SCAN_IMAGE(X,1,Z:Z+(MV_WIN-1))) / REAL(MV_WIN)
    END DO
         MV_MEAN_BSCAN(X,1,:) = TEMP(1,1,:)
    END DO
  
  END SUBROUTINE mv_meanz




  SUBROUTINE meany(B_SCAN_IMAGE, DIS, ROWS, A, B, C, TEMP2) 
  IMPLICIT NONE
  
  INTEGER :: Z !depth
  INTEGER :: X, Y 
  
  INTEGER, INTENT(IN) :: DIS
  INTEGER, INTENT(IN) :: ROWS
  INTEGER, INTENT(IN) :: A, B, C
 

  REAL*8, DIMENSION(DIS,C,ROWS), INTENT(IN) :: B_SCAN_IMAGE
  REAL*8, DIMENSION(DIS,1,ROWS) :: TEMP, TEMP2
  

    TEMP(:,1,:) = 0.0
    
    DO Y = A, B
         TEMP(:,1,:) = TEMP(:,1,:) + B_SCAN_IMAGE(:,Y,:)
    END DO
         
    TEMP2(:,1,:) = TEMP(:,1,:) / (B-A+1)

 
  END SUBROUTINE meany




  SUBROUTINE twmvmean(B_SCAN_IMAGE4, DIS, TRA, ROWS, W, B_SCAN_IMAGE5) 
  IMPLICIT NONE
  
  INTEGER :: X,Y,Z 
  INTEGER, INTENT(IN) :: W !the number of values for moving average
  INTEGER, INTENT(IN) :: DIS
  INTEGER, INTENT(IN) :: TRA
  INTEGER, INTENT(IN) :: ROWS
  
  INTEGER :: MV_DIS
  INTEGER :: MV_TRA
  INTEGER :: MV_X, MV_Y
  
  
  REAL*8, DIMENSION(DIS,TRA,ROWS), INTENT(IN) :: B_SCAN_IMAGE4
  REAL*8, DIMENSION(:,:,:), ALLOCATABLE ::B_SCAN_IMAGE5
  
  REAL*8, DIMENSION(1,1,ROWS) :: TEMP
  REAL*8 :: SUM_2DW_amp
  
  MV_DIS = (DIS - W) + 1
  MV_TRA = (TRA - W) +1
 
    DO Z = 1, ROWS
       DO MV_X = 1, MV_DIS
       DO MV_Y = 1, MV_TRA

          SUM_2DW_amp = 0.0
          DO X = MV_X, MV_X+W-1
          DO Y = MV_Y, MV_Y+W-1

             SUM_2DW_amp = SUM_2DW_amp + B_SCAN_IMAGE4(X,Y,Z)
                                             
          END DO
          END DO

          B_SCAN_IMAGE5(MV_X,MV_Y,Z) = ( SUM_2DW_amp ) / (W**2)

       END DO
       END DO
    END DO 

  END SUBROUTINE twmvmean

END MODULE MD_MVMEAN
