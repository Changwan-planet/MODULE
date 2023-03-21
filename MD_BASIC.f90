!PRINT_DATE: Thu 14 Apr 2022 04:26:54 PM KST
MODULE MD_BASIC
IMPLICIT NONE

CONTAINS
  SUBROUTINE f_symmetry(SIGNAL, ROWS, SMTR_SIGNAL) 
  IMPLICIT NONE

    INTEGER :: Z
    INTEGER, INTENT(IN) :: ROWS
    REAL*8, DIMENSION(ROWS), INTENT(IN) :: SIGNAL
    REAL*8, DIMENSION(ROWS) :: SMTR_SIGNAL

    DO Z = 1, ROWS/2
          SMTR_SIGNAL(Z) = SIGNAL(Z+ROWS/2)
    END DO

    DO Z = ROWS/2+1,ROWS
          SMTR_SIGNAL(Z) = SIGNAL(Z-ROWS/2)
    END DO

  END SUBROUTINE 

  SUBROUTINE flip_3d(CSCAN_IMAGE, LINE, TRA, ROWS, FLIP_CSCAN) 
  IMPLICIT NONE

    INTEGER :: X,F
    INTEGER, INTENT(IN) :: LINE,TRA,ROWS
    REAL*8, DIMENSION(LINE,TRA,ROWS), INTENT(IN) :: CSCAN_IMAGE
    REAL*8, DIMENSION(LINE,TRA,ROWS) :: FLIP_CSCAN

    DO X = 1, LINE
    DO F = 1, TRA
          FLIP_CSCAN(X,(TRA-F+1),:) = CSCAN_IMAGE(X,F,:)
    END DO
    END DO

  END SUBROUTINE 

  SUBROUTINE flip_2d(BSCAN_IMAGE, DIS, ROWS, FLIP_BSCAN) 
  IMPLICIT NONE

    INTEGER :: F
    INTEGER, INTENT(IN) :: DIS,ROWS
    REAL*8, DIMENSION(DIS,ROWS), INTENT(IN) :: BSCAN_IMAGE
    REAL*8, DIMENSION(DIS,ROWS) :: FLIP_BSCAN

    DO F = 1, DIS
          FLIP_BSCAN((DIS-F+1),:) = BSCAN_IMAGE(F,:)
    END DO

  END SUBROUTINE 

  SUBROUTINE flip_1D_INTEGER(SIGNAL, DIS, FLIP_SIGNAL) 
  IMPLICIT NONE

    INTEGER :: F
    INTEGER, INTENT(IN) :: DIS
    INTEGER, DIMENSION(DIS), INTENT(IN) :: SIGNAL
    INTEGER, DIMENSION(DIS) :: FLIP_SIGNAL

    DO F = 1, DIS
          FLIP_SIGNAL(DIS-F+1) = SIGNAL(F)
    END DO

  END SUBROUTINE 

  SUBROUTINE flip_1D_REAL(SIGNAL, DIS, FLIP_SIGNAL) 
  IMPLICIT NONE

    INTEGER :: F
    INTEGER, INTENT(IN) :: DIS
    REAL*8, DIMENSION(DIS), INTENT(IN) :: SIGNAL
    REAL*8, DIMENSION(DIS) :: FLIP_SIGNAL

    DO F = 1, DIS
          FLIP_SIGNAL(DIS-F+1) = SIGNAL(F)
    END DO

  END SUBROUTINE 

  SUBROUTINE background(signal,ROWS,e_1,e_2,bgr)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ROWS

  REAL*8, DIMENSION(1:ROWS), INTENT(IN) :: signal
  !Start from the 1 index. 
  !This is becasue I did not fix the number in the header.
  INTEGER, INTENT(IN) :: e_1
  INTEGER, INTENT(IN) :: e_2

  REAL*8 :: bgr_sum
  REAL*8 :: bgr

  INTEGER :: i
  INTEGER :: t

  bgr_sum = 0.0

  DO i = e_1,e_2,1

    bgr_sum = bgr_sum + signal(i)
    !print *, signal(i)
  END DO

  t = e_2 - e_1 + 1
  !    print *, "t=",t

  bgr =  bgr_sum / t
  !    print *, "bgr=",bgr

  END SUBROUTINE background


SUBROUTINE zero_padding(signal,ROWS,ROWS2,signal_ZEROPAD)
  INTEGER, INTENT(IN) :: ROWS 
  INTEGER, INTENT(IN) :: ROWS2
  REAL*8, DIMENSION(ROWS), INTENT(IN) :: signal
  REAL*8, DIMENSION(ROWS2) :: signal_ZEROPAD
  REAL*8, DIMENSION(ROWS) :: ZERO_PAD
  INTEGER :: Z, T, BAND, k, f_c 

  !=====ZERO_PADDING=====
  !ZERO_PAD = 1.0

  !GAUSSIAN 
  !BAND = 100
  !T = 3
  !DO Z = 1, ROWS
  !   ZERO_PAD(Z) = EXP(- (((Z-1)-(BAND * T)) / 2.0 * T) **2)
  !END DO 


  !BUTTERWORTH
  f_c = 450 
  k = 50
  DO Z = 1, ROWS
     ZERO_PAD(Z)  = SQRT(1.0/ (1.0 + ( (REAL(Z)/REAL(f_c))**(2.0 * k) ) ) ) 
  END DO 


  !ZERO PADDING
  DO Z = 1, ROWS2
     signal_ZEROPAD(Z) = signal(Z) * ZERO_PAD(Z)
  END DO 

END SUBROUTINE zero_padding




!reference Basic Earth Imaging writted by Jon Claerbout
!subroutine: adjnull, matmul
subroutine adjnull(adj, add, x, nx, y, ny)
integer  ix, iy,    adj, add,    nx,    ny
real                          x(nx), y(ny)
if(add == 0 ) then
       if(adj == 0) then
              do iy = 1, ny
                       y(iy) = 0
              end do 
       else 
              do ix = 1, nx
                       x(ix) =0 
              end do 
       end if 
end if 
!return; end 
end subroutine adjnull

!matrix multiply and its adjoint
subroutine matmult( adj, add, bb,         x,nx, y, ny)
integer ix, iy,     adj, add,               nx,    ny
real                          bb(ny, nx), x(nx), y(ny)
call adjnull(       adj, add,             x,nx,  y,ny)

do ix= 1, nx 
do iy= 1, ny
        if(add == 0 ) then
                       y(iy) = y(iy) + bb(iy, ix) * x(ix)
                       x(ix) = x(ix) + bb(iy, ix) * y(iy)
        end if 
end do 
end do

!return; end
end subroutine matmult

END MODULE MD_BASIC
