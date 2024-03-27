!Sun May 28 04:16:00 PM KST 2023
MODULE MD_KIRCHMIGR
IMPLICIT NONE

CONTAINS

subroutine simple_migration(permit, ROWS, DIS, r, s_n, ws_x, ws_z, hp_w, d, modl)

 real*4, intent(in) :: permit
 integer, intent(in) :: ROWS, DIS
 integer, intent(in) :: s_n
 integer, intent(in) :: ws_x, ws_z, hp_w
 
 real*8, intent(in) :: r

 integer :: t_n, w_c_x, w_c_z
 integer :: h_c_x, h_c_z
 real*8 :: c, T_0
 
 integer :: h_c,  h_s, h_e, v_s, v_e
 integer :: hp_s, hp_e 
 integer :: i, j, k, ii
 real :: x, depth, dt, t
 
 real*8, dimension(ROWS, DIS), intent(in) :: d
 !real*8, dimension(ROWS, DIS) :: modl
 real*8, dimension(ROWS,-ws_x:DIS+ws_x) :: modl


 c = 3.0*(10.0**8)

 T_0 = 2* r * sqrt(permit) / c 
 !T_0 = 2* r  / c 

 dt = T_0 / 4096.0 


 !moving  window according to highperbolic center
 do h_c_z =  1, ROWS, ws_z
 do h_c_x =  1, DIS, ws_x
 !print *, "h_c_x=",h_c_x, "d_n=",d_n
 
 !moving hipherbolic center in window
    h_s = h_c_x - ws_x/2.0 !hipher bolic start
    h_e = h_c_x + ws_x/2.0 -1 !hipher bolic end
 
    v_s = h_c_z  
    v_e = h_c_z + ws_z - 1
     
   !print *, "h_s=",h_s, "h_e=", h_e 
     
 do j = int(v_s), int(v_e)
          depth =  (j-s_n) * dt * c
 
 do i = int(h_s), int(h_e)   
   !print *, "j=",j, "i=",i, "k=", k 
     
         hp_s = i - hp_w/2.0
         hp_e = i + hp_w/2.0 
   
    do ii = int(hp_s), int(hp_e)
   
        x = abs(i - ii) * 0.5
        t = sqrt(depth**2 + x**2)/(c) * sqrt(permit)  
        !t = sqrt(depth**2 + x**2)/c   
        
        t_n = (j-s_n)/abs(j-s_n) * (t / dt) + s_n

        modl(j, i) = modl(j, i) + d(t_n,ii)
        !print *, "x=",x, "j=",j, "i=",i, "ii=", ii, "t_n=",t_n
        !print *, "x=",x, "i=",i, "t_n=",t_n, "h_c_z=",h_c_z, "depth=", depth

    end do 
        modl(j, i) = modl(j, i) / real(hp_e-hp_s+1)
     !print *, "                              " 
     
 end do
    !print *, "                              " 
 end do 
 
 end do    
 end do 

end subroutine

END MODULE MD_KIRCHMIGR
