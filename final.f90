
!/  FORTRAN PROGRAM TO COMPUTE A COUETTE-POISEUILLE FLOW WITH A MIXING LENGTH MODEL
!/  Y and U are the arrays used to store the co-ordinates of mesh and corresponding velocity
!   Viscosity= vis, density = rho, Pressure Gradient = dp_dx, Vw=Wall Velocity, T_vis=Turbulent viscosity
!/  u_tau1 and u_tau2 are the friction velocities for the lower and upper half of the channel respectively
!/  lm, lo, A are the mixing length model parameters
!/  Relax is the relaxation co-efficient
!/  U_scale is the velocity to scale the flow velocity. It is either equal to Vw or U_max
!/  U is the new flow velocity and U_old is the velocity calculated just before one step
!/  delta is the difference between U and U_old
!/  y_m is the array to store the co-ordinates at the middle points (N and S) and u_m is the corresponding velocity




program main
implicit none

  !!!!Declaring the variables

    Integer                             :: i, N, N1,N05, k, iterations,iter
    real*8                              :: Vw, H, rho, vis, dp_dx, DELI, DELJ, relax, u_tau1, u_tau2, delta,u_mean,u_max, u_scale
    Real*8, allocatable, dimension(:)   :: a, b, c, d, y, u, u_exact,  YPAD, YHAD, m, u_old, T_vis, lo, lm, u_m, y_m





!!!Reading Input data
open(unit=25, file='input.dat', action = 'read')
read(25,*) rho, vis, relax, N1
 close (25)


!!! Reading dp_dx, Vw ahd H from the keyboard
write(*,*) "Enter the Pressure Gradient"
read(*,*) dp_dx
write(*,*) "Enter the Wall Velocity"
read(*,*) Vw
write(*,*) "Enter the Height of the Channel"
read(*,*) H
write(*,*) "Enter the Number of the mesh points"
read(*,*) N1


!!Number of mesh points (same as N1 in the Mesh Generating Program)
N=N1
N05=(N1+1)/2


!!!Giving the dimensions to the allocatable arrays
allocate (a(N),b(N),c(N),d(N), y(N),m(N),YPAD(N),YHAD(N), u(N), u_exact(N), u_old(N), T_vis(N-1), lo(N-1), lm(N-1), u_m(N-1))


!!!!Opening the Mesh and reading the mesh parameters
open(unit=100, file='mesh.dat', action = 'read')
    do i=1, N
    read(100,*) m(i), YPAD(i), YHAD(i)
    enddo
    close(100)



!!!!! Adapting the Non-uniform mesh
 do i = 1, N
 y(i) = YPAD(i)*H
 end do




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        Turbulent Solution         !!!!!!!!!!!!!!!!!!!!!!!!!1!!!!!!!!!!!

! Initialising the initial velocity profile
do i= 1, N
u(i) = y(i)
end do


! Iteration loop for Turbulent viscosity
delta = 1.
iterations =0
do iter=1,10
u_old = u
iterations=iterations+1

call Turbulent_Viscosity(N,N05,rho,vis, H, y, u_old, T_vis,u_tau1,u_tau2)

do i= 2, N-1
    DELI=y(i)-y(i-1)
    DELJ=y(i+1)-y(i)
    a(i)=(vis+T_vis(i-1))/(DELI)
    c(i)=(vis+T_vis(i))/(DELJ)
    b(i) = -(a(i)+c(i))
    d(i)=0.5*dp_dx*(DELJ+DELI)
end do

 !Boundary conditions
    b(1) = 1.0
    c(1) = 1.0
    d(1) = 0.0
    a(N) = 1.0
    b(N) = 1.0
    d(N) = Vw*2.0

call thomas(a,b,c,d,u,N)
call Relaxation (N, u_old, u, relax)

!!!delta between u_old and u
do i =1, N
delta = delta+(abs(u_old(i)-u(i)))**2.
delta = sqrt(delta)/N
write(*,*) u_old(i), u(i), delta
end do

!u_old=u
 end do

 !!!!Mean Velocity
 !Get velocity at the middle points
 do i=1, N-1
 u_m(i) = (u(i)+u(i+1))/2.
 end do
 u_mean = 0
 do i=1, N-1
 u_mean = u_mean+u_m(i)*(y(i+1)-y(i))/H
 end do
 write(*,*) u_mean

!!!Getting the Maximum velocity
u_max= maxval(u)


 !SCALING HEIGHT AND VEOCITY
 !!! Important: For 3 cases of Anne Gilliot, scaling of velocity is done by u_mean instead of u_max and Vw)
  If(Vw.EQ.0.0) then
    U_scale=maxval(u)
  else
    U_scale=abs(Vw)
  END IF


   !!Writing the data
    write(*,*) iterations, u_tau1, u_tau2, u_max, u_mean

!!!Opening dat file for writing the scaled velocity profile
 open(unit=20, file='Velprofile.dat', action = 'write')
    do k=1, N
    write(*,*) y(k)/H, u(k)/U_scale
    end do
    close (20)



!!! Velocity profile without scaling
open(unit=15, file='iter.dat', action = 'write')
    do k=1, N
    write(15,*) y(k), u(k)
    end do
    close (15)

!!!!! Writing the turbulent viscosity
open(unit=17, file='turbvis.dat', action = 'write')
    do k=1, N
    write(17,*)  T_vis(i)
    end do
    close (17)

 !! Writing the velocity profile scaled by the fixed wall variables
  open(unit=30, file='fixedwall.dat', action = 'write')
    do i =1,N
    write(30,*) y(i)*u_tau1/vis, u(i)/u_tau1
    end do
    close(30)


    !Writing the U_old and the new velocity U and their difference,delta
    open(unit=40, file='error.dat', action = 'write')
    do i=1, N
    write(40,*)  u_old(i), u(i), delta
    end do
    close (40)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
subroutine thomas(a,b,c,d,u,N)

    implicit none
    integer :: N,i
    real*8, dimension(N)  :: a, b, c, d, u
real*8::delta
! Forward Elimination
    do i=2,N
           delta=a(i)/b(i-1)
        !a(i)=a(i)-delta*b(i-1)
        b(i)=b(i)-(delta)*c(i-1)
        d(i)=d(i)-(delta)*d(i-1)
    enddo

!Backward substitution
    u(N)=d(N)/b(N)
    do i= N-1,1,-1
    u(i)= (d(i)-c(i)*u(i+1))/b(i)
    enddo

end subroutine thomas

subroutine Analytical_Solution(N, Vw, dp_dx, H, vis, rho, y, u_exact)
implicit none
integer                         :: i, N
real*8                          :: H, Vw, dp_dx, rho, vis
real*8, dimension (N)           :: y, u_exact


do i= 1, N
u_exact(i) = (dp_dx/(2*vis))*y(i)*(y(i)-H)+Vw*y(i)/H
enddo

end subroutine Analytical_Solution


subroutine Turbulent_Viscosity(N, N05, rho, vis, H, y, u, T_vis, u_tau1, u_tau2)
implicit none
integer                     :: i, N, N05
real*8                      :: u_tau1, u_tau2, tau_w, A, vis, rho, H
real*8, dimension(N)        :: u, y
real*8, dimension(N-1)      ::  T_vis, du_dy, lo, lm, Ym



! Turbulent viscosity for the lower half of the channel
tau_w = vis*abs((u(2)-u(1))/(y(2)-y(1)))
u_tau1 = sqrt(tau_w/rho)
A = 26.*vis/u_tau1
write(*,*) 'utau1',u_tau1

do i = 1, N05
!GETTING THE MIDDLE POINTS
Ym(i) = (y(i) + y(i+1))/2.

du_dy(i) = abs((u(i+1)-u(i))/(y(i+1)-y(i)))
lo(i) = 0.5*H*(0.14-(0.08*(1-2*Ym(i)/H)**2)-(0.06*(1-2*Ym(i)/H)**4))
lm(i) = lo(i)*(1-exp(-Ym(i)/A))
T_vis(i) = lm(i)**2.*du_dy(i)
end do


!Turbulent viscosity for the upper half of the channel
tau_w = vis*abs((u(N)-u(N-1))/(y(N)-y(N-1)))
u_tau2 = sqrt(tau_w/rho)
A = 26.*vis/u_tau2
write(*,*) 'utau2',u_tau2

do i = N, N05+1, -1
Ym(i) = (y(i) + y(i+1))/2.

du_dy(i) = abs((u(i+1)-u(i))/(y(i+1)-y(i)))
lo(i) = 0.5*H*(0.14-(0.08*(1.-2.*(H-Ym(i))/H)**2)-(0.06*(1-2*Ym(i)/H)**4))
lm(i) = lo(i)*(1-exp(-(H-Ym(i))/A))
T_vis(i) = lm(i)**2.*du_dy(i)
end do
end subroutine Turbulent_Viscosity


subroutine Relaxation (N, u_old, u , relax)
implicit none
integer :: i, N
real*8 :: relax
real*8, dimension(N) :: u_old, u

do i = 1, N
u(i) = relax*u(i)+ (1.-relax)*u_old(i)
end do
end subroutine Relaxation


end program main

