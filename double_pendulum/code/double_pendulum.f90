! double_pendulum.f90 — RK4 integration of the double pendulum (Fortran 2008).
! Compile:  gfortran -O2 -o double_pendulum double_pendulum.f90
! Run:       ./double_pendulum [input_file]     (default: input.dat)
!
! Input file: three Fortran namelists (default: input.dat), read in order:
!   &pendulum1                     ! first arm
!       m1 = 0.01,  l1 = 1.0       ! mass [kg], rod length [m]
!       th1 = 0.002, w1 = 0.0      ! initial angle [rad], angular velocity [rad/s]
!   /
!   &pendulum2                     ! second arm
!       m2 = 1.0,   l2 = 1.0
!       th2 = 0.5,  w2 = 0.0
!   /
!   &simulation                    ! time integration & output options
!       dt = 0.001                 ! time step [s]
!       steps = 50000              ! number of RK4 steps
!       poincare = .false.         ! .true.: write only th1 zero-crossings
!       perturb = 0.0              ! non-zero: second trajectory with
!                                  !   th2 -> th2 + perturb -> results_pert.dat
!   /
! All groups must appear in this order in the file. Omitted keys keep
! their defaults. A malformed namelist or unknown key is a runtime error.
!
! Output: results.dat (and results_pert.dat if perturb /= 0)
!   # t  th1  th2  w1  w2  E

module pendulum
   implicit none
   private
   public :: params_type, state_type, rhs, rk4_step, total_energy, &
             run, P

   real(8), parameter :: G = 9.81_8

   type :: params_type
      real(8) :: m1 = 0.01_8, m2 = 1.0_8
      real(8) :: l1 = 1.0_8, l2 = 1.0_8
      real(8) :: th1 = 0.002_8, th2 = 0.5_8
      real(8) :: w1 = 0.0_8, w2 = 0.0_8
      real(8) :: dt = 0.001_8
      integer :: steps = 50000
      logical :: poincare = .false.
      real(8) :: perturb = 0.0_8
   end type

   type :: state_type
      real(8) :: th1, th2, w1, w2
   end type

   type(params_type) :: P    ! module parameters, set in the main program

contains

   function rhs(s) result(ds)
      type(state_type), intent(in) :: s
      type(state_type) :: ds
      real(8) :: delta, den1

      delta = s%th2 - s%th1
      den1  = P%l1 * ((P%m1 + P%m2) - P%m2 * cos(delta)**2)

      ds%th1 = s%w1
      ds%th2 = s%w2
      ds%w1 = (0.5_8 * P%m2 * P%l1 * s%w1**2 * sin(2*delta)      &
             + P%m2 * P%l2 * s%w2**2 * sin(delta)                &
             + P%m2 * G * sin(s%th2) * cos(delta)                &
             - (P%m1 + P%m2) * G * sin(s%th1)) / den1
      ds%w2 = (-P%l1 * ds%w1 * cos(delta)                        &
             - P%l1 * s%w1**2 * sin(delta) - G * sin(s%th2)) / P%l2
   end function rhs

   function rk4_step(s, dt) result(out)
      type(state_type), intent(in) :: s
      real(8), intent(in) :: dt
      type(state_type) :: out, k1, k2, k3, k4

      k1 = rhs(s)
      k2 = rhs(state_type(s%th1 + 0.5_8*dt*k1%th1, s%th2 + 0.5_8*dt*k1%th2, &
                          s%w1  + 0.5_8*dt*k1%w1,  s%w2  + 0.5_8*dt*k1%w2))
      k3 = rhs(state_type(s%th1 + 0.5_8*dt*k2%th1, s%th2 + 0.5_8*dt*k2%th2, &
                          s%w1  + 0.5_8*dt*k2%w1,  s%w2  + 0.5_8*dt*k2%w2))
      k4 = rhs(state_type(s%th1 + dt*k3%th1, s%th2 + dt*k3%th2, &
                          s%w1  + dt*k3%w1,  s%w2  + dt*k3%w2))

      out%th1 = s%th1 + dt*(k1%th1 + 2*k2%th1 + 2*k3%th1 + k4%th1)/6
      out%th2 = s%th2 + dt*(k1%th2 + 2*k2%th2 + 2*k3%th2 + k4%th2)/6
      out%w1  = s%w1  + dt*(k1%w1  + 2*k2%w1  + 2*k3%w1  + k4%w1 )/6
      out%w2  = s%w2  + dt*(k1%w2  + 2*k2%w2  + 2*k3%w2  + k4%w2 )/6
   end function rk4_step

   function total_energy(s) result(E)
      type(state_type), intent(in) :: s
      real(8) :: E, cd, T, V
      cd = cos(s%th2 - s%th1)
      T  = 0.5_8*P%m1*P%l1**2*s%w1**2                              &
         + 0.5_8*P%m2*(P%l1**2*s%w1**2 + P%l2**2*s%w2**2           &
                       + 2*P%l1*P%l2*s%w1*s%w2*cd)
      V  = -(P%m1 + P%m2)*G*P%l1*cos(s%th1) - P%m2*G*P%l2*cos(s%th2)
      E  = T + V
   end function total_energy

   ! Writes every state including the final one (steps+1 rows in total).
   subroutine run(s0, outfile)
      type(state_type), intent(in) :: s0
      character(*), intent(in) :: outfile
      type(state_type) :: s
      logical :: prev_negative
      integer :: i, u

      open(newunit=u, file=outfile, status='replace')
      write(u,'(a)') '# t th1 th2 w1 w2 E'
      s = s0
      prev_negative = (s%th1 < 0)
      do i = 0, P%steps
         if (.not. P%poincare .or. (prev_negative .and. s%th1 >= 0)) then
            write(u,'(6(1x,es23.15))') i*P%dt, s%th1, s%th2, s%w1, s%w2, &
                                       total_energy(s)
         end if
         prev_negative = (s%th1 < 0)
         if (i < P%steps) s = rk4_step(s, P%dt)
      end do
      close(u)
   end subroutine run

end module pendulum

program double_pendulum
   use pendulum
   implicit none
   character(256) :: input
   integer :: u
   type(state_type) :: s0, sp

   ! ---- read parameters from the three namelists, in file order ----
   real(8) :: m1 = 0.01_8, l1 = 1.0_8, th1 = 0.002_8, w1 = 0.0_8
   real(8) :: m2 = 1.0_8,  l2 = 1.0_8, th2 = 0.5_8,   w2 = 0.0_8
   real(8) :: dt = 0.001_8, perturb = 0.0_8
   integer :: steps = 50000
   logical :: poincare = .false.
   namelist /pendulum1/ m1, l1, th1, w1
   namelist /pendulum2/ m2, l2, th2, w2
   namelist /simulation/ dt, steps, poincare, perturb

   call get_command_argument(1, input)
   if (len_trim(input) == 0) input = 'input.dat'
   open(newunit=u, file=trim(input), status='old', action='read')
   read(u, nml=pendulum1)
   read(u, nml=pendulum2)
   read(u, nml=simulation)
   close(u)

   P = params_type(m1=m1, m2=m2, l1=l1, l2=l2, th1=th1, th2=th2, &
                   w1=w1, w2=w2, dt=dt, steps=steps,             &
                   poincare=poincare, perturb=perturb)

   s0 = state_type(P%th1, P%th2, P%w1, P%w2)
   call run(s0, 'results.dat')

   if (P%perturb /= 0) then   ! second trajectory for Lyapunov / butterfly plots
      sp = state_type(P%th1, P%th2 + P%perturb, P%w1, P%w2)
      call run(sp, 'results_pert.dat')
   end if

   print '(a,es16.8)', 'DONE. E0 = ', total_energy(s0)
end program double_pendulum
