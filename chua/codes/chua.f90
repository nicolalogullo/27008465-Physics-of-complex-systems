
! chua.f90 — RK4 integration of the (dimensionless) Chua oscillator.
! Compile:  gfortran -O2 -o chua chua.f90
! Run:      ./chua [input_file]     (default: input.dat)
!
! Input file: three namelists, in this order:
!   &circuit
!       alpha = 9.0, beta = 14.286, gamma = 0.05  ! dimensionless parameters
!       a = -1.143, b = -0.714                    ! Chua-diode slopes (inner, outer)
!   /
!   &initial
!       x0 = 0.1, y0 = 0.0, z0 = 0.0   ! initial state
!   /
!   &simulation
!       dt = 0.001        ! time step [dimensionless]
!       steps = 200000    ! number of RK4 steps
!       stride = 100      ! write one row every `stride` steps
!       poincare = .true. ! also write Poincare section (y crossing 0 upward)
!       perturb = 0.0     ! non-zero: second trajectory with x0 -> x0+perturb
!       analysis = .false.! .true.: print equilibria + eigenvalues, no integration
!   /
!
! Outputs (header lines start with '#'):
!   trajectory.txt  : # tau x y z   (one row every `stride` steps)
!   poincare.txt    : # x z         (states where y crosses 0 upward)
!   results_pert.txt: like trajectory.txt, only if perturb /= 0
!
! With analysis = .true. the program does NOT integrate: it computes the
! three equilibria, checks the existence conditions, prints the
! Jacobian eigenvalues at each, and stops — the numerical counterpart
! of the linear-stability analysis in theory.md (Exercise 1.3).

module chua_mod
   implicit none
   private
   public :: params_type, state_type, rhs, rk4_step, write_trajectory, &
             write_poincare, P

   type :: params_type
      real(8) :: alpha = 9.0_8, beta = 14.286_8, gamma = 0.05_8
      real(8) :: a = -1.143_8, b = -0.714_8
      real(8) :: dt = 0.001_8, perturb = 0.0_8
      integer :: steps = 200000, stride = 100
   end type

   type :: state_type
      real(8) :: x, y, z
   end type

   type(params_type) :: P    ! set once in the main program

contains

   ! Nonlinear-resistor characteristic (dimensionless):
   !   f(x) = b x + (a-b)/2 (|x+1| - |x-1|)
   ! inner slope a (|x| < 1), outer slope b (|x| > 1)
   function f(x) result(v)
      real(8), intent(in) :: x
      real(8) :: v
      v = P%b * x + 0.5_8 * (P%a - P%b) * (abs(x + 1) - abs(x - 1))
   end function f

   function rhs(s) result(ds)
      type(state_type), intent(in) :: s
      type(state_type) :: ds
      ds%x = P%alpha * (s%y - s%x - f(s%x))
      ds%y = s%x - s%y + s%z
      ds%z = -P%beta * s%y - P%gamma * s%z
   end function rhs

   function rk4_step(s, dt) result(out)
      type(state_type), intent(in) :: s
      real(8), intent(in) :: dt
      type(state_type) :: out, k1, k2, k3, k4

      k1 = rhs(s)
      k2 = rhs(state_type(s%x + 0.5_8*dt*k1%x, s%y + 0.5_8*dt*k1%y, &
                          s%z + 0.5_8*dt*k1%z))
      k3 = rhs(state_type(s%x + 0.5_8*dt*k2%x, s%y + 0.5_8*dt*k2%y, &
                          s%z + 0.5_8*dt*k2%z))
      k4 = rhs(state_type(s%x + dt*k3%x, s%y + dt*k3%y, s%z + dt*k3%z))

      out%x = s%x + dt*(k1%x + 2*k2%x + 2*k3%x + k4%x)/6
      out%y = s%y + dt*(k1%y + 2*k2%y + 2*k3%y + k4%y)/6
      out%z = s%z + dt*(k1%z + 2*k2%z + 2*k3%z + k4%z)/6
   end function rk4_step

   ! Full trajectory -> file with "# tau x y z" (every `stride` steps)
   subroutine write_trajectory(s0, outfile)
      type(state_type), intent(in) :: s0
      character(*), intent(in) :: outfile
      type(state_type) :: s
      integer :: i, u

      open(newunit=u, file=outfile, status='replace')
      write(u,'(a)') '# tau x y z'
      s = s0
      do i = 0, P%steps
         if (mod(i, P%stride) == 0) &
            write(u,'(4(1x,es23.15))') i*P%dt, s%x, s%y, s%z
         s = rk4_step(s, P%dt)
      end do
      close(u)
   end subroutine write_trajectory

   ! Poincare section: states where y crosses 0 upward -> "# x z"
   subroutine write_poincare(s0, outfile)
      type(state_type), intent(in) :: s0
      character(*), intent(in) :: outfile
      type(state_type) :: s
      logical :: prev_neg
      integer :: i, u

      open(newunit=u, file=outfile, status='replace')
      write(u,'(a)') '# x z'
      s = s0
      prev_neg = (s%y < 0)
      do i = 1, P%steps
         s = rk4_step(s, P%dt)
         if (prev_neg .and. s%y >= 0) write(u,'(2(1x,es23.15))') s%x, s%z
         prev_neg = (s%y < 0)
      end do
      close(u)
   end subroutine write_poincare

end module chua_mod

!---------------------------------------------------------------------
! Equilibrium / linear-stability analysis (analysis = .true.)
!---------------------------------------------------------------------
module chua_analysis
   use chua_mod
   implicit none
   private
   public :: analyze_equilibria

contains

   ! Jacobian at an equilibrium whose x* lies in the region with slope m:
   !   J = [ -alpha(1+m), alpha, 0 ; 1, -1, 1 ; 0, -beta, -gamma ]
   ! Characteristic polynomial: det(lambda I - J)
   !      = lambda^3 + c2 lambda^2 + c1 lambda + c0,
   ! with the closed forms below (derive them yourself: Exercise 1.3).
   subroutine equilibria_print(xe, ye, ze, m)
      real(8), intent(in) :: xe, ye, ze, m
      real(8) :: c2, c1, c0, lam(3)

      c2 = P%gamma + 1 + P%alpha*(1 + m)
      c1 = P%gamma + P%beta + P%alpha*(1 + m)*(P%gamma + 1) - P%alpha
      c0 = P%alpha*((1 + m)*(P%gamma + P%beta) - P%gamma)
      call solve_cubic(c2, c1, c0, lam)

      print '(a,3(1x,f10.6))', '  equilibrium (x,y,z) =', xe, ye, ze
      print '(a,1x,f6.3)',     '  f''(x*) =', m
      print '(a,3(1x,f12.6))', '  Jacobian eigenvalues =', lam
      print '(a,1x,f12.6)',    '  sum of eigenvalues   =', lam(1)+lam(2)+lam(3)
      print '(a,1x,f12.6)',    '  (must equal -trace J) =', -c2
      print *
   end subroutine equilibria_print

   ! Real roots of lambda^3 + c2 lambda^2 + c1 lambda + c0 (Cardano).
   ! If there is a complex-conjugate pair, only the real parts are shown
   ! (with a note); lam(3) is always sorted ascending.
   subroutine solve_cubic(c2, c1, c0, lam)
      real(8), intent(in) :: c2, c1, c0
      real(8), intent(out) :: lam(3)
      real(8) :: p, q, disc, aa, bb, r, phi, arg
      real(8), parameter :: PI = acos(-1._8)

      ! depressed cubic u^3 + p u + q = 0 with lambda = u - c2/3
      p = c1 - c2*c2/3
      q = c0 - c1*c2/3 + 2*c2**3/27
      disc = (q/2)**2 + (p/3)**3

      if (disc > 0) then                        ! one real root + complex pair
         aa = cube_root(-q/2 + sqrt(disc))
         bb = cube_root(-q/2 - sqrt(disc))
         lam(1) = aa + bb - c2/3
         lam(2) = -0.5_8*(aa + bb) - c2/3        ! real part of the pair
         lam(3) = lam(2)
         print '(a)', '  [complex pair: real parts shown]'
      else                                      ! three real roots
         r = 2*sqrt(-p/3)
         arg = (3*q)/(2*p)*sqrt(-3/p)
         phi = acos(max(-1._8, min(1._8, arg)))
         lam(1) = r*cos(phi/3)              - c2/3
         lam(2) = r*cos((phi + 2*PI)/3)      - c2/3
         lam(3) = r*cos((phi + 4*PI)/3)      - c2/3
         call sort3(lam)
      end if
   end subroutine solve_cubic

   real(8) function cube_root(v)
      real(8), intent(in) :: v
      cube_root = sign(abs(v)**(1._8/3), v)
   end function cube_root

   subroutine sort3(a)
      real(8), intent(inout) :: a(3)
      real(8) :: t
      if (a(1) > a(2)) then; t = a(1); a(1) = a(2); a(2) = t; end if
      if (a(2) > a(3)) then; t = a(2); a(2) = a(3); a(3) = t; end if
      if (a(1) > a(2)) then; t = a(1); a(1) = a(2); a(2) = t; end if
   end subroutine sort3

   subroutine analyze_equilibria()
      real(8) :: xpm, ye, ze

      ! Outer equilibria P+/P- (verify against theory.md, Exercise 1.3):
      !   x* = (b-a)(gamma+beta) / (beta + b*(gamma+beta))
      !   y* = gamma*x*/(beta+gamma),   z* = -beta*x*/(beta+gamma)
      ! NOTE the proportionality: (y*,z*) scale with x*, they are NOT
      ! the constants gamma/(beta+gamma) and -beta/(beta+gamma)!
      xpm = (P%b - P%a)*(P%gamma + P%beta)/(P%beta + P%b*(P%gamma + P%beta))
      ye = P%gamma*xpm/(P%beta + P%gamma)
      ze = -P%beta*xpm/(P%beta + P%gamma)

      print '(a)', '=== Chua oscillator: equilibria and stability ==='
      print '(a,f8.4,a,f8.4,a,f8.4)', &
            ' parameters: alpha =', P%alpha, '  beta =', P%beta, &
            '  gamma =', P%gamma
      print *
      print '(a)', 'P0 (central region, |x| < 1):'
      call equilibria_print(0._8, 0._8, 0._8, P%a)

      print '(a)', 'P+ (outer region, x > 1), exists only if x+ > 1:'
      print '(a,1x,f10.6)', '  x+ =', xpm
      if (xpm > 1) then
         call equilibria_print(xpm, ye, ze, P%b)
      else
         print '(a)', '  does NOT exist for these parameters (x+ <= 1)'
         print *
      end if

      print '(a)', 'P- (outer region, x < -1): symmetric image of P+'
      if (xpm > 1) then
         call equilibria_print(-xpm, -ye, -ze, P%b)
      else
         print '(a)', '  does NOT exist for these parameters'
      end if
   end subroutine analyze_equilibria

end module chua_analysis

!---------------------------------------------------------------------
program chua
   use chua_mod
   use chua_analysis
   implicit none
   character(256) :: fname
   type(state_type) :: s0, sp

   ! ---- namelist variables (defaults = canonical double scroll) ----
   real(8) :: alpha = 9.0_8, beta = 14.286_8, gamma = 0.05_8
   real(8) :: a = -1.143_8, b = -0.714_8
   real(8) :: x0 = 0.1_8, y0 = 0.0_8, z0 = 0.0_8
   real(8) :: dt = 0.001_8, perturb = 0.0_8
   integer :: steps = 200000, stride = 100
   logical :: poincare = .true., analysis = .false.
   namelist /circuit/    alpha, beta, gamma, a, b
   namelist /initial/    x0, y0, z0
   namelist /simulation/ dt, steps, stride, poincare, perturb, analysis

   call get_command_argument(1, fname)
   if (len_trim(fname) == 0) fname = 'input.dat'
   open(newunit=u, file=trim(fname), status='old', action='read')
   read(u, nml=circuit)
   read(u, nml=initial)
   read(u, nml=simulation)
   close(u)

   P = params_type(alpha=alpha, beta=beta, gamma=gamma, a=a, b=b,   &
                  dt=dt, perturb=perturb, steps=steps, stride=stride)

   if (analysis) then
      call analyze_equilibria()
      stop
   end if

   s0 = state_type(x0, y0, z0)
   call write_trajectory(s0, 'trajectory.txt')
   print '(a)', 'Data saved: trajectory.txt'

   if (poincare) then
      call write_poincare(s0, 'poincare.txt')
      print '(a)', 'Data saved: poincare.txt'
   end if
   if (perturb /= 0) then
      sp = state_type(x0 + perturb, y0, z0)
      call write_trajectory(sp, 'results_pert.txt')
      print '(a)', 'Data saved: results_pert.txt'
   end if

end program chua
