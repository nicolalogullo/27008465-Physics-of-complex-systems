
! kicked_rotor.f90 — quantum kicked rotor and its classical limit.
! Compile:  gfortran -O2 -o kicked_rotor kicked_rotor.f90
! Run:      ./kicked_rotor [input_file]     (default: input.dat)
!
! Hamiltonian:  H = p^2/2 + K cos(x) * sum_n delta(t - n)
! Classical:   standard map  p' = p + K sin(x),  x' = x + p'  (mod 2pi)
! Quantum:     Floquet operator
!                  U = exp(-i p^2/(2 hbar)) * exp(-i K cos(x)/hbar)
!              applied with a split-step method: the free part is
!              diagonal in the momentum basis |m> (p = m*hbar, integer m),
!              the kick is diagonal in the position basis x_j = 2pi j/N.
!              The two bases are connected by a unitary DFT (O(N^2), no
!              external library needed). The split-step is EXACT for the
!              delta-kicked Hamiltonian (no Trotter error).
!
! Input file: three namelists, in this order:
!   &rotor
!       K = 5.0            ! kick strength
!       hbar = 1.0         ! effective Planck constant
!       N = 128            ! momentum basis size (power of 2 recommended)
!   /
!   &initial
!       m0 = 0              ! quantum: initial state = |m0>
!       x0 = 0.1, p0 = 0.0  ! classical: single-trajectory initial condition
!       seed = 12345        ! LCG seed for the classical ensemble
!   /
!   &simulation
!       nkicks = 2000       ! number of kicks (Floquet iterations)
!       stride = 10         ! write one row every `stride` kicks
!       ensemble = 2000     ! number of classical particles (diffusion)
!       quantum = .true., classical = .true.
!   /
!
! Outputs (header lines start with '#'):
!   quantum.txt    : # n p2 norm      (quantum <p^2> per kick: the KEY file)
!   momentum.txt   : # m prob         (final |psi_m|^2 distribution)
!   classical.txt  : # n p2_avg       (ensemble <p^2>: classical diffusion)
!   map.txt        : # n x p          (single classical trajectory:
!                                     the standard-map phase portrait)

module rotor_mod
   implicit none
   private
   public :: params_type, allocate_all, run_quantum, run_classical, P

   type :: params_type
      real(8) :: K = 5.0_8, hbar = 1.0_8
      integer :: N = 128
      integer :: m0 = 0
      real(8) :: x0 = 0.1_8, p0 = 0.0_8
      integer :: seed = 12345
      integer :: nkicks = 2000, stride = 10, ensemble = 2000
      logical :: quantum = .true., classical = .true.
   end type

   type(params_type) :: P

   ! quantum workspace (module-level: allocated once in allocate_all)
   complex(8), allocatable :: psi(:), psix(:), F(:, :), Fi(:, :)
   real(8), allocatable ::  free(:), kick(:)

   ! classical ensemble workspace
   real(8), allocatable :: cx(:), cp(:)

   ! deterministic LCG (reproducible runs, no external RNG)
   integer(8) :: lcg_state = 12345_8

contains

   real(8) function urand()
      lcg_state = mod(1103515245_8 * lcg_state + 12345_8, 2147483648_8)
      urand = real(lcg_state, 8) / 2147483648.0_8
   end function urand

   subroutine allocate_all()
      integer :: j, m, mv
      real(8) :: sq, pi

      pi = acos(-1._8)
      sq = 1.0_8 / sqrt(real(P%N, 8))

      allocate (psi(P%N), psix(P%N), free(P%N), kick(P%N))
      allocate (F(P%N, P%N), Fi(P%N, P%N))
      allocate (cx(max(1, P%ensemble)), cp(max(1, P%ensemble)))

      ! unitary DFT:  F(j,m) = exp(-2 pi i (j-1)(m-1)/N) / sqrt(N)
      !               Fi = conjugate transpose of F  (F^-1 = F^dag)
      ! NOTE: only F^-1 = F^dag matters here, because the kick phases are
      ! real (cos x): any overall sign convention cancels.
      do j = 1, P%N
         do m = 1, P%N
            F(j, m) = exp(-2._8 * (0._8, 1._8) * pi * real((j - 1) * (m - 1), 8) / P%N) * sq
            Fi(m, j) = conjg(F(j, m))
         end do
      end do

      ! momentum values:  mv = m - 1 - N/2   (integer),  p = mv * hbar
      ! free evolution: exp(-i p^2/(2 hbar)) = exp(-i mv^2 hbar / 2)
      do m = 1, P%N
         mv = m - 1 - P%N / 2
         free(m) = -0.5_8 * P%hbar * mv * mv
         psi(m) = merge(1.0_8, 0.0_8, mv == P%m0)   ! initial state |m0>
      end do

      ! kick phases: x_j = 2 pi (j-1)/N,  exp(-i K cos(x_j)/hbar)
      do j = 1, P%N
         kick(j) = -P%K * cos(2._8 * pi * real(j - 1, 8) / P%N) / P%hbar
      end do

      lcg_state = int(P%seed, 8)
      ! classical ensemble: x uniform in [0, 2pi), all p = p0
      do m = 1, max(1, P%ensemble)
         cx(m) = 2._8 * pi * urand()
         cp(m) = P%p0
      end do
   end subroutine allocate_all

   ! one EXACT Floquet step:  psi -> e^{-i p^2/2hbar} F^-1 e^{-iK cos x /hbar} F psi
   subroutine qstep()
      integer :: j, m
      complex(8) :: acc
      do j = 1, P%N                       ! to x basis
         acc = (0._8, 0._8)
         do m = 1, P%N
            acc = acc + F(j, m) * psi(m)
         end do
         psix(j) = cmplx(cos(kick(j)), sin(kick(j)), 8) * acc   ! kick phase
      end do
      do m = 1, P%N                       ! back to momentum basis
         acc = (0._8, 0._8)
         do j = 1, P%N
            acc = acc + Fi(m, j) * psix(j)
         end do
         psi(m) = cmplx(cos(free(m)), sin(free(m)), 8) * acc    ! free phase
      end do
   end subroutine qstep

   subroutine run_quantum()
      integer :: n, m, mv, u
      real(8) :: p2, nrm

      open (newunit=u, file='quantum.txt', status='replace')
      write (u, '(a)') '# n p2 norm'
      do n = 0, P%nkicks
         if (mod(n, P%stride) == 0) then
            p2 = 0._8; nrm = 0._8
            do m = 1, P%N
               mv = m - 1 - P%N / 2
               nrm = nrm + abs(psi(m))**2
               p2 = p2 + abs(psi(m))**2 * (mv * P%hbar)**2
            end do
            write (u, '(i0, 2(1x, es23.15))') n, p2, nrm
         end if
         if (n < P%nkicks) call qstep()
      end do
      close (u)

      open (newunit=u, file='momentum.txt', status='replace')
      write (u, '(a)') '# m prob'
      do m = 1, P%N
         mv = m - 1 - P%N / 2
         write (u, '(i0, 1x, es23.15)') mv, abs(psi(m))**2
      end do
      close (u)
      print '(a)', 'Data saved: quantum.txt, momentum.txt'
   end subroutine run_quantum

   ! ---- classical side ----

   ! standard map step (note: x uses the UPDATED p)
   subroutine cmap_step(x, p)
      real(8), intent(inout) :: x, p
      real(8) :: pi
      pi = acos(-1._8)
      p = p + P%K * sin(x)
      x = mod(x + p + 2._8 * pi * 10._8, 2._8 * pi)   ! wrap into [0, 2pi)
   end subroutine cmap_step

   subroutine run_classical()
      integer :: n, i, u, mv
      real(8) :: p2avg, xs, ps, pi
      pi = acos(-1._8)

      ! (a) ensemble diffusion -> classical.txt
      open (newunit=u, file='classical.txt', status='replace')
      write (u, '(a)') '# n p2_avg'
      do n = 0, P%nkicks
         if (mod(n, P%stride) == 0) then
            p2avg = 0._8
            do i = 1, P%ensemble
               p2avg = p2avg + cp(i)**2
            end do
            p2avg = p2avg / P%ensemble
            write (u, '(i0, 1x, es23.15)') n, p2avg
         end if
         do i = 1, P%ensemble
            call cmap_step(cx(i), cp(i))
         end do
      end do
      close (u)

      ! (b) single trajectory -> map.txt (standard-map phase portrait)
      xs = P%x0; ps = P%p0
      open (newunit=u, file='map.txt', status='replace')
      write (u, '(a)') '# n x p'
      do n = 0, P%nkicks
         write (u, '(i0, 2(1x, es23.15))') n, merge(xs - 2._8 * pi, xs, xs > pi), ps
         call cmap_step(xs, ps)
      end do
      close (u)
      print '(a)', 'Data saved: classical.txt, map.txt'
   end subroutine run_classical

end module rotor_mod

!---------------------------------------------------------------------
program kicked_rotor
   use rotor_mod
   implicit none
   character(256) :: fname
   integer :: u

   ! ---- namelist variables (defaults = the reference run of the theory) ----
   real(8) :: K = 5.0_8, hbar = 1.0_8
   integer :: N = 128, m0 = 0
   real(8) :: x0 = 0.1_8, p0 = 0.0_8
   integer :: seed = 12345, nkicks = 2000, stride = 10, ensemble = 2000
   logical :: quantum = .true., classical = .true.
   namelist /rotor/      K, hbar, N
   namelist /initial/    m0, x0, p0, seed
   namelist /simulation/ nkicks, stride, ensemble, quantum, classical

   call get_command_argument(1, fname)
   if (len_trim(fname) == 0) fname = 'input.dat'
   open (newunit=u, file=trim(fname), status='old', action='read')
   read (u, nml=rotor)
   read (u, nml=initial)
   read (u, nml=simulation)
   close (u)

   P = params_type(K=K, hbar=hbar, N=N, m0=m0, x0=x0, p0=p0, seed=seed, &
                   nkicks=nkicks, stride=stride, ensemble=ensemble,     &
                   quantum=quantum, classical=classical)

   if (mod(N, 2) /= 0) then
      print '(a)', 'ERROR: N must be even (momentum grid symmetric in m=0).'
      stop 1
   end if

   call allocate_all()
   if (P%quantum) call run_quantum()
   if (P%classical) call run_classical()

   print '(a)', 'DONE.'
end program kicked_rotor
