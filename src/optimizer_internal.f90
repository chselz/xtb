! This file is part of xtb.
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

!> Quasi-Newton minimisation in redundant internal coordinates.
!>
!> Fortran re-implementation of the minimum search of Sella 2.6.0
!> (E. D. Hermes, K. Sargsyan, H. N. Najm, J. Zador, J. Chem. Theory Comput.
!> 18, 6974 (2022); https://github.com/zadorlab/sella, LGPL-3.0), restricted
!> to what is needed here:
!>
!> * redundant primitives (bonds, angles, linear bends, dihedrals) of the
!>   covalent graph (xtb_internals_redundant),
!> * diagonal Fischer-Almloef guess Hessian, projected onto the range of B,
!> * TS-BFGS update of the full redundant Hessian, never reset,
!> * quasi-Newton step (|eigenvalues|) in the non-redundant subspace with a
!>   level shift that limits the largest change of a single internal
!>   coordinate ("max internal step"),
!> * Sella's trust-radius rules for minimisations,
!> * iterative back-transformation to Cartesian coordinates.
!>
!> Additions beyond Sella:
!> * covalently disconnected fragments get translation and linearised
!>   rotation coordinates (TRIC-like, re-referenced every cycle) instead of
!>   Sella's weld bonds,
!> * the largest Cartesian atom displacement per step is capped (0.5 bohr).
!>
!> Not supported (falls back to ANCopt): fixed atoms, periodic systems.
module xtb_optimizer_internal
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_mctc_param, only : covalent_radius_d3
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_restart, only : TRestart
   use xtb_type_calculator, only : TCalculator
   use xtb_type_data, only : scc_results
   use xtb_io_writer, only : writeMolecule
   use xtb_mctc_fileTypes, only : fileType
   use xtb_type_neighbourlist, only : TNeighbourList, init
   use xtb_internals_graph, only : graph_type, init
   use xtb_internals_type, only : coord_bond, coord_angle, coord_dihedral, coord_linbend
   use xtb_internals_redundant, only : redundant_type, init
   use xtb_bmatrix, only : get_bmatrix
   implicit none
   private

   public :: intopt, intopt_preferred

   real(wp), parameter :: pi = 3.14159265358979323846_wp
   !> Sella defaults for minima (rs = 'mis', method = 'qn')
   real(wp), parameter :: delta0 = 0.1_wp       ! Angstrom / rad
   real(wp), parameter :: delta_min = 1.0e-4_wp
   real(wp), parameter :: sigma_inc = 1.15_wp, sigma_dec = 0.90_wp
   real(wp), parameter :: rho_inc = 1.035_wp, rho_dec = 100.0_wp
   !> Largest Cartesian displacement of a single atom per step (bohr)
   real(wp), parameter :: max_atom_displ = 0.5_wp
   !> Initial curvature of the fragment translation/rotation coordinates (Eh/bohr^2)
   real(wp), parameter :: frag_curvature = 0.05_wp

contains

!> Heuristic engine choice for the automatic mode (--opt): use intopt for
!> molecules made of main-group elements only (H, B-F, Si-Cl, Ge-Br, Se,
!> Sn-I) with at most six covalent neighbours per atom, up to maxat atoms
!> and without fixed atoms or periodicity; ANCopt otherwise.
!> Benchmark basis (GFN-FF, 108 systems): intopt needed fewer cycles for
!> 73 of the 79 metal-free systems (incl. NCI dimers and solvated clusters)
!> but was clearly worse for nearly all systems containing metal atoms
!> (carbonyls, sandwich complexes, metal and ionic clusters, M-crown).
function intopt_preferred(mol, maxat, reason) result(yes)
   use xtb_fixparam, only : fixset
   type(TMolecule), intent(in) :: mol
   integer, intent(in) :: maxat
   character(len=:), allocatable, intent(out) :: reason
   logical :: yes
   type(TNeighbourList) :: neigh
   type(graph_type) :: graph
   character(len=16) :: buf

   yes = .false.
   if (mol%npbc > 0) then
      reason = "periodic system"
   else if (fixset%n > 0) then
      reason = "fixed atoms"
   else if (mol%n > maxat) then
      write(buf, '(i0)') maxat
      reason = "more than " // trim(buf) // " atoms"
   else if (any(.not.organic(mol%at))) then
      reason = "metal or noble-gas atoms present"
   else
      call init(neigh, mol%n)
      call neigh%generate_covalent(mol%at, mol%xyz)
      call init(graph, neigh)
      if (maxval(graph%neighs(:mol%n)) > 6) then
         reason = "atoms with more than six covalent neighbours"
      else
         reason = "main-group covalent topology"
         yes = .true.
      end if
   end if
end function intopt_preferred


subroutine intopt(env, ilog, mol, chk, calc, egap, et, maxiter, maxcycle_in, etot, g, sigma, &
      & tight, pr, fail, iter_needed)
   use xtb_optimizer, only : ancopt, get_optthr
   use xtb_fixparam, only : fixset
   character(len=*), parameter :: source = "optimizer_intopt"
   type(TEnvironment), intent(inout) :: env
   integer, intent(in) :: ilog
   type(TMolecule), intent(inout) :: mol
   type(TRestart), intent(inout) :: chk
   class(TCalculator), intent(inout) :: calc
   real(wp), intent(inout) :: egap
   real(wp), intent(in) :: et
   integer, intent(in) :: maxiter, maxcycle_in, tight
   real(wp), intent(inout) :: etot
   real(wp), intent(inout) :: g(3, mol%n)
   real(wp), intent(inout) :: sigma(3, 3)
   logical, intent(in) :: pr
   logical, intent(out) :: fail
   integer, intent(out), optional :: iter_needed

   type(redundant_type) :: prim
   type(scc_results) :: res
   integer :: nat, n3, np, nq, nf, iter, maxopt, ncomp, i, prlevel
   integer, allocatable :: frag(:), cn(:)
   integer :: nmet
   real(wp) :: ethr, gthr, acc, energy, eold, echng, gnorm, delta, rho, dfpred, smax, dmax
   real(wp), allocatable :: q(:), qold(:), gq(:), gqold(:), h(:,:), h0(:), w(:)
   real(wp), allocatable :: bq(:,:), basis(:,:), binvt(:,:), s(:), xold(:,:), sact(:)
   real(wp), allocatable :: bl(:,:), blold(:,:), hf(:,:)
   logical :: have_step, conv, ok
   character(len=*), parameter :: scifmt = '(10x,":",3x,a,e22.7,1x,a,1x,":")'
   character(len=*), parameter :: dblfmt = '(10x,":",3x,a,f18.7,5x,a,1x,":")'
   character(len=*), parameter :: intfmt = '(10x,":",3x,a,i18,      10x,":")'
   integer(8) :: c0, c1, crate
   real(wp) :: tm(6)

   fail = .false.
   nat = mol%n
   n3 = 3*nat
   call get_optthr(nat, tight, ethr, gthr, maxopt, acc)
   if (maxcycle_in > 0) maxopt = maxcycle_in

   ! primitives of the covalent graph (+ fragment coordinates)
   call build_primitives(mol%at, mol%xyz, prim, ncomp, frag, cn)
   np = prim%ncoords
   if (pr) then
      nmet = count(.not.organic(mol%at))
      write(env%unit, '(1x,a,i6,i5,i7,i6,i5,i5,i5)') &
         & "intopt topology: nat ncomp nprim nmetal maxcn ncn>4 ncn0", &
         & nat, ncomp, np, nmet, maxval(cn), count(cn > 4), count(cn == 0)
   end if
   if (ncomp > 1) then
      call fragment_rows(mol%xyz, frag, ncomp, bl)
   else
      allocate(bl(0, n3))
   end if
   nq = np + size(bl, 1)
   if (fixset%n > 0 .or. mol%npbc > 0 .or. nq < 1) then
      if (pr) write(env%unit, '(/,1x,a)') &
         & "intopt: fixed atoms or PBC -> using ANCopt"
      call ancopt(env, ilog, mol, chk, calc, egap, et, maxiter, maxcycle_in, etot, g, sigma, &
         & tight, pr, fail, iter_needed)
      return
   end if

   allocate(q(nq), qold(nq), gq(nq), gqold(nq), h(nq, nq), h0(nq), w(nq), bq(nq, n3), &
      & s(nq), sact(nq), xold(3, nat))
   ! weights for the max-internal-step criterion: lengths in Angstrom, angles in rad
   w = 1.0_wp/aatoau
   do i = 1, np
      if (prim%kind(i) /= coord_bond) w(i) = 1.0_wp
   end do

   ! initial Hessian: P diag(h0) P on the range of B
   call bmat_all(prim, bl, mol%xyz, bq)
   call nonredundant_basis(bq, basis, binvt, nf)
   if (np > 0) call fa_guess(prim, mol%at, mol%xyz, h0(:np))
   h0(np+1:) = frag_curvature
   h = 0.0_wp
   do i = 1, nq
      h(i, i) = h0(i)
   end do
   ! h = U (U^T diag(h0) U) U^T
   h = mm(basis, mm(mm(basis, spread(h0, 2, nf)*basis, 'T', 'N'), basis, 'N', 'T'))

   if (pr) then
      write(env%unit, '(/,10x,51("."))')
      write(env%unit, '(10x,":",17x,a,18x,":")') "SETUP (intopt)"
      write(env%unit, '(10x,":",49("."),":")')
      write(env%unit, intfmt) "covalent fragments", ncomp
      write(env%unit, intfmt) "fragment T/R crds.", size(bl, 1)
      write(env%unit, intfmt) "primitive internal", np
      write(env%unit, intfmt) "degrees of freedom", nf
      write(env%unit, intfmt) "max. optcycles    ", maxopt
      write(env%unit, '(10x,":",49("."),":")')
      write(env%unit, scifmt) "energy convergence", ethr, "Eh  "
      write(env%unit, scifmt) "grad. convergence ", gthr, "Eh/a"
      write(env%unit, dblfmt) "max. atom displ.  ", max_atom_displ, "a0  "
      write(env%unit, '(10x,51("."))')
      write(env%unit, '(/,a)') "  cycle      energy         change     gnorm      step     trust   dxmax  linmrg"
   end if

   delta = delta0
   xold = mol%xyz
   blold = bl
   s = 0.0_wp
   have_step = .false.
   eold = etot
   dfpred = 0.0_wp
   prlevel = 0
   conv = .false.
   tm = 0.0_wp
   call system_clock(c0, crate)
   do iter = 1, maxopt
      ! energy and gradient
      g = 0.0_wp
      calc%accuracy = acc
      call calc%singlepoint(env, mol, chk, prlevel, iter == 1, energy, g, sigma, egap, res)
      call env%check(fail)
      if (fail) then
         call env%error('SCF not converged, aborting...', source)
         return
      end if
      call writeMolecule(mol, ilog, format=fileType%xyz, energy=res%e_total, &
         & gnorm=res%gnorm, number=iter)
      call system_clock(c1); tm(1) = tm(1) + real(c1 - c0, wp)/real(crate, wp); c0 = c1

      ! internal coordinates, gradient and basis at this geometry; the
      ! fragment rotations are re-referenced to the current geometry
      if (size(bl, 1) > 0) call fragment_rows(mol%xyz, frag, ncomp, bl)
      call bmat_all(prim, bl, mol%xyz, bq)
      call nonredundant_basis(bq, basis, binvt, nf)
      call eval_all(prim, bl, mol%xyz, q)
      gq = matmul(binvt, reshape(g, [n3]))
      gnorm = norm2(g)
      echng = energy - eold
      ! Hessian projected onto the current non-redundant basis
      hf = mm(basis, mm(h, basis), 'T', 'N')
      call system_clock(c1); tm(2) = tm(2) + real(c1 - c0, wp)/real(crate, wp); c0 = c1

      if (have_step) then
         ! secant pair from the step actually taken (fragment rows in the
         ! reference of the previous geometry, where they are exactly linear)
         sact(:np) = q(:np) - qold(:np)
         call wrap_dihedrals(prim, sact)
         if (nq > np) sact(np+1:) = matmul(blold, reshape(mol%xyz - xold, [n3]))
         call ts_bfgs_update(h, hf, basis, sact, gq - gqold)
         ! trust radius (Sella rules for minima)
         if (abs(dfpred) > 1.0e-14_wp*max(abs(eold), 1.0_wp)) then
            rho = echng/dfpred
            smax = maxval(abs(w*s))
            if (rho < 1.0_wp/rho_dec .or. rho > rho_dec) then
               delta = max(sigma_dec*smax, delta_min)
            else if (rho > 1.0_wp/rho_inc .and. rho < rho_inc) then
               delta = max(sigma_inc*smax, delta)
            end if
         end if
      end if

      call system_clock(c1); tm(3) = tm(3) + real(c1 - c0, wp)/real(crate, wp); c0 = c1
      if (pr) write(env%unit, '(i7,f18.10,es12.3,es10.2,es10.2,f9.4,2f8.3)') &
         & iter, energy, echng, gnorm, merge(maxval(abs(w*s)), 0.0_wp, have_step), delta, &
         & maxval(norm2(mol%xyz - xold, dim=1)), lin_margin(prim, mol%xyz)*180.0_wp/pi

      ! convergence (xtb criteria; energy must not have increased)
      if (have_step .and. abs(echng) < ethr .and. gnorm < gthr .and. echng < 1.0e-10_wp) then
         conv = .true.
         exit
      end if

      ! quasi-Newton step restricted by the largest internal-coordinate change
      call system_clock(c1); tm(6) = tm(6) + real(c1 - c0, wp)/real(crate, wp); c0 = c1
      call qn_step(h, hf, gq, basis, w, delta, s, dfpred)
      call system_clock(c1); tm(4) = tm(4) + real(c1 - c0, wp)/real(crate, wp); c0 = c1

      ! realise the step in Cartesian coordinates
      xold = mol%xyz
      qold = q
      gqold = gq
      eold = energy
      blold = bl
      call back_transform(prim, bl, mol%xyz, q, s, binvt, ok)
      ! limit the largest Cartesian atom displacement
      do i = 1, 3
         dmax = maxval(norm2(mol%xyz - xold, dim=1))
         if (dmax <= 1.05_wp*max_atom_displ) exit
         s = s*(max_atom_displ/dmax)
         mol%xyz = xold
         call back_transform(prim, bl, mol%xyz, q, s, binvt, ok)
      end do
      dfpred = dot_product(gq, s) + 0.5_wp*dot_product(s, matmul(h, s))
      call mol%update
      have_step = .true.
      call system_clock(c1); tm(5) = tm(5) + real(c1 - c0, wp)/real(crate, wp); c0 = c1

   end do

   etot = energy
   ! return the last evaluated geometry (its energy and gradient are in etot, g)
   if (.not.conv .and. have_step) then
      mol%xyz = xold
      call mol%update
   end if
   if (present(iter_needed)) iter_needed = min(iter, maxopt)
   if (pr) then
      if (conv) then
         write(env%unit, '(/,3x,a,i5,a)') "*** GEOMETRY OPTIMIZATION CONVERGED AFTER ", iter, " ITERATIONS ***"
      else
         write(env%unit, '(/,3x,a,i5,a)') "*** FAILED TO CONVERGE GEOMETRY OPTIMIZATION IN ", iter-1, " ITERATIONS ***"
      end if
   end if
   if (pr) write(env%unit, '(1x,a,6f9.2)') "intopt timings (s): SP, coords, update, step, backtr, misc:", tm
   fail = .not.conv
end subroutine intopt


!> Redundant primitives of the covalent graph and its connected components
!> (covalent fragments, ncomp of them)
subroutine build_primitives(at, xyz, prim, ncomp, frag, cn)
   integer, intent(in) :: at(:)
   real(wp), intent(in) :: xyz(:,:)
   type(redundant_type), intent(out) :: prim
   integer, intent(out) :: ncomp
   !> fragment index (1..ncomp) of every atom
   integer, allocatable, intent(out) :: frag(:)
   !> covalent coordination number of every atom
   integer, allocatable, intent(out), optional :: cn(:)
   type(TNeighbourList) :: neigh
   type(graph_type) :: graph
   integer :: nat, i, j, k
   integer, allocatable :: lab(:)

   nat = size(at)
   call init(neigh, nat)
   call neigh%generate_covalent(at, xyz)
   call init(graph, neigh)
   if (present(cn)) cn = graph%neighs(:nat)
   ! connected components (label propagation over the adjacency list)
   allocate(lab(nat))
   lab = [(i, i = 1, nat)]
   k = 1
   do while (k > 0)
      k = 0
      do i = 1, nat
         do j = 1, graph%neighs(i)
            if (lab(graph%i_neigh(j, i)) < lab(i)) then
               lab(i) = lab(graph%i_neigh(j, i)); k = k + 1
            else if (lab(graph%i_neigh(j, i)) > lab(i)) then
               lab(graph%i_neigh(j, i)) = lab(i); k = k + 1
            end if
         end do
      end do
   end do
   ncomp = 0
   allocate(frag(nat))
   do i = 1, nat
      if (lab(i) == i) then
         ncomp = ncomp + 1
         frag(i) = ncomp
      else
         frag(i) = frag(lab(i))
      end if
   end do
   call init(prim, graph, xyz)
end subroutine build_primitives


!> Orthonormal basis of the range of B (non-redundant subspace, nq x nf) and
!> binvt = (B^+)^T = G^+ B with G = B B^T (maps Cartesian gradients to internals)
subroutine nonredundant_basis(bq, basis, binvt, nf)
   real(wp), intent(in) :: bq(:,:)
   real(wp), allocatable, intent(out) :: basis(:,:)
   real(wp), allocatable, intent(out) :: binvt(:,:)
   integer, intent(out) :: nf
   real(wp), allocatable :: m(:,:), mw(:), v(:,:)
   integer :: n3, i, k, info
   real(wp) :: thr

   ! the non-zero spectrum of G = B B^T equals that of B^T B (3N x 3N):
   ! B^T B = V L V^T  =>  U = B V L^-1/2 spans range(B),  G^+ B = U L^-1/2 V^T
   n3 = size(bq, 2)
   m = mm(bq, bq, 'T', 'N')
   allocate(mw(n3))
   call symeig(m, mw, info)
   thr = 1.0e-8_wp*maxval(mw)
   nf = count(mw > thr)
   allocate(v(n3, nf))
   k = 0
   do i = 1, n3
      if (mw(i) > thr) then
         k = k + 1
         v(:, k) = m(:, i)/sqrt(mw(i))
      end if
   end do
   basis = mm(bq, v)                     ! nq x nf, orthonormal
   k = 0
   do i = 1, n3
      if (mw(i) > thr) then
         k = k + 1
         v(:, k) = m(:, i)
      end if
   end do
   k = 0
   deallocate(m)
   allocate(m(size(bq, 1), nf))
   do i = 1, n3
      if (mw(i) > thr) then
         k = k + 1
         m(:, k) = basis(:, k)/sqrt(mw(i))
      end if
   end do
   binvt = mm(m, v, 'N', 'T')            ! nq x 3N
end subroutine nonredundant_basis


!> C = op(A) op(B) via dgemm
function mm(a, b, ta, tb) result(c)
   real(wp), intent(in) :: a(:,:), b(:,:)
   character(len=1), intent(in), optional :: ta, tb
   real(wp), allocatable :: c(:,:)
   character(len=1) :: xa, xb
   integer :: m, n, k
   external :: dgemm
   xa = 'N'; xb = 'N'
   if (present(ta)) xa = ta
   if (present(tb)) xb = tb
   if (xa == 'N') then
      m = size(a, 1); k = size(a, 2)
   else
      m = size(a, 2); k = size(a, 1)
   end if
   if (xb == 'N') then
      n = size(b, 2)
   else
      n = size(b, 1)
   end if
   allocate(c(m, n))
   if (m*n == 0) return
   if (k == 0) then
      c = 0.0_wp
      return
   end if
   call dgemm(xa, xb, m, n, k, 1.0_wp, a, size(a, 1), b, size(b, 1), 0.0_wp, c, m)
end function mm


!> Quasi-Newton step s (nq) with level shift alpha >= 0 such that the largest
!> weighted internal-coordinate change does not exceed delta. dfpred is the
!> predicted energy change of the step.
subroutine qn_step(h, hfin, gq, basis, w, delta, s, dfpred)
   real(wp), intent(in) :: h(:,:), hfin(:,:), gq(:), basis(:,:), w(:), delta
   real(wp), intent(out) :: s(:), dfpred
   real(wp), allocatable :: hf(:,:), lam(:), vg(:), u(:,:)
   real(wp) :: alo, ahi, a
   integer :: nf, it, info

   nf = size(basis, 2)
   hf = hfin
   allocate(lam(nf))
   call symeig(hf, lam, info)
   lam = abs(lam)
   lam = max(lam, 1.0e-10_wp)
   vg = matmul(matmul(gq, basis), hf)
   u = mm(basis, hf)              ! nq x nf
   s = step_of(0.0_wp)
   if (maxval(abs(w*s)) > delta) then
      ! bisection on log(alpha): constraint value decreases with alpha
      alo = 0.0_wp
      ahi = max(1.0e-6_wp, maxval(lam))
      do while (maxval(abs(w*step_of(ahi))) > delta)
         alo = ahi
         ahi = 4.0_wp*ahi
      end do
      do it = 1, 100
         a = 0.5_wp*(alo + ahi)
         if (maxval(abs(w*step_of(a))) > delta) then
            alo = a
         else
            ahi = a
         end if
         if (ahi - alo < 1.0e-10_wp*max(1.0_wp, ahi)) exit
      end do
      s = step_of(ahi)
   end if
   dfpred = dot_product(gq, s) + 0.5_wp*dot_product(s, matmul(h, s))
contains
   function step_of(alpha) result(st)
      real(wp), intent(in) :: alpha
      real(wp) :: st(size(gq))
      st = -matmul(u, vg/(lam + alpha))
   end function step_of
end subroutine qn_step


!> TS-BFGS rank-two update for one secant pair (Sella, Bofill-type).
!> |H| is evaluated in the non-redundant subspace (basis), where s lies.
subroutine ts_bfgs_update(h, hf, basis, s, y)
   real(wp), intent(inout) :: h(:,:)
   !> projected Hessian basis^T h basis, updated consistently
   real(wp), intent(inout) :: hf(:,:)
   real(wp), intent(in) :: basis(:,:), s(:), y(:)
   real(wp), allocatable :: v(:,:), lam(:), j(:), absbs(:), x(:), u(:), pu(:), pj(:)
   real(wp) :: xs, jts
   integer :: n, info, k

   n = size(s)
   if (norm2(s) < 1.0e-8_wp) return
   v = hf
   allocate(lam(size(v, 1)))
   call symeig(v, lam, info)
   if (info /= 0) return
   j = y - matmul(h, s)
   absbs = matmul(basis, matmul(v, abs(lam)*matmul(matmul(s, basis), v)))
   x = dot_product(s, y)*y + dot_product(s, absbs)*absbs
   xs = dot_product(x, s)
   if (abs(xs) < tiny(1.0_wp)) then
      ! PSB fallback
      u = s/dot_product(s, s)
   else
      u = x/xs
   end if
   jts = dot_product(j, s)
   do k = 1, n
      h(:, k) = h(:, k) + u*j(k) + j*u(k) - jts*u*u(k)
   end do
   pu = matmul(u, basis); pj = matmul(j, basis)
   do k = 1, size(hf, 2)
      hf(:, k) = hf(:, k) + pu*pj(k) + pj*pu(k) - jts*pu*pu(k)
   end do
end subroutine ts_bfgs_update


!> Iterative back-transformation: find x with q(x) = q0 + s (chord iterations
!> with the generalised inverse at the start geometry). On failure the
!> first-order step is kept.
subroutine back_transform(prim, bl, xyz, q0, s, binvt, ok)
   type(redundant_type), intent(in) :: prim
   real(wp), intent(in) :: bl(:,:)
   real(wp), intent(inout) :: xyz(:,:)
   real(wp), intent(in) :: q0(:), s(:), binvt(:,:)
   logical, intent(out) :: ok
   real(wp), allocatable :: x(:), x1(:), qt(:), q(:), r(:), dx(:)
   real(wp) :: rnorm, rold
   integer :: it, n3

   n3 = size(xyz)
   x = reshape(xyz, [n3])
   qt = q0 + s
   allocate(q(size(q0)))
   dx = matmul(s, binvt)
   x = x + dx
   x1 = x
   ok = .false.
   rold = huge(1.0_wp)
   do it = 1, 50
      call eval_all(prim, bl, reshape(x, shape(xyz)), q)
      r = qt - q
      call wrap_dihedrals(prim, r)
      rnorm = norm2(r)
      if (rnorm > rold .or. rnorm /= rnorm) exit
      rold = rnorm
      dx = matmul(r, binvt)
      x = x + dx
      if (maxval(abs(dx)) < 1.0e-7_wp) then
         ok = .true.
         exit
      end if
   end do
   if (.not.ok) then
      ! keep the last improving iterate if it is sane, else the first-order step
      if (rold < huge(1.0_wp) .and. all(x == x)) then
         ok = .true.
      else
         x = x1
      end if
   end if
   xyz = reshape(x, shape(xyz))
end subroutine back_transform


!> Fischer-Almloef diagonal guess (as in Sella), atomic units
subroutine fa_guess(prim, at, xyz, h0)
   type(redundant_type), intent(in) :: prim
   integer, intent(in) :: at(:)
   real(wp), intent(in) :: xyz(:,:)
   real(wp), intent(out) :: h0(:)
   integer, allocatable :: nb(:)
   integer :: ic, i, j, k, l
   real(wp) :: rab, rbc, ra, rb, rc, rcab, rcbc

   allocate(nb(size(at)), source=0)
   do ic = 1, prim%ncoords
      if (prim%kind(ic) /= coord_bond) cycle
      nb(prim%atoms(1, ic)) = nb(prim%atoms(1, ic)) + 1
      nb(prim%atoms(2, ic)) = nb(prim%atoms(2, ic)) + 1
   end do
   do ic = 1, prim%ncoords
      i = prim%atoms(1, ic); j = prim%atoms(2, ic)
      select case (prim%kind(ic))
      case (coord_bond)
         rab = norm2(xyz(:, i) - xyz(:, j))
         h0(ic) = 0.3601_wp*exp(-1.944_wp*(rab - rcov(at(i)) - rcov(at(j))))
      case (coord_angle, coord_linbend)
         k = prim%atoms(3, ic)
         rab = norm2(xyz(:, i) - xyz(:, j)); rbc = norm2(xyz(:, k) - xyz(:, j))
         rcab = rcov(at(i)) + rcov(at(j)); rcbc = rcov(at(k)) + rcov(at(j))
         h0(ic) = 0.089_wp + 0.11_wp*exp(-0.44_wp*(rab + rbc - rcab - rcbc)) &
            & / (rcab*rcbc)**(-0.42_wp)
      case (coord_dihedral)
         k = prim%atoms(3, ic)
         rbc = norm2(xyz(:, j) - xyz(:, k))
         rcbc = rcov(at(j)) + rcov(at(k))
         l = nb(j) + nb(k) - 2
         h0(ic) = 0.0015_wp + 14.0_wp*real(max(l, 0), wp)**0.57_wp &
            & *exp(-2.85_wp*(rbc - rcbc))/(rbc*rcbc)**4.0_wp
      case default
         h0(ic) = 0.05_wp
      end select
   end do
contains
   !> unscaled covalent radius (D3 radii are Pyykkoe radii scaled by 4/3), bohr
   elemental real(wp) function rcov(z)
      integer, intent(in) :: z
      rcov = 0.75_wp*covalent_radius_d3(z)
   end function rcov
end subroutine fa_guess


!> Values of the primitives (same conventions as xtb_internals_redundant)
!> Main-group 'organic' elements (H, B-F, Si-Cl, Ge-Br, Se, Sn-I)
elemental logical function organic(z)
   integer, intent(in) :: z
   organic = z == 1 .or. (z >= 5 .and. z <= 9) .or. (z >= 14 .and. z <= 17) &
      & .or. (z >= 32 .and. z <= 35) .or. (z >= 50 .and. z <= 53)
end function organic


!> Translation and (linearised) rotation coordinates of every fragment,
!> referenced to the geometry xyz (TRIC-like, Wang & Song, JCP 144, 214108
!> (2016)), as constant B-matrix rows: q = bl x.
!>   T_f = (1/n_f) sum_i x_i
!>   R_f = r_g I^+ sum_i d_i x x_i,  d_i = x_i^0 - c_f^0,
!>   I = sum_i (|d_i|^2 1 - d_i d_i^T),  r_g = sqrt(sum_i |d_i|^2 / n_f)
!> For a rigid rotation x_i = x_i^0 + theta x d_i this gives R_f = r_g theta
!> (length units); R_f is translation invariant. One-atom fragments only get
!> translations, linear fragments a rank-two rotation block.
subroutine fragment_rows(xyz, frag, nfrag, bl)
   real(wp), intent(in) :: xyz(:,:)
   integer, intent(in) :: frag(:), nfrag
   real(wp), allocatable, intent(out) :: bl(:,:)
   integer :: nat, f, i, a, k, nrow, nf, info, row
   real(wp) :: c(3), d(3), im(3,3), ev(3), ip(3,3), rg, e(3), cr(3)

   nat = size(frag)
   nrow = 0
   do f = 1, nfrag
      nrow = nrow + merge(6, 3, count(frag == f) > 1)
   end do
   allocate(bl(nrow, 3*nat), source=0.0_wp)
   row = 0
   do f = 1, nfrag
      nf = count(frag == f)
      c = 0.0_wp
      do i = 1, nat
         if (frag(i) == f) c = c + xyz(:, i)
      end do
      c = c/real(nf, wp)
      do a = 1, 3
         do i = 1, nat
            if (frag(i) == f) bl(row + a, 3*(i-1) + a) = 1.0_wp/real(nf, wp)
         end do
      end do
      row = row + 3
      if (nf == 1) cycle
      im = 0.0_wp
      rg = 0.0_wp
      do i = 1, nat
         if (frag(i) /= f) cycle
         d = xyz(:, i) - c
         rg = rg + dot_product(d, d)
         do a = 1, 3
            im(a, a) = im(a, a) + dot_product(d, d)
         end do
         im = im - spread(d, 2, 3)*spread(d, 1, 3)
      end do
      rg = sqrt(rg/real(nf, wp))
      ! pseudo-inverse of the inertia-like tensor
      call symeig(im, ev, info)
      ip = 0.0_wp
      do k = 1, 3
         if (ev(k) > 1.0e-8_wp*maxval(ev)) &
            & ip = ip + spread(im(:, k), 2, 3)*spread(im(:, k), 1, 3)/ev(k)
      end do
      ! d/dx_i,b of (d_i x x_i)_g = (d_i x e_b)_g
      do i = 1, nat
         if (frag(i) /= f) cycle
         d = xyz(:, i) - c
         do k = 1, 3
            e = 0.0_wp; e(k) = 1.0_wp
            cr = cross(d, e)
            bl(row+1:row+3, 3*(i-1) + k) = rg*matmul(ip, cr)
         end do
      end do
      row = row + 3
   end do
end subroutine fragment_rows


!> Smallest distance (rad) of any bending angle, or any interior angle of a
!> dihedral, from 0 or pi; dihedrals and angles become singular there.
pure function lin_margin(prim, xyz) result(m)
   type(redundant_type), intent(in) :: prim
   real(wp), intent(in) :: xyz(:,:)
   real(wp) :: m
   integer :: ic
   m = pi
   do ic = 1, prim%ncoords
      select case (prim%kind(ic))
      case (coord_angle)
         m = min(m, margin(prim%atoms(1, ic), prim%atoms(2, ic), prim%atoms(3, ic)))
      case (coord_dihedral)
         m = min(m, margin(prim%atoms(1, ic), prim%atoms(2, ic), prim%atoms(3, ic)), &
            & margin(prim%atoms(2, ic), prim%atoms(3, ic), prim%atoms(4, ic)))
      end select
   end do
contains
   pure real(wp) function margin(i, j, k)
      integer, intent(in) :: i, j, k
      real(wp) :: v1(3), v2(3), a
      v1 = xyz(:, i) - xyz(:, j); v2 = xyz(:, k) - xyz(:, j)
      a = acos(max(-1.0_wp, min(1.0_wp, dot_product(v1, v2)/(norm2(v1)*norm2(v2)))))
      margin = min(a, pi - a)
   end function margin
end function lin_margin


pure subroutine eval_internals(prim, xyz, q)
   type(redundant_type), intent(in) :: prim
   real(wp), intent(in) :: xyz(:,:)
   real(wp), intent(out) :: q(:)
   integer :: ic, i, j, k, l
   real(wp) :: v1(3), v2(3), b1(3), b2(3), b3(3), n1(3), n2(3), m(3), co

   do ic = 1, prim%ncoords
      i = prim%atoms(1, ic); j = prim%atoms(2, ic)
      select case (prim%kind(ic))
      case (coord_bond)
         q(ic) = norm2(xyz(:, i) - xyz(:, j))
      case (coord_angle)
         k = prim%atoms(3, ic)
         v1 = xyz(:, i) - xyz(:, j); v2 = xyz(:, k) - xyz(:, j)
         co = dot_product(v1, v2)/(norm2(v1)*norm2(v2))
         q(ic) = acos(max(-1.0_wp, min(1.0_wp, co)))
      case (coord_linbend)
         k = prim%atoms(3, ic)
         v1 = xyz(:, i) - xyz(:, j); v1 = v1/norm2(v1)
         v2 = xyz(:, k) - xyz(:, j); v2 = v2/norm2(v2)
         q(ic) = dot_product(prim%frame(:, ic), v1 + v2)
      case (coord_dihedral)
         k = prim%atoms(3, ic); l = prim%atoms(4, ic)
         b1 = xyz(:, j) - xyz(:, i)
         b2 = xyz(:, k) - xyz(:, j)
         b3 = xyz(:, l) - xyz(:, k)
         n1 = cross(b1, b2)
         n2 = cross(b2, b3)
         m = cross(n1, b2/norm2(b2))
         q(ic) = atan2(dot_product(m, n2), dot_product(n1, n2))
      end select
   end do
end subroutine eval_internals


pure function cross(a, b) result(c)
   real(wp), intent(in) :: a(3), b(3)
   real(wp) :: c(3)
   c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
end function cross


!> All coordinates: primitives followed by the linear fragment rows
subroutine eval_all(prim, bl, xyz, q)
   type(redundant_type), intent(in) :: prim
   real(wp), intent(in) :: bl(:,:), xyz(:,:)
   real(wp), intent(out) :: q(:)
   integer :: np
   np = prim%ncoords
   if (np > 0) call eval_internals(prim, xyz, q(:np))
   if (size(bl, 1) > 0) q(np+1:) = matmul(bl, reshape(xyz, [size(xyz)]))
end subroutine eval_all


subroutine bmat_all(prim, bl, xyz, bq)
   type(redundant_type), intent(in) :: prim
   real(wp), intent(in) :: bl(:,:), xyz(:,:)
   real(wp), intent(out) :: bq(:,:)
   integer :: np
   np = prim%ncoords
   if (np > 0) call get_bmatrix(prim, xyz, bq(:np, :))
   if (size(bl, 1) > 0) bq(np+1:, :) = bl
end subroutine bmat_all


pure subroutine wrap_dihedrals(prim, dq)
   type(redundant_type), intent(in) :: prim
   real(wp), intent(inout) :: dq(:)
   integer :: ic
   do ic = 1, prim%ncoords
      if (prim%kind(ic) == coord_dihedral) dq(ic) = dq(ic) - 2.0_wp*pi*anint(dq(ic)/(2.0_wp*pi))
   end do
end subroutine wrap_dihedrals


subroutine symeig(a, w, info)
   real(wp), intent(inout) :: a(:,:)
   real(wp), intent(out) :: w(:)
   integer, intent(out) :: info
   real(wp), allocatable :: work(:)
   integer, allocatable :: iwork(:)
   integer :: n, lw, liw
   external :: dsyevd
   n = size(a, 1)
   allocate(work(1), iwork(1))
   call dsyevd('V', 'U', n, a, n, w, work, -1, iwork, -1, info)
   lw = int(work(1)); liw = iwork(1)
   deallocate(work, iwork); allocate(work(lw), iwork(liw))
   call dsyevd('V', 'U', n, a, n, w, work, lw, iwork, liw, info)
end subroutine symeig

end module xtb_optimizer_internal
