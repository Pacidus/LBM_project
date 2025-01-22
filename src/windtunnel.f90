! Author : Pacidus
! Started at: 31.10.2020 00:17:57
! Rework Sunday 19.01.2025 23:33

!===============================================================================
!                               Fluid Constants.
!===============================================================================
Module Fconst
  Implicit None
  Real(8), Parameter :: cs  = 3.4d2         ! m·s⁻¹   Speed of sound
  !Real(8), Parameter :: rho = 1.177d0      ! kg·m⁻³  Fluid density
  Real(8), Parameter :: rho = 3.639d-1      ! kg·m⁻³  Fluid density troposphere
  !Real(8), Parameter :: nu  = 1            ! m²·s⁻¹  Kinematic viscosity
  Real(8), Parameter :: nu  = 3.5d-5        ! m²·s⁻¹  Kinematic viscosity tropo
End Module Fconst


!===============================================================================
!                            Simulation Constants.
!===============================================================================
Module Sconst
  Implicit None
  Real(8), Parameter :: lx  = 2.5d0         ! m     Length of the box
  Real(8), Parameter :: ly  = 2d0           ! m     Height of the box
  Real(8), Parameter :: t   = 5d0           ! s     Total time of the simulation
  Real(8), Parameter :: sp  = 1.5d1         ! m·s⁻¹ Inlet Speed
End Module Sconst


!===============================================================================
!                         Discretisation Constants.
!===============================================================================
Module Dconst
  Use Fconst, Only: cs
  Use Sconst, Only: lx, ly, t
  Implicit None
  Real(8), Parameter :: dt    = 1d-6      ! s Timestep
  Real(8), Parameter :: dl    = cs*dt       ! m Spatial step
  Integer(4), Parameter :: H  = Int(ly/dl)  ! ø Number of height-step
  Integer(4), Parameter :: L  = Int(lx/dl)  ! ø Number of length-step
  Integer(4), Parameter :: NP = Int(t/dt)   ! ø Number of timestep (partial)
  Integer(4), Parameter :: st = 200         ! ø Number of snapshot
  integer(4), Parameter :: d  = NP/st       ! ø Itterations between of snapshot
  integer(4), Parameter :: ds = 10          ! ø space between saved point
  Integer(4), Parameter :: N  = st*(d + 1)  ! ø Number of timestep (total)
End Module Dconst

!===============================================================================
!                         Computationnal Constants.
!===============================================================================
Module Cconst
  Use Fconst, Only: cs, nu
  Use Dconst, Only: dt
  Use Sconst, Only: sp
  Implicit None
  Real(8), Parameter :: sqrt3 = Dsqrt(3d0)  ! ø       No reason except aesthetic
  Real(8), Parameter :: c     = cs*sqrt3    ! m·s⁻¹   Lattice speed
  Real(8), Parameter :: cs2   = cs*cs       ! m²·s⁻²  Sound speed squared
  Real(8), Parameter :: ics2  = 1d0/cs2     ! s²·m⁻²  Inverse of cs2
  Real(8), Parameter :: clean = nu*ics2/dt  ! ø       No reason except aesthetic
  Real(8), Parameter :: tau   = clean+5d-1  ! ø       Characteristic timescale
  Real(8), Parameter :: itau  = 1d0 / tau   ! ø       Inverse of tau
  Real(8), Parameter :: mitau = 1d0 - itau  ! ø       One minus inverse of tau
  Real(8), Parameter :: zh    = cs/(cs-sp)  ! ø       Zou/He parameter
  Real(8), Parameter :: sp2   = sp*sp       ! m²·s⁻²  Outlet speed squared
End Module Cconst


!===============================================================================
!                             Lattice Constants.
!===============================================================================
Module Lconst
  Use Cconst, Only: c
  Implicit None
  Integer(1), Parameter :: D        = 2     ! D2  Number of spatial dimension
  Integer(1), Parameter :: Q        = 9     ! Q9  Number of speed discretization
  Integer(8), Parameter :: ed(Q,D)  = &     ! ø   Directions of the speeds
    reshape(&
    [0,  1,  0, -1,  0,  1, -1, -1,  1,&
     0,  0,  1,  0, -1,  1,  1, -1, -1], [Q,D])
  Real(8), Parameter :: e(Q,D) = c*dble(ed) ! m·s⁻¹   Speeds vectors
  Real(8), Parameter :: w(Q)   = &          ! ø       D2Q9 Weights
    [4d0, 1d0, 1d0, 1d0, 1d0, .25d0, .25d0, .25d0, .25d0]/9d0
End Module Lconst

!===============================================================================
!                           Simulation Variables.
!===============================================================================
Module Var
  Implicit None
  Logical, Allocatable :: Obj(:, :)         ! ø       Mask of the Object
  Logical, Allocatable :: nObj(:, :)        ! ø       Mask of the Empty space
  Real(8), Allocatable :: Feq(:,:,:)        ! kg·m⁻³  Density speed distribution
  Real(8), Allocatable :: Usqr(:,:)         ! m·s⁻¹   Macroscopic speed squared
  Real(8), Allocatable :: F(:,:,:)          ! kg·m⁻³  Density speed distribution
  Real(8), Allocatable :: U(:,:,:)          ! m·s⁻¹   Macroscopic speed
  Real(8), Allocatable :: Rho(:,:)          ! kg·m⁻³  Macroscopic density
End Module Var


!===============================================================================
!                               The main program.
!===============================================================================
Program LBM2D
  Use Var
  Use Dconst, Only: L, H, N, d, st
  Implicit None
  Integer(4) :: dN, ni = 0
  Call Allocatall
  Call InitF

  Do dN = 0, N
    If (Modulo(dN, d) == 0) Then
      Call Savefile(ni)
      Write(*, *) ni, "/", st
      ni = ni + 1
    End If

    CALL IOlet
    CALL Bound
    CALL CMacro
    CALL CFeq
    CALL Collide
    CALL Stream
  End Do

  Call Deallocatall
End Program LBM2D


!===============================================================================
!                               Subroutines
!===============================================================================

!===============================================================================
!                             Inlet & Outlet
!           Using Zou/He boundary condition to implement Dirichlet
!===============================================================================
Subroutine IOlet
  Use Cconst, Only: cs, ics2, zh, sp2
  Use Lconst, Only: Q, e, w
  Use Var, Only: F, Rho, U
  Use Dconst, Only: L, H
  Use Sconst, Only: sp
  Implicit None
  Real(8) :: sca(H), Feq1(H), Feq2(H), R(H)
  Integer(4) :: j!, i

  F([7,4,8],L,:) = F([7,4,8],L-1,:)

  ! Inflow condition Left wall
  sca(:) = (F(4,1,:)+F(1,1,:)+F(2,1,:)+2*(F(7,1,:)+F(4,1,:)+F(8,1,:)))
  R(:) = sca*zh
  Rho(1,:) = R(:)

  sca(:) = e(8,1)*sp
  Feq1 = w(8)*(1d0+(sca+.5d0*((sca*sca*ics2)-sp2))*ics2)
  sca = e(6,1)*sp
  Feq2 = w(6)*(1d0+(sca+.5d0*((sca*sca*ics2)-sp2))*ics2)
  F(6,1,:) = F(8,1,:) + r*(Feq2 - Feq1)

  sca = e(4,1)*sp
  Feq1 = w(4)*(1d0+(sca+.5d0*((sca*sca*ics2)-sp2))*ics2)
  sca = e(2,1)*sp
  Feq2 = w(2)*(1d0+(sca+.5d0*((sca*sca*ics2)-sp2))*ics2)
  F(2,1,:) = F(4,1,:) + r*(Feq2 - Feq1)

  sca = e(7,1)*sp
  Feq1 = w(7)*(1d0+(sca+.5d0*((sca*sca*ics2)-sp2))*ics2)
  sca = e(9,1)*sp
  Feq2 = w(9)*(1d0+(sca+.5d0*((sca*sca*ics2)-sp2))*ics2)
  F(9,1,:) = F(7,1,:) + r*(Feq2 - Feq1)
End Subroutine


!===============================================================================
!                           Boundary with collision
!===============================================================================
Subroutine Bound
  Use Dconst, Only: L, H
  Use Var, Only: F, Obj
  Implicit None
  Real(8) :: temp(L, H)

  Where(Obj)
    temp = F(4, :, :)
    F(4, :, :) = F(2, :, :)
    F(2, :, :) = temp

    temp = F(5, :, :)
    F(5, :, :) = F(3, :, :)
    F(3, :, :) = temp

    temp = F(8, :, :)
    F(8, :, :) = F(6, :, :)
    F(6, :, :) = temp

    temp = F(9, :, :)
    F(9, :, :) = F(7, :, :)
    F(7, :, :) = temp
  End Where
End Subroutine

!===============================================================================
!                               Compute Macro
!===============================================================================
Subroutine CMacro
  Use Var, Only: F, Rho, U, Usqr
  Use Dconst, Only: L, H
  Use Lconst, Only: Q, e
  Implicit None
  Integer(1) :: k
  Real(8) :: R(L, H)
  Rho = 0
  U = 0
  Do k = 1, Q
    Rho(:,:) = Rho(:,:) + F(k,:,:)
    If(e(k,1) /= 0) Then
      U(1,:,:) = U(1,:,:) + F(k,:,:) * e(k,1)
    End If
    If(e(k,2) /= 0) Then
      U(2,:,:) = U(2,:,:) + F(k,:,:) * e(k,2)
    End If
  End Do
  U(1,:,:)  = U(1,:,:) / Rho(:,:)
  U(2,:,:)  = U(2,:,:) / Rho(:,:)
  Usqr(:,:) = U(1,:,:) * U(1,:,:) + U(2,:,:) * U(2,:,:)
End Subroutine


!===============================================================================
!                               Compute Feq
!===============================================================================
Subroutine CFeq
  Use Var, Only: Feq, U, Usqr, Rho, nObj
  Use Lconst, Only: Q, e, w
  Use Dconst, Only: L, H
  Use Cconst, Only: ics2
  Implicit None
  Real(8) :: Ue(L,H)
  Integer(1) :: k
  
  Feq(1,:,:) = w(1) * Rho(:,:) * ( 1d0 - ( .5d0 * Usqr(:,:) * ics2))
  
  Do k = 2, Q
      Ue = e(k,1)*U(1,:,:)+e(k,2)*U(2,:,:)
      Feq(k,:,:) = w(k)*Rho(:,:)*(1d0+(Ue(:,:)+.5d0*((Ue(:,:)*Ue(:,:)*ics2)-Usqr(:,:)))*ics2)
  End Do
 End Subroutine


!===============================================================================
!                             Collide Step
!===============================================================================
Subroutine Collide
  Use Var, Only: F, Feq, nObj
  Use Cconst, Only: itau, mitau
  Use Dconst, Only: L, H
  Use Lconst, Only: Q
  Implicit None
  Integer(8) :: k
  Do k=1, Q
      F(k, :, :) = mitau * F(k, :, :) + itau * Feq(k, :, :)
  End Do
End Subroutine


!===============================================================================
!                             Stream Step
!===============================================================================
Subroutine Stream
  Use Lconst, Only: Q, ed
  Use Dconst, Only: L, H
  Use Var, Only: F, nObj
  Implicit None
  Integer(8) :: i, j
  Integer(1) :: k
  
  Do k = 2, Q
    If (ed(k,1) == 1) Then
      F(k,2:L,:) = F(k,:L-1,:)
    Else If (ed(k,1) == -1) Then
      F(k,:L-1,:) = F(k,2:L,:)
    End If
    If (ed(k,2) == 1) Then
      F(k,:,2:H) = F(k,:,:H-1)
    Else If (ed(k, 2) == -1) Then
      F(k,:,1:H-1) = F(k,:,2:H)
    End If
  End Do 
End Subroutine

!===============================================================================
!                           Initialisation of F
!===============================================================================
Subroutine InitF
  Use var, Only: U, Rho, F, Feq, Usqr
  Use Cconst, Only: sp2
  Use Sconst, Only: sp
  Rho = 1
  U(1,:,:) = sp
  U(2,:,:) = 0
  Usqr = sp2
  Call CFeq
  F = Feq
End Subroutine


!===============================================================================
!              Convert the image of the object from svg to pgm
!     and import it into a matrix with 1 and 0 (1 stand for the object)
!===============================================================================
Subroutine Object
  Use Dconst, Only: L, H
  Use Var, Only: Obj, nObj
  Implicit None
  Character(len=210) :: command             ! The bash command
  Integer(1) :: Ob(L,H)
  Integer(4) :: i, j                        ! Itterator

  command = '("rsvg-convert -w ",i0," -h ",i0," objet.svg -o objet.png")'
  Write(command, command) L, H

  Call execute_command_line(command, wait=.True.)

  command = "magick -compress none objet.png -alpha extract -threshold 0%&
  & -negate objet.pbm"
  Call execute_command_line(command, wait=.true.)

  Open(10, file = 'objet.pbm', action='read')

  Read(10, *) command
  Read(10, *) command

  Do i=1, H
    Read(10, *) Ob(:, i)
  End Do
  Close(10)
  Allocate(Obj(L,H))
  Allocate(nObj(L,H))
  Obj = .False.
  nObj = (Ob == 0)
  Open(10, file = 'nObj', action='write')
  Open(11, file = 'Obj', action='write')
  Do i=1, L
    Do j=1, H
      If( Ob(i,j) == 1 ) Then
        Obj(i,j) = ANY(Ob(i-1:i+1,j-1:j+1) == 0)
        nObj(i,j) = Obj(i,j)
      End If
    End Do
    Write(10,*) nObj(i,:)
    Write(11,*) Obj(i,:)
  End Do
  Close(10)
  Close(11)

End Subroutine


!===============================================================================
!              Save the macroscopic values such as speed and density
!===============================================================================
Subroutine Savefile(dN)
  Use Dconst, Only: ds
  Use Var, Only: U, Rho
  Use Dconst, Only: H
  Implicit None
  Character(len=24) :: Namefile
  Integer(4) :: dN, i

  Write(Namefile, '("./U/",i5.5,".bin")') dN
  Open(10, file=Namefile, action="Write", form="unformatted")
  Write(10) U(::ds,::ds,:)
  Close(10)

  Write(Namefile, '("./P/",i5.5,".bin")')  dN
  Open(10, file=Namefile, action="Write", form="unformatted")
  Write(10) Rho(::ds,::ds)
  Close(10)
End Subroutine



!===============================================================================
!                             Allocate all matrixes
!===============================================================================
Subroutine Allocatall
  Use Var, Only: F, Feq, U, Usqr, Rho
  Use Dconst, Only: L, H
  Use Lconst, Only: Q, D

  Allocate(F(Q,L,H))
  Allocate(Feq(Q,L,H))
  Allocate(U(D,L,H))
  Allocate(Usqr(L,H))
  Allocate(Rho(L,H))
  Call Object
End Subroutine


!===============================================================================
!                           Deallocate all matrixes
!===============================================================================
Subroutine Deallocatall
  Use Var

  Deallocate(F)
  Deallocate(Feq)
  Deallocate(U)
  Deallocate(Usqr)
  Deallocate(Rho)
  Deallocate(Obj)
End Subroutine
