!  Author : Pacidus
! Started at: 31.10.2020

Module Param
  Implicit None
  Integer(1), Parameter :: D      = 2                       ! D2
  Integer(1), Parameter :: Q      = 9                       ! Q9
  Real(8), Parameter :: c         = 340._8                  ! m·s⁻¹
  Real(8), Parameter :: lx        = 4._8                    ! m
  Real(8), Parameter :: ly        = 2._8                    ! m
  Real(8), Parameter :: dl        = 1D-2                    ! m
  Real(8), Parameter :: dt        = dl/c                    ! s
  Real(8), Parameter :: t         = 10._8                   ! s
  Real(8), Parameter :: o1        = 3._8/c                  ! s·m⁻¹
  Real(8), Parameter :: o2        = 9._8/(2*c*c)            ! s²·m⁻²
  Real(8), Parameter :: o3        = 3._8/(2*c*c)            ! s²·m⁻²
  Real(8), Parameter :: rho       = 1.177_8                 ! kg·m⁻³
  Real(8), Parameter :: nu        = 1.57D-5                 ! m²·s⁻¹
  Real(8), Parameter :: tau       = (nu/(c*c*dt)) + 5D-1    ! ø
  Real(8), Parameter :: itau      = 1._8/tau                ! ø
  Real(8), Parameter :: mitau     = 1 - itau                ! ø
  Integer(4), Parameter :: H      = int(ly/dl)              ! ø
  Integer(4), Parameter :: L      = int(lx/dl)              ! ø
  Integer(4), Parameter :: P      = int(t/dt)               ! ø
  Integer(1), Parameter :: e(Q,D) = transpose(reshape(&     ! ø
    (/0, 0,&
      1, 0,&
      0, 1,&
     -1, 0,&
      0,-1,&
      1, 1,&
     -1, 1,&
     -1,-1,&
      1,-1/), (/D,Q/)))
  Real(8), Parameter :: w(Q) = &                            ! ø
    (/4._8/9._8,&
      1._8/9._8,&
      1._8/9._8,&
      1._8/9._8,&
      1._8/9._8,&
      1._8/36._8,&
      1._8/36._8,&
      1._8/36._8,&
      1._8/36._8/)
End Module Param

Module Functions
  
  Contains

    Function Object()
      Use Param
      Implicit None
      Integer(1) :: Object(L,H)
      Character(Len = 140), Parameter  :: Parbash = '(&
        "convert -background white -compress None objet.svg -resize ",&
        i0,&
        "x",&
        i0,&
        "^ -gravity center -extent ",i0,"x",i0," -threshold 99% objet.pbm")'
      Character(len=140) :: command
      Character(LEN=20) :: trash
    	Integer :: i
      
      Write(command, Parbash) L, H, L, H
      Call execute_command_line(command, wait=.true.)
    
      Open(10, file = 'objet.pbm', status = 'old', action='read')
      Read(10, *) trash
      Read(10, *) trash
    
      Do i=1, H
        Read(10, *) Object(:, i)
    	End Do
    	
    	Close(10)
    End Function

End Module Functions

Module Subroutines

  Contains
  
    subroutine InitF(F)
      Use Param
      Implicit None
      Real(8) :: F(Q,L,H)
      Integer(1) :: i
      
      F = rho
      
      Do i=1, Q
          F(i,:,:) = F(i,:,:)*w(i) 
      End Do
      
    end subroutine
    
    subroutine Stream(F)
      Use Param
      Implicit None
      Real(8) :: F(Q,L,H)
      Integer(1) :: i, j

      Do i=1, Q
        Do j=1, D
          F(i,:,:) = Cshift(F(i,:,:), e(i,j),j)
        End Do
      End Do

    end subroutine
    
End Module Subroutines

Module Var
  Use Param
  Implicit None
  Integer(1) :: Obj(L, H)
  Real(8) :: F(Q,L,H)
End Module Var

Program LBM2D
  Use Param
  Use Functions
  Use Subroutines
  Use Var
  Implicit None
  Integer(8) :: dP = 0
  
  Call InitF(F)

  Obj = Object()
  Do dP = 0, P  
    Call Stream(F)
  End Do
End Program LBM2D
