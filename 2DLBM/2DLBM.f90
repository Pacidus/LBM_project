!  Author : Pacidus
! Started at: 31.10.2020 00:17:57

Module Param
  ! Initializes all the constants
  Implicit None
  Integer(1), Parameter :: D      = 2                       ! D2
  Integer(1), Parameter :: Q      = 9                       ! Q9
  Real(4), Parameter :: lx        = 4._8                    ! m
  Real(4), Parameter :: ly        = 2._8                    ! m
  Real(4), Parameter :: dl        = 10D-3                    ! m
  Real(4), Parameter :: dt        = 1D-3                    ! s
  Real(4), Parameter :: c         = dl/dt                   ! m·s⁻¹
  Real(4), Parameter :: t         = 1._8                    ! s
  Real(4), Parameter :: rho       = 1.177_8                 ! kg·m⁻³
  Real(4), Parameter :: nu        = 1.57D-5                 ! m²·s⁻¹
  Real(4), Parameter :: tau       = (nu/(c*c*dt)) + dt/2._8 ! ø
  Real(4), Parameter :: itau      = 1._8/tau                ! ø
  Real(4), Parameter :: mitau     = 1._8 - itau             ! ø
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
  Integer(1), Parameter :: Me(Q,Q) = matmul(e,transpose(e))
  Real(4), Parameter :: w(Q) = &                            ! ø
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
  ! Initializes all the functions
  Implicit None
  Contains
  
    Function R(F)
      Use Param
      Implicit None
      Real(4) :: R(L, H)
      Real(4) :: F(Q, L, H)
      R = Sum(F, 1)
    End Function

    Function eu(F)
      Use Param
      Implicit None
      Real(4) :: eu(Q, L, H)
      Real(4) :: F(Q, L, H)
      Integer(1) :: i, j
      
      eu = 0
      
      !Omp Parallel Do
        Do i = 2, Q
          Do j = 1, Q
            If (Me(i, j) /= 0) Then
              eu(i, :, :) = eu(i, :, :) + F(j, :, :)*Me(i, j)
            End If
          End Do
        End Do
      !Omp End Parallel Do
    End Function
     
    Function uu(F)
      Use Param
      Implicit None
      Real(4) :: uu(L, H)
      Real(4) :: u1(L, H) = 0
      Real(4) :: u2(L, H) = 0
      Real(4) :: Fi(L, H)
      Real(4) :: F(Q, L, H)
      Integer(1) :: i

      Do i = 2, Q
        Fi = F(i, :, :)
        If (e(i, 1) /= 0) Then
          u1 = u1 + Fi*e(i, 1)
        End If
        If (e(i, 2) /= 0) Then
          u2 = u2 + Fi*e(i, 2)
        End If
      End Do

      uu = u1*u1+u2*u2
     End Function
     
    Function Feq(F)
      Use Param
      Implicit None
      Real(4) :: Feq(Q, L, H)
      Real(4) :: F(Q, L, H)
      Real(4) :: steu(Q, L, H)
      Real(4) ::  st2eu(Q, L, H)
      Real(4) :: stuu(L, H)
      Real(4) :: stR(L, H)
      Integer(1) :: i
      
      stuu = uu(F)
      steu = eu(F)
      st2eu = steu*steu
      stR = R(F)
      
      !$Omp Parallel Do
        Do i = 1, Q
          Feq(i,:,:) = w(i)*(3*(steu(i,:,:) + (3*st2eu(i,:,:) - stuu*5D-1)/stR) + stR)
        End Do
      !$Omp End Parallel Do
     End Function
End Module Functions

Module Subroutines
  ! Initializes all the Subroutines
  Implicit None
  Contains
    
    Subroutine img(stuu, dP)
      Use Param
      Implicit None
      Integer(4) :: dP
      Real(4) :: stuu(L,H)
      Character(len=20) :: Namefile
      
      Write(Namefile, '("./img/",i0,".pgm")') dP

      
      Open(10, file = Namefile, action='Write')
      Write(10, '(a)') "P2"
      Write(10, '(i0," ",i0)') Shape(stuu)
      Write(10, '(i0)') 1000
      If (Maxval(stuu) - Minval(stuu) >= 1D-5) Then
        Write(10, '(i0)') Int((stuu - Minval(stuu))*1000._8/(Maxval(stuu)-Minval(stuu)))
      Else
        stuu = 0
        Write(10, '(i0)') Int(stuu)
      End If
      Close(10)
    End Subroutine
    
    Subroutine Object(Obj)
      Use Param
      Implicit None
      Integer(1) :: Obj(L,H)
      Character(len=140) :: command
      Character(LEN=20) :: trash
      Integer(4) :: i
      
      Write(command, '("convert -background white -compress None objet.svg -resize ",&
        &i0,"x",i0,&
        &"^ -gravity center -extent ",&
        &i0,"x",i0,&
        &" -threshold 99% objet.pbm")') L, H, L, H

      Call execute_command_line(command, wait=.true.)

      Open(10, file = 'objet.pbm', status = 'old', action='read')
      Read(10, *) trash
      Read(10, *) trash
    
      Do i=1, H
        Read(10, *) Obj(:, i)
      End Do
    	
      Close(10)
 
    End Subroutine
    
    Subroutine InitF(F)
      Use Param
      Implicit None
      Real(4) :: F(Q,L,H) 
      Integer(1) :: i
      
      Do i=1, Q
          F(i,:,:) = w(i) 
      End Do
    End Subroutine
    
    Subroutine Stream(F)
      Use Param
      Implicit None
      Real(4) :: F(Q,L,H)
      Integer(1) :: i
      !$Omp Parallel Do
        Do i=2, Q
          If (e(i,1) /= 0) Then
            F(i,:,:) = Cshift(F(i,:,:), e(i,1), 1)
          End if
          If (e(i,2) /= 0) Then
            F(i,:,:) = Cshift(F(i,:,:), e(i,2), 2)
          End if
        End Do
      !$Omp End Parallel Do
    End Subroutine
    
    Subroutine Collision(F)
      Use Param
      Use Functions
      Implicit None
      Real(4) :: F(Q,L,H)
      F = F*mitau - itau*Feq(F)
    End Subroutine
    
End Module Subroutines

Module Var
  ! Initializes all the Variables
  Use Param
  Implicit None
  Integer(1) :: Obj(L, H)
  Real(4) :: F(Q,L,H)
End Module Var

Program LBM2D
  Use Param
  Use Functions
  Use Subroutines
  Use Var
  Implicit None
  Integer(4) :: i, j, dP = 0
  
  Call Object(Obj)
  Call InitF(F)
  Do i = 1, L
    Do j = 1, H
      F(:, i, j) = 10._8*exp(-((i-(L/2))**2+(j-(H/2))**2)/200._8) 
    End Do
  End Do
  
  Do dP = 0, P
    If (modulo(dP,10) == 0) Then
      Call img(uu(F)/(R(F)*R(F)), Int(dP/10))
    End If
    Call Stream(F)
    Call Collision(F)
  End Do
End Program LBM2D
