!  Author : Pacidus
! Started at: 31.10.2020 00:17:57

Module Param
  ! Initializes all the constants
  Implicit None
  Private
  Real(8), Parameter :: lx         = 4._8                           ! m
  Real(8), Parameter :: ly         = 2._8                           ! m
  Real(8), Parameter :: dt         = 5D-5                           ! s
  Real(8), Public, Parameter :: cs = 340._8                         ! m·s⁻¹
  Real(8), Parameter :: c          = cs*Dsqrt(3._8)                 ! m·s⁻¹
  Real(8), Parameter :: dl         = cs*dt                          ! m
  Real(8), Parameter :: t          = 1D4                            ! s
  Real(8), Parameter :: rho        = 1.177_8                        ! kg·m⁻³
  Real(8), Parameter :: nu         = 1.57D-5                        ! m²·s⁻¹
  Real(8), Parameter :: tau        = (nu/(dt*cs**2)) + 5D-1         ! ø
  Real(8), Public, Parameter :: itau       = 1._8/tau               ! ø
  Real(8), Public, Parameter :: mitau      = 1._8 - itau            ! ø
  Integer(1), Public, Parameter :: D       = 2                      ! D2
  Integer(1), Public, Parameter :: Q       = 9                      ! Q9
  Integer(4), Public, Parameter :: H       = int(ly/dl)             ! ø
  Integer(4), Public, Parameter :: L       = int(lx/dl)             ! ø
  Integer(4), Public, Parameter :: P       = int(t/dt)              ! ø
  Integer(1), Public, Parameter :: ed(Q,D) = transpose(reshape(&    ! ø
    (/0, 0,&
      1, 0,&
      0, 1,&
     -1, 0,&
      0,-1,&
      1, 1,&
     -1, 1,&
     -1,-1,&
      1,-1/), (/D,Q/)))
  Real(8), Public, Parameter :: e(Q,D) = c*ed
  Real(8), Public, Parameter :: w(Q) = &                            ! ø
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
  Private
  
  Public :: Feq, R, U
  
  Contains
  
    Function R(F)
      Use Param, Only: Q, L, H, w
      Implicit None
      Real(8) :: R(L, H)
      Real(8) :: F(Q, L, H)
      Integer(1) :: i
      R = 0
      Do i = 1, Q
        R = R + F(i,:,:)
      End Do
    End Function

    Function U(F, Rv)
      Use Param, Only: D, Q, L, H, e
      Implicit None
      Real(8) :: U(D, L, H)
      Real(8) :: Rv(L, H)
      Real(8) :: F(Q, L, H)
      Integer(1) :: i

      U = 0
      
      !Omp Parallel Do
        Do i = 2, Q
          If (e(i, 1) /= 0) Then
            U(1, :, :) = U(1, :, :) + F(i, :, :)*e(i, 1)
          End If
          If (e(i, 2) /= 0) Then
            U(2, :, :) = U(2, :, :) + F(i, :, :)*e(i, 2)
          End If
        End Do
      !Omp End Parallel Do

      U(1, :, :) = U(1, :, :)/Rv
      U(2, :, :) = U(2, :, :)/Rv
    End Function
     
     Function Feq(F)
      Use Param, Only: D, Q, L, H, w, cs, e
      Implicit None
      Real(8) :: Feq(Q, L, H)
      Real(8) :: F(Q, L, H)
      Real(8) :: stu(D, L, H)
      Real(8) :: stR(L, H)
      Real(8) :: steu(L, H)
      Real(8) :: stuu(L, H)
      Integer(1) :: i
      
      stR = R(F)
      stu = U(F, stR)
      stuu = stu(1, :, :)*stu(1, :, :) + stu(2, :, :)*stu(2, :, :)
      
      Feq(1,:,:) = w(1)*stR*(1 - stuu/(2*cs**2))
      !$Omp Parallel Do
        Do i = 2, Q
          steu = e(i, 1)*stu(1, :, :) + e(i, 2)*stu(2, :, :)
          Feq(i,:,:) = w(i)*stR*(cs**2 + steu + ((steu*steu)/(cs**2) - stuu)/2)/(cs**2)
        End Do
      !$Omp End Parallel Do
     End Function
End Module Functions

Module Subroutines
  ! Initializes all the Subroutines
  Implicit None
  Private
  
  Public :: Initf
  Public :: img
  Public :: Stream
  Public :: Collision
  Public :: Object
  
  Contains
    
    Subroutine img(stuu, dP)
      Use Param, Only: L, H
      Implicit None
      Integer(4) :: dP
      Real(8) :: stuu(L,H)
      Character(len=20) :: Namefile
      
      Write(Namefile, '("./img/",i0,".pgm")') dP

      
      Open(10, file = Namefile, action='Write')
      Write(10, '(a)') "P2"
      Write(10, '(i0," ",i0)') Shape(stuu)
      Write(10, '(i0)') 1000
      Print *, Minval(stuu), Maxval(stuu)
      If (Maxval(stuu) - Minval(stuu) >= 1D-5) Then
        Write(10, '(i0)') Int((stuu - Minval(stuu))*1000._8/(Maxval(stuu)-Minval(stuu)))
      Else
        stuu = 0
        Write(10, '(i0)') Int(stuu)
      End If
      Close(10)
    End Subroutine
    
    Subroutine Object(Obj)
      Use Param, Only: L, H
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
      Use Param, Only: Q, L, H, w
      Implicit None
      Real(8) :: F(Q,L,H) 
      Integer(1) :: i
      
      Do i=1, Q
          F(i,:,:) = w(i) 
      End Do
    End Subroutine
    
    Subroutine Stream(F)
      Use Param, Only: Q, L, H, ed
      Implicit None
      Real(8) :: F(Q,L,H)
      Integer(1) :: i
      !$Omp Parallel Do
        Do i=2, Q
          If (ed(i,1) /= 0) Then
            F(i,:,:) = Cshift(F(i,:,:), ed(i,1), 1)
          End if
          If (ed(i,2) /= 0) Then
            F(i,:,:) = Cshift(F(i,:,:), ed(i,2), 2)
          End if
        End Do
      !$Omp End Parallel Do
    End Subroutine
    
    Subroutine Collision(F)
      Use Param, Only: D, Q, L, H, itau, Mitau
      Use Functions, Only: Feq, U, R
      Implicit None
      Real(8) :: F(Q,L,H)
      Real(8) :: eqF(Q,L,H)
      Real(8) :: S(D, L,H)
      Integer(1) :: i
      eqF = Feq(F)
      S = 0.00001*U(F, R(F))
      Do i = 1, Q
        F(i, :, :) = (F(i, :, :)*mitau) + (itau*eqF(i, :, :)) - (S(1, :, :)**2+S(2, :, :)**2)
      End Do
    End Subroutine
    
End Module Subroutines

Module Var
  ! Initializes all the Variables
  Use Param, Only: Q, L, H
  Implicit None
  Private
  Integer(1), Public :: Obj(L, H)
  Real(8), Public :: F(Q,L,H)
End Module Var

Program LBM2D
  Use Param, Only: P, L, H, itau, w, Q
  Use Functions, Only: R
  Use Subroutines
  Use Var
  Implicit None
  Integer(4) :: i, j, k, dP = 0
  
  Print *, L, H, itau
  
  Call Object(Obj)
  Call InitF(F)
  Do i = 1, L
    Do j = 1, H
      Do k = 1, Q
        F(k, i, j) = w(k)*(10 + 1*exp(-((i-(L/2))**2+(j-(H/2))**2)/1000._8)) 
      End Do
    End Do
  End Do
  
  Do dP = 0, P
    If (modulo(dP,10) == 0) Then
      Call img(R(F), Int(dP/10))
    End If
    Call Stream(F)
    Call Collision(F)
  End Do
End Program LBM2D
