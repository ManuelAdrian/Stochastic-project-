PROGRAM SDE_example

IMPLICIT NONE

integer, parameter :: dp = selected_real_kind(6, 37)
integer, parameter :: dp1 = selected_real_kind(15, 307)

integer,parameter :: N = 2**8, R1 = 2.0d0, L = N/R1, N1 = 10000.0d0

real(kind=dp) :: r,sigma,S0,dt,Dt1,Xtemp,Winc,Esp1,Var0	

real(kind=dp1) :: T

real :: pi, temp, mean, sd, dW(N), Xe(N1,L), array(N1,N), Esp(L), Var(L)
integer :: i,j,i1

100     FORMAT(1X,8000E17.6E3) !Printing format
!110     FORMAT(A,F5.3)
CHARACTER(15) filename, filename1, filename2

!$$$$$$$$$$$$ Distribucion Normal $$$$$$$$$$$$!
mean = 0.0d0
sd = 1.0d0
pi = 4.0*ATAN(1.0)

DO i1 = 1, N1

call random_number(array(i1,:))

DO i = 1, N-1, 2
temp = sd * SQRT(-2.0*LOG(array(i1,i))) * COS(2*pi*array(i1,i+1)) + mean
array(i1,i+1) = sd * SQRT(-2.0*LOG(array(i1,i))) * SIN(2*pi*array(i1,i+1)) + mean
array(i1,i) = temp
END DO

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!

dW = SQRT(dt)*array(i1,:)

r = 0.05d0
sigma = 0.2d0
S0 = 1.0d0
T = 1.0d0
dt = T/N
Dt1 = R1*dt
Xtemp = S0

DO j = 1, L
Winc = sum(dW(R1*(j-1)+1: R1*j))
Xtemp = Xtemp + Dt1*r*Xtemp + sigma*Xtemp*Winc
Xe(i1,j) = Xtemp
END DO

END DO

filename = 'solu_aprox.dat'
OPEN(UNIT=45,FILE=filename,position="APPEND")

DO i1 = 1, N1
write(45,100) (Xe(i1,j),j=1,L)
END DO

CLOSE(45)

DO j = 1, L
Esp1 = 0.0d0  
DO i1 = 1, N1
Esp1 = Esp1 + Xe(i1,j)
END DO
Esp(j) = Esp1/N1
END DO

DO j = 1, L
  Var0 = 0.0d0
  DO i1 = 1, N1
  Var0 = Var0 + (Xe(i1,j)-Esp(j))**2
  END DO
  Var(j) = Var0/N1
END DO

filename1 = 'esperanza.dat'
filename2 = 'varianza.dat'
OPEN(UNIT=46,FILE=filename1,position="APPEND")
OPEN(UNIT=47,FILE=filename2,position="APPEND")

DO j =1, L
write(46,100) Esp(j)
write(47,100) Var(j)
END DO

CLOSE(46)
CLOSE(47)

END PROGRAM SDE_example