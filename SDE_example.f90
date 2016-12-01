PROGRAM SDE_example

IMPLICIT NONE

integer, parameter :: dp = selected_real_kind(6, 37)
integer, parameter :: dp1 = selected_real_kind(15, 307)

real(kind=dp) :: r,sigma,S0,dt,Dt1,Xtemp,Winc,Wtemp,dttemp,Var,Esp,Var0	!,W

real(kind=dp1) :: T

integer,parameter :: N = 2**8, R1 = 2.0d0, L = N/R1
real :: array(N), pi, temp, mean, sd, dW(N), Strue(N), Xe(L)
integer :: i,j
integer :: values(1:8), k
integer, dimension(:), allocatable :: seed

100     FORMAT(1X,8000E17.6E3) !Printing format
110     FORMAT(A,F5.3)
CHARACTER(15) filename, filename1, filename2

r = 0.05d0
sigma = 0.2d0
S0 = 1.0d0
T = 1.0d0
dt = T/N
Dt1 = R1*dt
Xtemp = S0
dttemp = dt

!$$$$$$$$$$$$ Distribucion Normal $$$$$$$$$$$$!
mean = 0.0d0
sd = 1.0d0
pi = 4.0*ATAN(1.0)

call date_and_time(values=values)

call random_seed(size=k)
allocate(seed(1:k))
seed(:) = values(8)
call random_seed(put=seed)
call random_number(array)

DO i = 1, N-1, 2
temp = sd * SQRT(-2.0*LOG(array(i))) * COS(2*pi*array(i+1)) + mean
array(i+1) = sd * SQRT(-2.0*LOG(array(i))) * SIN(2*pi*array(i+1)) + mean
array(i) = temp
END DO

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!

dW = SQRT(dt)*array
Wtemp = 0.0d0

filename = 'solu_exact.dat'
OPEN(UNIT=45,FILE=filename,position="APPEND")
write(45,100) S0

DO j = 1, N
Wtemp = Wtemp + dW(j)
dttemp = dttemp + dt

Strue(j) = S0*exp((r-0.5d0*(sigma**2))*dttemp+sigma*Wtemp)
write(45,100) Strue(j)
END DO

CLOSE(45)

filename1 = 'solu_aprox.dat'
OPEN(UNIT=46,FILE=filename1,position="APPEND")
write(46,100) Xtemp

DO j = 1, L
Winc = sum(dW(R1*(j-1)+1: R1*j))
Xtemp = Xtemp + Dt1*r*Xtemp + sigma*Xtemp*Winc
Xe(j) = Xtemp
write(46,100) Xtemp
END DO

CLOSE(46)

Esp = (sum(Xe) + S0)/(L+1)
Var0 = (S0-Esp)**2

DO j = 1,L
Var0 = Var0 + (Xe(j)-Esp)**2 
END DO

Var = Var0/(L+1)

filename2 = 'estimad.dat'
OPEN(UNIT=47,FILE=filename2,position="APPEND")
write(47,110) 'valor esperado = ',Esp
write(47,110) 'varianza = ', Var
CLOSE(47)

END PROGRAM SDE_example