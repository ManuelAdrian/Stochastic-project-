clc
clear all

H_0 = 500;
A_1 = 700;
A_2 = 800;

mu_h = 1/(70*365);

z_1 = 0.04;
tz_1 = 0.1;
tz_2 = 0.11;
tpi_h = 0.0009;
pi_h = 0.0009;
pi_v = 0.03;
tpi_v = 0.49;

mu_1 = 1/(6*365);
mu_2 = 1/(2.5*365);
a_1 = 0.7714; 
a_2 = 0.7;
e_1 = 0.0067;
e_2 = 0.007;
K_1 = 35*(H_0+A_1);
K_2 = 150*A_2;

r_1 = a_1 - e_1;
r_2 = a_2 - e_2;

w_1 = 0.9;
w_2 = 0.9;
eta_1 = 0.1;


%%%%%%%% Parámetros adimensionalizados %%%%%%%%

tmu_h = mu_h/r_1;
tmu_1 = mu_1/r_1;
tmu_2 = mu_2/r_1;
te_1 = e_1/r_1;
te_2 = e_2/r_1;
delta_1 = a_1/r_1;
delta_2 = a_2/r_1;
r = r_2/r_1;

sigma_h = pi_v/r_1;
alpha_1 = tpi_h*K_1/(r_1*H_0);
tsigma_2 = (tpi_v/r_1)*(A_2/(A_2+(1-eta_1)*w_1*A_1));
alpha_2 = pi_h*K_2/(r_1*(1-eta_1)*w_1*A_1+r_1*A_2);
tsigma_1 = (tpi_v/r_1)*((1-eta_1)*w_1*A_1/(A_2+(1-eta_1)*w_1*A_1));
sigma_11 = (tpi_v*eta_1*w_1)/(r_1*(1-w_1)+eta_1*w_1*r_1);
sigma_12 = (tpi_v*(1-w_1)*w_2)/(r_1*(1-w_1)+eta_1*w_1*r_1);
talpha_11 = ((pi_h*eta_1*w_1)/(eta_1*w_1*A_1+(1-w_1)*A_1))*(K_1/r_1);
talpha_12 = ((pi_h*(1-w_1)*w_2)/(eta_1*w_1*A_1+(1-w_1)*A_1))*(K_1/r_1);
talpha_2 = (pi_h*(1-eta_1)*w_1/(A_2+(1-eta_1)*w_1*A_1))*(K_2/r_1);

%%%%%%%%%%%%%%%%%%%%%% R0s %%%%%%%%%%%%%%%%%%%%%%

%R_hv1 = sqrt(sigma_h*alpha_1/(delta_1*tmu_h));
%R_v1a1 = sqrt(sigma_1*talpha_1/(delta_1*tmu_1));
%R_a1v2 = sqrt(tsigma_1*talpha_2/(delta_2*tmu_1));
%R_v2a2 = sqrt(tsigma_2*alpha_2/(delta_2*tmu_2));

%R_1 = sqrt(R_hv1^2 + R_v1a1^2);
%R_v1a1v2 = sqrt(R_v1a1^2 + R_a1v2^2);
%R_2 = sqrt(R_a1v2^2 + R_v2a2^2);

%R_0 = sqrt(((R_1^2 + R_2^2)/2) + sqrt(((R_1^2 - R_2^2)/2)^2 + R_v1a1^2 * R_a1v2^2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = 100000;
N = 10^8;
dt = T/N;

j1 = 1;

x_0 = [0.1 0.25 0.1 0.45 0.25 0.1 0.5];
nx = size(x_0,2);

R = 4;
Dt = R*dt;
L = N/R; % L EM steps of size Dt = R*dt

beta1 = 0.01;%0.01;
beta2 = 0.05;%0.01;
beta3 = 0.05;%0.05;

for i=1:nx
    x_temp(i) = x_0(i);
end


for i2=1:100
    dW1 = sqrt(dt)*randn(1,N);
    dW2 = sqrt(dt)*randn(1,N);
    dW3 = sqrt(dt)*randn(1,N);
    
    for j = 1:L
        Wdisc1 = sum(dW1(R*(j-1)+1:R*j));
        Wdisc2 = sum(dW2(R*(j-1)+1:R*j));
        Wdisc3 = sum(dW3(R*(j-1)+1:R*j));
        
        a = [-(z_1*alpha_1*x_temp(3) + tmu_h) , -((tz_1*talpha_11 + tz_2*talpha_12)*x_temp(3) + tz_1*talpha_2*x_temp(6) + tmu_1) , -(z_1*sigma_h*x_temp(1) + (tz_1*sigma_11 + tz_2*sigma_12)*x_temp(2) + te_1 + x_temp(4)) , (1 - x_temp(4)) , -(tz_2*alpha_2*x_temp(6) + tmu_2) , -(tz_1*tsigma_1*x_temp(2) + tz_2*tsigma_2*x_temp(5) + te_2 + r*x_temp(7)) , r*(1 - x_temp(7))];
        b = [z_1*alpha_1*x_temp(3) , (tz_1*talpha_11 + tz_2*talpha_12)*x_temp(3) + tz_1*talpha_2*x_temp(6) , (z_1*sigma_h*x_temp(1) + (tz_1*sigma_11 + tz_2*sigma_12)*x_temp(2))*x_temp(4) , 0 , tz_2*alpha_2*x_temp(6) , (tz_1*tsigma_1*x_temp(2) + tz_2*tsigma_2*x_temp(5))*x_temp(7) , 0];
        
        for i=1:nx
            if a(i) == 0
                c1(i) = 1;
                c2(i) = 0;
            else
                c1(i) = 0;
                c2(i) = 1;
            end
        end
        
        A1 = zeros(nx);
        A2 = zeros(nx);
        
        for i=1:nx
            A1(i,i) = exp(Dt*a(i));
            A2(i,i) = ((A1(i,i) - 1)/a(i))*c2(i) + Dt*c1(i);
        end
        
        G = [beta1*alpha_1*x_temp(3)*(1 - x_temp(1)) 0 0 ; 0 beta2*(talpha_11*x_temp(3) + talpha_2*x_temp(6))*(1 - x_temp(2)) beta3*talpha_12*x_temp(3)*(1 - x_temp(2)) ; beta1*sigma_h*x_temp(1)*(x_temp(4) - x_temp(3)) beta2*sigma_11*x_temp(2)*(x_temp(4) - x_temp(3)) beta3*sigma_12*x_temp(2)*(x_temp(4) - x_temp(3)) ; 0 0 0 ; 0 0 beta3*alpha_2*x_temp(6)*(1 - x_temp(5)) ; 0 beta2*tsigma_1*x_temp(3)*(x_temp(7) - x_temp(6)) beta3*tsigma_2*x_temp(5)*(x_temp(7) - x_temp(6)) ; 0 0 0];
        F = A1*x_temp' + A2*b' + G*[Wdisc1 Wdisc2 Wdisc3]';
        
        for i=1:nx
            x_temp(i) = F(i);    
        end
        
        if mod(Dt*j,2500) == 0
            X1(i2,j1) = x_temp(1);      %(iteraciones , dias)
            X2(i2,j1) = x_temp(2);
            X3(i2,j1) = x_temp(3);
            X4(i2,j1) = x_temp(4);
            X5(i2,j1) = x_temp(5);
            X6(i2,j1) = x_temp(6);
            X7(i2,j1) = x_temp(7);
        
            j1 = j1+1;
        end
    end
    
    j1 = 1;
%h2a = plot([0:Dt:T],[x1_0,X1], 'LineWidth',.5); h2a.Color(4)=0.3;
%hold on
end
