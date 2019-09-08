clc;clear;clear all; close
dtr = pi/180;

%%
r1 = 0.33;
r2 = 0.055;
r3 = 0.283;
r3p = 0.17;
r4 = 0.322;
r5p = 0;
r6 = 0.165;
r7 = 0.061;
rw = 0;
t =  0.057;
tp = 0.013;

gr = 2.27;
sr = 1.71;
IJ = 0.9;%kg-m*m
IF = 0;%kg-m*m

m2 = 53;
m3 = 7.1;
m4 = 10.59;
m5 = 47.93;
m6 = 50.03;
m7 = 60.07;

L2 = 0.0247;
L3 = 0.1287;
L4 = 0.1467;
L5 = 0;
L6 = 0.0749;
L7 = 0.0277;

phi2 = 0;
phi3 = 13*dtr;

phi4 = 0;
phi5 = 0;
phi6 = 0;
phi7 = 0;
pressure_angle = 20*dtr;
% 
I2 = 0.651;
I3 = 0.047;
I4 = 0.887;
I6 = 0.313;
I7 = 0.955;
% I2 = 0;
% I3 = 0;
% I4 = 0;
% I6 = 0;
% I7 = 0;

% %
spring_initial = 0.27;  %m
stiffness = 37000;      %N/m
damping = 11000; %Ns/m
% spring_initial = 0.27;  %m
% stiffness = 0;      %N/m
% damping = 0; %Ns/m

g = 9.81;
%%
theta_1 = (360-107.47)*dtr;
theta_2i = 0;
theta_3p = 27*dtr;
theta_7i = 287.45*dtr;
theta_5 = 270*dtr;

theta_2 = zeros(360,1);
theta_7 = zeros(360,1);


theta_6 = zeros(360,1);
theta_3 = zeros(360,1);
r5 = zeros(360,1);
theta_4 = zeros(360,1);
theta_d = zeros(360,1);

h3 = zeros(360,1);
hr5 = zeros(360,1);
h4 = zeros(360,1);
h6 = zeros(360,1);

h3p = zeros(360,1);
hr5p = zeros(360,1);
h4p = zeros(360,1);
h6p = zeros(360,1);


%% Displacement
interval = 1*dtr;
A = 2*(r1*r6*cos(theta_1)-r2*r6*cos(theta_2)+r6*r7*cos(theta_7));
B = 2*(r1*r6*sin(theta_1)-r2*r6*sin(theta_2)+r6*r7*sin(theta_7));
C = r3^2-r1^2-r2^2-r6^2-r7^2+2*r1*r2*cos(theta_1-theta_2)-2*r1*r7*cos(theta_1-theta_7)+2*r2*r7*cos(theta_2-theta_7);

for i = 1:360
    theta_2(i) = theta_2i+(i-1)*interval;
    theta_7(i) = theta_7i-(i-1)*interval;
    A = 2*(r1*r6*cos(theta_1)-r2*r6*cos(theta_2(i))+r6*r7*cos(theta_7(i)));
    B = 2*(r1*r6*sin(theta_1)-r2*r6*sin(theta_2(i))+r6*r7*sin(theta_7(i)));
    C = r3^2-r1^2-r2^2-r6^2-r7^2+2*r1*r2*cos(theta_1-theta_2(i))-2*r1*r7*cos(theta_1-theta_7(i))+2*r2*r7*cos(theta_2(i)-theta_7(i));
    
    
    
    if  (2*atan2(B+((B.*B+A.*A-C.*C).^0.5),A+C) + 360*dtr)/dtr < 360
        theta_6(i) = 2*atan2(B+((B.*B+A.*A-C.*C).^0.5),A+C ) + 360*dtr;
    else
        theta_6(i) = 2*atan2(B+((B.*B+A.*A-C.*C).^0.5),A+C) -360*dtr;
    end
    
    theta_3(i) = -(acos((r1*cos(theta_1)-r2*cos(theta_2(i))+r6*cos(theta_6(i))+r7*cos(theta_7(i)))/r3))+360*dtr;
    
    D = 2*(-r1*cos(theta_1-theta_5)+r5p*cos(theta_5)-r6*cos(theta_5-theta_6(i))-r7*cos(theta_5-theta_7(i)));
    E = r1^2-r4^2+r5p^2+r6^2+r7^2-2*r1*r5p*cos(theta_1)+2*r1*r6*cos(theta_1-theta_6(i))+2*r1*r7*cos(theta_1-theta_7(i))...
        +2*r6*r7*cos(theta_6(i)-theta_7(i))-2*r5p*r6*cos(theta_6(i))-2*r5p*r7*cos(theta_7(i));
    
    r5(i) =(-D+sqrt(D^2-4*E))/2;
    
    theta_4(i) = acos((-r1*cos(theta_1)+r5(i)*cos(theta_5)+r5p-r6*cos(theta_6(i))-r7*cos(theta_7(i)))/r4);
    
    if 0 < theta_4(i) && theta_4(i) < 90*dtr
        theta_4(i) = -theta_4(i)+360*dtr;
    else
        theta_4(i) = 360*dtr-theta_4(i);
    end
end

% figure_chienli2(300,260,0.23,0.25,0.7,0.6)
%% Angular velocity
omega =[43.9396];
 for k = 1:1

omega_2 = [0 0 omega(k)];
omega_7 = [0 0 -omega(k)];
omega_6 = zeros(360,3);
omega_3 = zeros(360,3);
omega_4 = zeros(360,3);
r5v =  zeros(360,3);
X = zeros(360,2);

for i = 1:360
    AA = [r3*sin(theta_3(i)) -r6*sin(theta_6(i));-r3*cos(theta_3(i)) r6*cos(theta_6(i))];
    BB = [-r2*sin(theta_2(i))-r7*sin(theta_7(i));r2*cos(theta_2(i))+r7*cos(theta_7(i))];
    X(i,:) = AA\BB;
    h3(i) = X(i,1);
    h6(i) = X(i,2);
    omega_3(i,3) = h3(i)*omega_2(1,3);
    omega_6(i,3) = h6(i)*omega_2(1,3);
    
    
    AA = [-r4*sin(theta_4(i)) 0;r4*cos(theta_4(i)) 1];
    BB = [r6*h6(i)*sin(theta_6(i))-r7*sin(theta_7(i));-r6*h6(i)*cos(theta_6(i))+r7*cos(theta_7(i))];
    X(i,:) = AA\BB;
    h4(i) = X(i,1);
    hr5(i) = X(i,2);
    omega_4(i,3) = h4(i)*omega_2(1,3);
    r5v(i,2) = -hr5(i)*omega_2(1,3);
end
%% Angular acceleration

alpha_6 = zeros(360,3);
alpha_3 = zeros(360,3);
alpha_4 = zeros(360,3);
r5a = zeros(360,3);

for i = 1:360
    AA = [r3*sin(theta_3(i)) -r6*sin(theta_6(i));-r3*cos(theta_3(i)) r6*cos(theta_6(i))];
    BB = [h6(i)^2*r6*cos(theta_6(i))+r7*cos(theta_7(i))-r2*cos(theta_2(i))-h3(i)^2*r3*cos(theta_3(i));
        h6(i)^2*r6*sin(theta_6(i))+r7*sin(theta_7(i))-r2*sin(theta_2(i))-h3(i)^2*r3*sin(theta_3(i))];
    X(i,:) = AA\BB;
    h3p(i) = X(i,1);
    h6p(i) = X(i,2);
    alpha_3(i,3) = h3p(i)*omega_2(1,3)^2;
    alpha_6(i,3) = h6p(i)*omega_2(1,3)^2;
    
    AA = [-r4*sin(theta_4(i)) 0;r4*cos(theta_4(i)) 1];
    BB = [r4*h4(i)^2*cos(theta_4(i))+r6*h6p(i)*sin(theta_6(i))+r6*h6(i)^2*cos(theta_6(i))+r7*cos(theta_7(i));
        r4*h4(i)^2*sin(theta_4(i))-r6*h6p(i)*cos(theta_6(i))+r6*h6(i)^2*sin(theta_6(i))+r7*sin(theta_7(i))];
    X(i,:) = AA\BB;
    h4p(i) = X(i,1);
    hr5p(i) = X(i,2);
    
    alpha_4(i,3) = h4p(i)*omega_2(1,3)^2;
    r5a(i,2) = -hr5p(i)*omega_2(1,3)^2;
end



%% Linkage displacement
PG = [0 -0.75 0];
PH = [0.22 0.44 0];

PO2 = [0 0 0];
PO7 = [PO2(1,1) + r1*cos(theta_1) ,PO2(1,2) + r1*sin(theta_1) , 0];

PA = zeros(360,3);
PB = zeros(360,3);
PD = zeros(360,3);
PR5 = zeros(360,3);
PC = zeros(360,3);
PS = zeros(360,3);    %the joint between damper and r3`


for i = 1:360
    for j = 1:2
        PA(i,1) = PO2(1,1) + r2*cos(theta_2(i));
        PA(i,2) = PO2(1,2) + r2*sin(theta_2(i));
        
        PB(i,1) = PA(i,1) + r3*cos(theta_3(i));
        PB(i,2) = PA(i,2) + r3*sin(theta_3(i));
        
        PD(i,1) = PO7(1,1) + r7*cos(theta_7(i));
        PD(i,2) = PO7(1,2) + r7*sin(theta_7(i));
        
        PR5(i,1) = PO2(1,1);
        PR5(i,2) = PO2(1,2) - r5(i);
        
        PC(i,1) = PB(i,1) + r4*cos(theta_4(i));
        PC(i,2) = PB(i,2) + r4*sin(theta_4(i));
        
        PS(i,1) = PA(i,1) + r3p*cos(theta_3(i) + theta_3p);
        PS(i,2) = PA(i,2) + r3p*sin(theta_3(i) + theta_3p);
        
        %  XC(k) = XA(k) + r3*cos(theta_3(k))+r4*cos(theta_4(k));
    end
end

%% Displacement vector

vector_O2_A = zeros(360,3);
vector_A_B = zeros(360,3);
vector_B_C = zeros(360,3);
vector_O7_D = zeros(360,3);
vector_D_B = zeros(360,3);

for i = 1:360
    for j = 1:2
        vector_O2_A(i,j) = PA(i,j)-PO2(1,j);
        vector_A_B(i,j) = PB(i,j)-PA(i,j);
        vector_B_C(i,j) = PC(i,j)-PB(i,j);
        vector_O7_D(i,j) = PD(i,j)-PO7(1,j);
        vector_D_B(i,j) = PB(i,j)-PD(i,j);
        
    end
end

%% Linkage mass-center displacement

PG_r2 = zeros(360,3);
PG_r3 = zeros(360,3);
PG_r4 = zeros(360,3);
PG_r5 = zeros(360,3);
PG_r6 = zeros(360,3);
PG_r7 = zeros(360,3);


vector_O2_Gr2 = zeros(360,3);
vector_A_Gr3 = zeros(360,3);
vector_B_Gr4 = zeros(360,3);
vector_O7_Gr7 = zeros(360,3);
vector_D_Gr6 = zeros(360,3);
vector_B_C = zeros(360,3);

for i = 1 :360
    PG_r2(i,1) = PO2(1,1) + L2*cos(theta_2(i)+phi2);
    PG_r2(i,2) = PO2(1,2) + L2*sin(theta_2(i)+phi2);
    PG_r3(i,1) = PA(i,1) + L3*cos(theta_3(i)+phi3);
    PG_r3(i,2) = PA(i,2) + L3*sin(theta_3(i)+phi3);
    PG_r4(i,1) = PB(i,1) + L4*cos(theta_4(i)+phi4);
    PG_r4(i,2) = PB(i,2) + L4*sin(theta_4(i)+phi4);
    
    PG_r5(i,1) = PB(i,1) + r4*cos(theta_4(i));
    PG_r5(i,2) = PB(i,2) + r4*sin(theta_4(i));
    
    PG_r6(i,1) = PD(i,1) + L6*cos(theta_6(i)+phi6);
    PG_r6(i,2) = PD(i,2) + L6*sin(theta_6(i)+phi6);
    PG_r7(i,1) = PO7(1,1) + L7*cos(theta_7(i)+phi7);
    PG_r7(i,2) = PO7(1,2) + L7*sin(theta_7(i)+phi7);
    for j = 1:2
        vector_O2_Gr2(i,j) = PG_r2(i,j)-PO2(1,j);
        vector_A_Gr3(i,j) = PG_r3(i,j)-PA(i,j);
        vector_B_Gr4(i,j) = PG_r4(i,j)-PB(i,j);
        vector_O7_Gr7(i,j) = PG_r7(i,j)-PO7(1,j);
        vector_D_Gr6(i,j) = PG_r6(i,j)-PD(i,j);
        
    end
end

%% Linkage velocity

VO2 = [0 0 0];
VO7 = [0 0 0];
VA = zeros(360,3);
VB = zeros(360,3);
VD = zeros(360,3);


for i = 1:360
    
    VA(i,:) = cross(omega_2(1,:),vector_O2_A(i,:));
    VB(i,:) = cross(omega_3(i,:),vector_A_B(i,:)) + VA(i,:);
    VD(i,:) = cross(omega_7(1,:),vector_O7_D(i,:));
    
end

%%  Linkage mass velocity
Vr2 = zeros(360,3);
Vr3 = zeros(360,3);
Vr4 = zeros(360,3);
Vr6 = zeros(360,3);
Vr7 = zeros(360,3);


for i = 1:360
    Vr2(i,:) = cross(omega_2(1,:),vector_O2_Gr2(i,:));
    Vr3(i,:) = cross(omega_3(i,:),vector_A_Gr3(i,:))+VA(i,:);
    Vr4(i,:) = cross(omega_4(i,:),vector_B_Gr4(i,:))+VB(i,:);
    Vr6(i,:) = cross(omega_6(i,:),vector_D_Gr6(i,:))+VD(i,:);
    Vr7(i,:) = cross(omega_7(1,:),vector_O7_Gr7(i,:));
end

%%  Linkage acceleration

AA = zeros(360,3);
AB = zeros(360,3);
AC = zeros(360,3);
AD = zeros(360,3);

for i = 1:360
    AA(i,:) = cross(omega_2(1,:),VA(i,:));
    AB(i,:) = AA(i,:) + cross(omega_3(i,:),cross(omega_3(i,:),vector_A_B(i,:))) + cross(alpha_3(i,:),vector_A_B(i,:));
    AC(i,:) = r5a(i,:);
    AD(i,:) = cross(omega_7(1,:),VD(i,:));
end

%% Linkage mass acceleration
Ar2 = zeros(360,3);
Ar3 = zeros(360,3);
Ar4 = zeros(360,3);
Ar6 = zeros(360,3);
Ar7 = zeros(360,3);
for i = 1:360
    Ar2(i,:) = cross(omega_2(1,:),Vr2(i,:));
    Ar3(i,:) = AA(i,:) + cross(omega_3(i,:),cross(omega_3(i,:),vector_A_Gr3(i,:))) + cross(alpha_3(i,:),vector_A_Gr3(i,:));
    
    Ar4(i,:) = AB(i,:) + cross(omega_4(i,:),cross(omega_4(i,:),vector_B_Gr4(i,:))) + cross(alpha_4(i,:),vector_B_Gr4(i,:));
    
    Ar6(i,:) = AD(i,:) + cross(omega_6(i,:),cross(omega_6(i,:),vector_D_Gr6(i,:))) + cross(alpha_6(i,:),vector_D_Gr6(i,:));
    
    Ar7(i,:) = cross(omega_7(1,:),Vr7(i,:));
    
end



%% spring_length and damper velocity
spring_length = zeros(360,1);
rd =  zeros(360,1);
rdv = zeros(360,3);

for i = 1:360
    rd(i) = sqrt((PS(i,1)-PH(1,1))^2+(PS(i,2)-PH(1,2))^2);
    spring_length(i) = (PO2(1,2)-PG(1,2)) -r5(i);
    
    
    theta_d(i) = atan2(PH(1,2)-PS(i,2),PH(1,1)-PS(i,1));
    
    AA = [cos(theta_d(i)) -rd(i)*sin(theta_d(i));sin(theta_d(i)) rd(i)*cos(theta_d(i))];
    BB = [omega_2(1,3)*r2*sin(theta_2(i))+r3p*omega_3(i,3)*sin(theta_3(i) + theta_3p);-omega_2(1,3)*r2*cos(theta_2(i))-r3p*omega_3(i,3)*cos(theta_3(i) + theta_3p)];
    
    X(i,:) = AA\BB;
    
    rdv(i,2) = X(i,1);
end

dead_center = spring_length(1);
for i = 1:360
    if dead_center > spring_length(i)
        dead_center = spring_length(i);
        
    end
end




%% 18*18 matrix



for i = 1:360;
    O2x(i) = PO2(1,1)-PG_r2(i,1);       %點 桿 座標
    O2y(i) = PO2(1,2)-PG_r2(i,2);
    A2x(i) = PA(i,1)-PG_r2(i,1);
    A2y(i) = PA(i,2)-PG_r2(i,2);
    G2x(i) = PO2(1,1)+ (r1/2)*cos(theta_1) - PG_r2(i,1);
    G2y(i) = PO2(1,2)+ (r1/2)*sin(theta_1) - PG_r2(i,2);
    
    A3x(i) = PA(i,1)-PG_r3(i,1);
    A3y(i) = PA(i,2)-PG_r3(i,2);
    S3x(i) = PS(i,1)-PG_r3(i,1);
    S3y(i) = PS(i,2)-PG_r3(i,2);
    B3x(i) = PB(i,1)-PG_r3(i,1);
    B3y(i) = PB(i,2)-PG_r3(i,2);
    
    B4x(i) = PB(i,1)-PG_r4(i,1);
    B4y(i) = PB(i,2)-PG_r4(i,2);
    C4x(i) = PC(i,1)-PG_r4(i,1);
    C4y(i) = PC(i,2)-PG_r4(i,2);
    
    B6x(i) = PB(i,1)-PG_r6(i,1);
    B6y(i) = PB(i,2)-PG_r6(i,2);
    D6x(i) = PD(i,1)-PG_r6(i,1);
    D6y(i) = PD(i,2)-PG_r6(i,2);
    
    O7x(i) = PO7(1,1)-PG_r7(i,1);
    O7y(i) = PO7(1,2)-PG_r7(i,2);
    D7x(i) = PD(i,1)-PG_r7(i,1);
    D7y(i) = PD(i,2)-PG_r7(i,2);
    G7x(i) = PO2(1,1) + (r1/2)*cos(theta_1) - PG_r7(i,1);
    G7y(i) = PO2(1,2) + (r1/2)*sin(theta_1) - PG_r7(i,2);
    
    
    Fs(i) = stiffness*(spring_length(i)-spring_initial);     %tension force,force down
    Fc(i) = -rdv(i,2)*damping;
     
    
%     Fg2x = cos(theta_1-90*dtr-pressure_angle);
%     Fg2y = sin(theta_1-90*dtr-pressure_angle);
%     Fg7x = -cos(theta_1-90*dtr-pressure_angle);
%     Fg7y = -sin(theta_1-90*dtr-pressure_angle);
    
    Fg2x = cos(theta_1+90*dtr-pressure_angle);
    Fg2y = sin(theta_1+90*dtr-pressure_angle);
    Fg7x = -cos(theta_1+90*dtr-pressure_angle);
    Fg7y = -sin(theta_1+90*dtr-pressure_angle);
%     
%     g1(i) = Fg2x*G2y(i)+Fg2y*G2x(i);
%     g2(i) = -Fg7x*G7y(i)-Fg7y*G7x(i);

    g1(i) = Fg2x*(r1/2*sin(theta_1+pi)+L2*sin(theta_2(i)))-Fg2y*(r1/2*cos(theta_1+pi)+L2*cos(theta_2(i)));   % WHY
    g2(i) = -Fg7x*(r1/2*sin(theta_1+pi)-L7*sin(theta_7(i)))-Fg2y*(r1/2*cos(theta_1+pi)-L7*cos(theta_7(i)));  % WHY

  
    if spring_length(i)-dead_center <= 0.005 && r5v(i,2) < 0
        Fw(i) = 89000;
       
    else
        Fw(i) = 0;
    end
    
%%force    1       2       3       4       5       6       7       8       9       10      11      12      13      14      15      16      17      18          Force
force = [  -1      0       1       0       0       Fg2x    0       0       0       0       0       0       0       0       0       0       0       0     ;     %F12x     1   Link2x
           0       1       0       -1      0       Fg2y    0       0       0       0       0       0       0       0       0       0       0       0     ;     %F12y     2   Link2y
           O2y(i)  O2x(i)  -A2y(i) -A2x(i) 1       g1(i)   0       0       0       0       0       0       0       0       0       0       0       0     ;     %F32x     3   Link2z
           0       0       -1      0       0       0       1       0       0       0       0       0       0       0       0       0       0       0     ;     %F32y     4   Link3x
           0       0       0       1       0       0       0       1       0       0       0       0       0       0       0       0       0       0     ;     %T        5   Link3y
           0       0       A3y(i)  A3x(i)  0       0       -B3y(i) B3x(i)  0       0       0       0       0       0       0       0       0       0     ;     %Fg       6   Link3z
           0       0       0       0       0       0       -1      0       1       0       0       0       1       0       0       0       0       0     ;     %F43x     7   Link4x
           0       0       0       0       0       0       0       -1      0       1       0       0       0       -1      0       0       0       0     ;     %F43y     8   Link4y
           0       0       0       0       0       0       B4y(i)  -B4x(i) -C4y(i) C4x(i)  0       0       -B4y(i) -B4x(i) 0       0       0       0     ;     %F54x     9   Link4z
           0       0       0       0       0       0       0       0       -1      0       1       1       0       0       0       0       0       0     ;     %F54y     10  Link5x
           0       0       0       0       0       0       0       0       0       -1      0       0       0       0       0       0       0       0     ;     %N15      11  Link5y
           0       0       0       0       0       0       0       0       0       0       t       -tp     0       0       0       0       0       0     ;     %N15p     12  Link5z
           0       0       0       0       0       0       0       0       0       0       0       0       -1      0       -1      0       0       0     ;     %F64x     13  Link6x
           0       0       0       0       0       0       0       0       0       0       0       0       0       1       0       1       0       0     ;     %F64y     14  Link6y
           0       0       0       0       0       0       0       0       0       0       0       0       B6y(i)  B6x(i)  D6y(i)  D6x(i)  0       0     ;     %F67x     15  Link6z
           0       0       0       0       0       Fg7x    0       0       0       0       0       0       0       0       1       0       -1      0     ;     %F67y     16  Link7x
           0       0       0       0       0       Fg7y    0       0       0       0       0       0       0       0       0       -1      0       1     ;     %F17x     17  Link7y
           0       0       0       0       0       g2(i)   0       0       0       0       0       0       0       0       -D7y(i) -D7x(i) O7y(i)  O7x(i)];    %F17y     18  Link7z
           


    maG2x(i) = Ar2(i,1)*m2;                                                                                     
    maG2y(i) = Ar2(i,2)*m2+m2*g;                                                                                
    moment2(i) = 0;                                                                                             
    maG3x(i) = Ar3(i,1)*m3-Fc(i)*cos(theta_d(i)+180*dtr);                                                       
    maG3y(i) = Ar3(i,2)*m3-Fc(i)*sin(theta_d(i)+180*dtr)+m3*g;                                                  
    moment3(i) = I3*alpha_3(i,3)+Fc(i)*cos(theta_d(i)+180*dtr)*S3y(i)-Fc(i)*sin(theta_d(i)+180*dtr)*S3x(i);     
    maG4x(i) = Ar4(i,1)*m4;                                                                                     
    maG4y(i) = Ar4(i,2)*m4+m4*g;                                                                                
    moment4(i) = I4*alpha_4(i,3);                                                                               
    maG5x(i) = 0;                                                                                               
    maG5y(i) = r5a(i,2)*m5+Fs(i)-Fw(i)+m5*g;                                                                    
    moment5(i) = 0;                                                                                      
    maG6x(i) = Ar6(i,1)*m6;                                                                                    
    maG6y(i) = Ar6(i,2)*m6+m6*g;                                                                                
    moment6(i) = I6*alpha_6(i,3);                                                                               
    maG7x(i) = Ar7(i,1)*m7;                                                                                    
    maG7y(i) = Ar7(i,2)*m7+m7*g;                                                                                
    moment7(i) = 0;   
    %等效力
%     
%     maG2x(i) = 0;                                                                                    
%     maG2y(i) = 0+m2*g;                                                                                
%     moment2(i) = 0;                                                                                             
%     maG3x(i) = -Fc(i)*cos(theta_d(i)+180*dtr);                                                       
%     maG3y(i) = -Fc(i)*sin(theta_d(i)+180*dtr)+m3*g;                                                  
%     moment3(i) = I3*alpha_3(i,3)+Fc(i)*cos(theta_d(i)+180*dtr)*S3y(i)-Fc(i)*sin(theta_d(i)+180*dtr)*S3x(i);     
%     maG4x(i) = 0;                                                                                     
%     maG4y(i) = 0+m4*g;                                                                                
%     moment4(i) = I4*alpha_4(i,3);                                                                               
%     maG5x(i) = 0;                                                                                               
%     maG5y(i) = 0+Fs(i)-Fw(i)+m5*g;                                                                    
%     moment5(i) = 0;                                                                                      
%     maG6x(i) = 0;                                                                                     
%     maG6y(i) = 0+m6*g;                                                                                
%     moment6(i) = I6*alpha_6(i,3);                                                                               
%     maG7x(i) = 0;                                                                                     
%     maG7y(i) =0+m7*g;                                                                                
%     moment7(i) = 0;   
    
    acceleration_matrix = [ 
        maG2x(i);
        maG2y(i);
        moment2(i);
        maG3x(i);
        maG3y(i);
        moment3(i);
        maG4x(i);
        maG4y(i);
        moment4(i);
        maG5x(i);
        maG5y(i);
        moment5(i);
        maG6x(i);
        maG6y(i);
        moment6(i);
        maG7x(i);
        maG7y(i);
        moment7(i)];
    
    Force_matrix = force\acceleration_matrix;
    
    F12x(i) = Force_matrix(1);
    F12y(i) = Force_matrix(2);
    F32x(i) = Force_matrix(3);
    F32y(i) = Force_matrix(4);
    T(i)    = Force_matrix(5);
    Fg(i)   = Force_matrix(6);
    F43x(i) = Force_matrix(7);
    F43y(i) = Force_matrix(8);
    F54x(i) = Force_matrix(9);
    F54y(i) = Force_matrix(10);
    N15(i)  = Force_matrix(11);
    N15p(i) = Force_matrix(12);
    F64x(i) = Force_matrix(13);
    F64y(i) = Force_matrix(14);
    F67x(i) = Force_matrix(15);
    F67y(i) = Force_matrix(16);
    F17x(i) = Force_matrix(17);
    F17y(i) = Force_matrix(18);
    
end


k = 0.03; 
delta_ei = 0;

% Calculate Δei
 Tavg=mean(T);% average torque
 
 
 for  i = 1:359
    delta_e = ((T(i)-Tavg)+(T(i+1)-Tavg))*dtr/2;
    
    delta_ei = delta_ei+delta_e;
    
    e(i+1) = delta_ei;
    
 end
% 
E=max(e)-min(e);
% Calculate Is
Is=E/k/(omega(1)^2);
% Calculate Ic 
Ic=I2;
% Calculate If (inertia of flywheel)
IF=(Is-Ic)/(gr^2)-IJ;





%% resultant
for i=1:360
%     F12(i) = sqrt(F12x(i)^2+F12y(i)^2);
%     F32(i) = sqrt(F32x(i)^2+F32y(i)^2);
%     F43(i) = sqrt(F43x(i)^2+F43y(i)^2);
%     F54(i) = sqrt(F54x(i)^2+F54y(i)^2);
%     F64(i) = sqrt(F64x(i)^2+F64y(i)^2);
%     F67(i) = sqrt(F67x(i)^2+F67y(i)^2);
%     F17(i) = sqrt(F17x(i)^2+F17y(i)^2);
end

%% rms
% rms_F12 = 0;
% rms_F32 = 0;
% rms_F43 = 0;
% rms_F54 = 0;
% rms_F64 = 0;
% rms_F67 = 0;
% rms_F17 = 0;
% rms_N15 = 0;
% rms_N15p = 0;
% rms_Fg = 0;

for i=1:360
%     rms_F12 = rms_F12 + F12(i)^2;
%     rms_F32 = rms_F32 + F32(i)^2;
%     rms_F43 = rms_F43 + F43(i)^2;
%     rms_F54 = rms_F54 + F54(i)^2;
%     rms_F64 = rms_F64 + F64(i)^2;
%     rms_F67 = rms_F67 + F67(i)^2;
%     rms_F17 = rms_F17 + F17(i)^2;
%     rms_N15 = rms_N15 + N15(i)^2;
%     rms_N15p = rms_N15p + N15p(i)^2;
%     rms_Fg = rms_Fg + Fg(i)^2;
      
end


%     rms_F12 = sqrt(rms_F12/360);
%     rms_F32 = sqrt(rms_F32/360);
%     rms_F43 = sqrt(rms_F43/360);
%     rms_F54 = sqrt(rms_F54/360);
%     rms_F64 = sqrt(rms_F64/360);
%     rms_F67 = sqrt(rms_F67/360);
%     rms_F17 = sqrt(rms_F17/360);
%     rms_N15 = sqrt(rms_N15/360);
%     rms_N15p = sqrt(rms_N15p/360);
%     rms_Fg = sqrt(rms_Fg/360);
      



%%
%''^{}
% 
% plot(theta_2/pi*180,atan2(F17y,F17x)/dtr,'--')
% xlim([0,360]);
% title('\theta_{17} vs \theta_2');
% xlabel('\theta_2 (\circ)'); 
% ylabel('\theta_{17} (\circ)');
% 
% 
% plot(theta_2/pi*180,T,'--');
% xlim([0,360]);
% title('T_{} vs \theta_2');
% xlabel('\theta_2 (\circ)'); 
% ylabel('T_{} (N)');
% 
% 
% 
% set(gca,'Fontsize',12,'FontName','Times New Roman','Xtick',60*(0:6),'FontWeight','Bold');
% zzz = legend('\omega_2 = 10 rad/s','\omega_2 = 60 rad/s','\omega_2 = 120 rad/s');
% set(zzz,'box','off','FontSize',8,'FontName','Times New Roman','location','best')

%% Excel

% peak_max(k,:) = [max(F12)    max(F32)    max(F43)    max(F54)    max(F64)    max(F67)    max(F17)    max(N15)    max(N15p)    max(T)      max(Fg)];  
% 
% peak_min(k,:) = [min(F12)    min(F32)    min(F43)    min(F54)    min(F64)    min(F67)    min(F17)    min(N15)    min(N15p)    min(T)      min(Fg)]; 
% 
% rms(k,:) = [rms_F12    rms_F32    rms_F43    rms_F54    rms_F64    rms_F67    rms_F17    rms_N15    rms_N15p      rms_Fg];
% 
% 
% fluctuation(k,:) = [max(F12)-min(F12)    max(F32)-min(F32)    max(F43)-min(F43)    max(F54)- min(F54)    max(F64)-min(F64)...
%                     max(F67)-min(F67)    max(F17)-min(F17)    max(N15)-min(N15)    max(N15p)-min(N15p)   max(T)-min(T)      max(Fg)-min(Fg)]; 
% 
% 
% T_mean = mean(T);
% 
% T_mean_matrix(k,:) = [T_mean];



%% POWER Equation

for i = 1:360
    
    
%     energy_Velocity(i) = m2*(Ar2(i,1)*Vr2(i,1)+Ar2(i,2)*Vr2(i,2)) + m3*(Ar3(i,1)*Vr3(i,1)+Ar3(i,2)*Vr3(i,2))+m4*(Ar4(i,1)*Vr4(i,1)+Ar4(i,2)*Vr4(i,2))...
%                         +m5*(r5a(i,2)*r5v(i,2)) + m6*(Ar6(i,1)*Vr6(i,1)+Ar6(i,2)*Vr6(i,2)) + m7*(Ar7(i,1)*Vr7(i,1)+Ar7(i,2)*Vr7(i,2));                 
%     energy_Rotation(i) = I3*omega_3(i,3)*alpha_3(i,3) + I4*omega_4(i,3)*alpha_4(i,3) + I6*omega_6(i,3)*alpha_6(i,3);
%     
%     energy_gravity(i) = -(m2*g*Vr2(i,2)+m3*g*Vr3(i,2)+m4*g*Vr4(i,2)+m5*g*r5v(i,2)+m6*g*Vr6(i,2)+m7*g*Vr7(i,2));    %with minus,because gravity direction must be equal to velocity's 
%     energy_dmaping(i) = Fc(i)*(rdv(i,2));
%     energy_spring(i) = stiffness*(spring_initial - spring_length(i))*r5v(i,2);                     %%???    
%     energy_Fw (i)= Fw(i)*r5v(i,2);
%     T_power(i) = (energy_Velocity(i) + energy_Rotation(i) - energy_gravity(i) - energy_dmaping(i) - energy_spring(i) - energy_Fw(i) )/omega_2(1,3);
    
   
end

%  plot(theta_2/pi*180,T_power);

 end