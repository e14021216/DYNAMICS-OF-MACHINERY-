clc;clear all; close all
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
t_ =  0.057;
tp = 0.013;

gr = 2.27;
sr = 1.71;

global a m2 m6 m7 phi2 phi6 phi7 I2 I6 I7 L2 L6 L7


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

IJ = 0.9;%kg-m*m
IF = 0;%kg-m*m
I2 = 0.651;
I3 = 0.047;
I4 = 0.887;
I6 = 0.313;
I7 = 0.955;

density = 7800;
thickness = 0.08;

spring_initial = 0.27;  %m

theta_1 = (360-107.47)*dtr;
theta_2i = 0;
theta_3p = 27*dtr;
theta_7i = 287.45*dtr;
theta_5 = 270*dtr;

PG = [0 -0.75 0];
PH = [0.22 0.44 0];

IJ = 0.9;
I2 = 0.651 + (IJ + IF)*(gr^2);
PO2 = [0 0 0];
PO7 = [PO2(1,1) + r1*cos(theta_1) ,PO2(1,2) + r1*sin(theta_1) , 0];

g = 9.81;
dead_center = 100;
%% Calculate 
m4b = m4*(0.322-0.1467)/0.322;

mL2x = -(m3*r2-L3/r3*cos(phi3)*r2*m3);
mL2y = L3/r3*sin(phi3)*r2*m3;

mL6x = -m4b*r6-m3*L3/r3*r6*cos(phi3);
mL6y = -m3*L3/r3*r6*sin(phi3);

mlficw2_y=mL2y-m2*L2*sin(phi2);
mlficw2_x=mL2x-m2*L2*cos(phi2);

mlficw6_y=mL6y-m6*L6*sin(phi6);
mlficw6_x=mL6x-m6*L6*cos(phi6);

rcw2 = ((mlficw2_y^2+mlficw2_x^2)/((pi*density*thickness))^2)^(1/6);
mrangle_2 = atan2(mlficw2_y,mlficw2_x);

rcw6 = ((mlficw6_y^2+mlficw6_x^2)/((pi*density*thickness))^2)^(1/6);
mrangle_6 = atan2(mlficw6_y,mlficw6_x);

mcw2=(pi*density*thickness)*rcw2^2;
mcw6=(pi*density*thickness)*rcw6^2;

mL7x = -m3*L3/r3*r7*cos(phi3)-(m6 + mcw6)*r7-m4b*r7;
mL7y = -m3*L3/r3*r7*sin(phi3);

mlficw7_y=mL7y-m7*L7*sin(phi7);
mlficw7_x=mL7x-m7*L7*cos(phi7);

rcw7 = ((mlficw7_y^2+mlficw7_x^2)/((pi*density*thickness))^2)^(1/6);
mrangle_7 = atan2(mlficw7_y,mlficw7_x);

mcw7 = (pi*density*thickness)*rcw7^2;

L2_newx = ( m2*L2*cos(phi2)+mcw2*rcw2*cos(mrangle_2) )/(m2+mcw2);
L2_newy = ( m2*L2*sin(phi2)+mcw2*rcw2*sin(mrangle_2) )/(m2+mcw2);
L2_new = (L2_newx^2 + L2_newy^2)^0.5;

L6_newx = ( m6*L6*cos(phi6)+mcw6*rcw6*cos(mrangle_6) )/(m6+mcw6);
L6_newy = ( m6*L6*sin(phi6)+mcw6*rcw6*sin(mrangle_6) )/(m6+mcw6);
L6_new = (L6_newx^2 + L6_newy^2)^0.5;

L7_newx = ( m7*L7*cos(phi7)+mcw7*rcw7*cos(mrangle_7) )/(m7+mcw7);
L7_newy = ( m7*L7*sin(phi7)+mcw7*rcw7*sin(mrangle_7) )/(m7+mcw7);
L7_new = (L7_newx^2 + L7_newy^2)^0.5;

%%
stiffness = 37000; damping = 11000;

Case_1 = input('Choose a? (a) a = 0.0027 (b) a = 0.0054 :','S');


switch Case_1
    case 'a'
        a = 0.0027;
    case 'b'
        a = 0.0054;
end




% Case_2 = input('Add spring and damper? a[x/x] b[S/x] c[x/S] d[S/D]: ','S');
% 
% 
% switch Case_2
%     case 'a'
%         stiffness = 0; damping = 0;
%     case 'b'
%         stiffness = 37000; damping = 0;
%     case 'c'
%         stiffness = 0; damping = 11000;
%     case 'd'
%         stiffness = 37000; damping = 11000;
% end

Case = input('Choose condition? (a) none (b) add link2 (c) add link6 link7 (d) add link2 link6 link7:','S');
switch Case
    case 'a'
        
    case 'b'         
        I2 = I2 + 0.5*mcw2*rcw2^2 + mcw2*((rcw2*cos(mrangle_2)-L2_newx)^2+(rcw2*sin(mrangle_2)-L2_newy)^2) + m2*((L2*cos(phi2)-L2_newx)^2+(L2*sin(phi2)-L2_newy)^2);
        L2 = L2_new;                                                        
        phi2 = atan2(L2_newy,L2_newx);
        m2 = m2 + mcw2;
    case 'c'
       
        I6 = I6 + 0.5*mcw6*rcw6^2 + mcw6*( (rcw6*cos(mrangle_6)-L6_newx)^2+(rcw6*sin(mrangle_6)-L6_newy)^2 ) +m6*( (L6*cos(phi6)-L6_newx)^2+(L6*sin(phi6)-L6_newy)^2 );
        I7 = I7 + 0.5*mcw7*rcw7^2 + mcw7*( (rcw7*cos(mrangle_7)-L7_newx)^2+(rcw7*sin(mrangle_7)-L7_newy)^2 ) +m7*( (L7*cos(phi7)-L7_newx)^2+(L7*sin(phi7)-L7_newy)^2 );
        L6 = L6_new;
        L7 = L7_new;       
        phi6 = atan2(L6_newy,L6_newx);
        phi7 = atan2(L7_newy,L7_newx);
        m6 = m6 + mcw6;        
        m7 = m7 + mcw7;
    case 'd'               
        I2 = I2 + 0.5*mcw2*rcw2^2 + mcw2*( (rcw2*cos(mrangle_2)-L2_newx)^2+(rcw2*sin(mrangle_2)-L2_newy)^2 ) +m2*( (L2*cos(phi2)-L2_newx)^2+(L2*sin(phi2)-L2_newy)^2 );
        I6 = I6 + 0.5*mcw6*rcw6^2 + mcw6*( (rcw6*cos(mrangle_6)-L6_newx)^2+(rcw6*sin(mrangle_6)-L6_newy)^2 ) +m6*( (L6*cos(phi6)-L6_newx)^2+(L6*sin(phi6)-L6_newy)^2 );
        I7 = I7 + 0.5*mcw7*rcw7^2 + mcw7*( (rcw7*cos(mrangle_7)-L7_newx)^2+(rcw7*sin(mrangle_7)-L7_newy)^2 ) +m7*( (L7*cos(phi7)-L7_newx)^2+(L7*sin(phi7)-L7_newy)^2 );
        L2 = L2_new;
        L6 = L6_new;
        L7 = L7_new;
        phi2 = atan2(mlficw2_y,mlficw2_x);
        phi6 = atan2(mlficw6_y,mlficw6_x);
        phi7 = atan2(mlficw7_y,mlficw7_x);
        m2 = m2 + mcw2;
        m6 = m6 + mcw6;
        m7 = m7 + mcw7;
        
end


%% solve ode
t_origin = 0;
t_end = 15;
interval = 0.001;

options = odeset('AbsTol',10^-6,'RelTol',10^-6);
[t,y] = ode45('overwatch',t_origin:interval:t_end,[0,0],options);

theta_2_temp = y(:,1)'; 
theta_2 = rem(y(:,1),2*pi)';

omega_2 = zeros(length(t),3);
omega_7 = zeros(length(t),3);
omega_2(:,3) = y(:,2);
omega_7(:,3) = -y(:,2);

t = t';




%% displacement


for i = 1:1:length(t)

    theta_7(i) = theta_7i-(theta_2(i)-0);
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




    PA(i,1) = PO2(1) + r2*cos(theta_2(i));
    PA(i,2) = PO2(2) + r2*sin(theta_2(i));
    PA(i,3) = 0;
    
    PB(i,1) = PA(i,1) + r3*cos(theta_3(i));
    PB(i,2) = PA(i,2) + r3*sin(theta_3(i));
    PB(i,3) = 0;
    
    PD(i,1) = PO7(1) + r7*cos(theta_7(i));
    PD(i,2) = PO7(2) + r7*sin(theta_7(i));
    PD(i,3) = 0;
    
    PR5(i,1) = PO2(1);
    PR5(i,2) = PO2(2) - r5(i);
    PR5(i,3) = 0;
    
    PC(i,1) = PB(i,1) + r4*cos(theta_4(i));
    PC(i,2) = PB(i,2) + r4*sin(theta_4(i));
    PC(i,3) = 0;
    
    PS(i,1) = PA(i,1) + r3p*cos(theta_3(i) + theta_3p);
    PS(i,2) = PA(i,2) + r3p*sin(theta_3(i) + theta_3p);
    PS(i,3) = 0;

end

for i = 1:1:length(t)
    %% h
    
    AA_ = [r3*sin(theta_3(i)) -r6*sin(theta_6(i));-r3*cos(theta_3(i)) r6*cos(theta_6(i))];
    BB = [-r2*sin(theta_2(i))-r7*sin(theta_7(i));r2*cos(theta_2(i))+r7*cos(theta_7(i))];
    X(i,:) = AA_\BB;
    h3(i) = X(i,1);
    h6(i) = X(i,2);
    omega_3(i,3) = h3(i)*omega_2(i,3);
    omega_6(i,3) = h6(i)*omega_2(i,3);
    
    
    AA_ = [-r4*sin(theta_4(i)) 0;r4*cos(theta_4(i)) 1];
    BB = [r6*h6(i)*sin(theta_6(i))-r7*sin(theta_7(i));-r6*h6(i)*cos(theta_6(i))+r7*cos(theta_7(i))];
    X(i,:) = AA_\BB;
    h4(i) = X(i,1);
    hr5(i) = X(i,2);
    omega_4(i,3) = h4(i)*omega_2(i,3);
    r5v(i,2) = -hr5(i)*omega_2(i,3);
    
    
    %% hp
    
    AA_ = [r3*sin(theta_3(i)) -r6*sin(theta_6(i));-r3*cos(theta_3(i)) r6*cos(theta_6(i))];
    BB = [h6(i)^2*r6*cos(theta_6(i))+r7*cos(theta_7(i))-r2*cos(theta_2(i))-h3(i)^2*r3*cos(theta_3(i));
        h6(i)^2*r6*sin(theta_6(i))+r7*sin(theta_7(i))-r2*sin(theta_2(i))-h3(i)^2*r3*sin(theta_3(i))];
    X(i,:) = AA_\BB;
    h3p(i) = X(i,1);
    h6p(i) = X(i,2);
   
    AA_ = [-r4*sin(theta_4(i)) 0;r4*cos(theta_4(i)) 1];
    BB = [r4*h4(i)^2*cos(theta_4(i))+r6*h6p(i)*sin(theta_6(i))+r6*h6(i)^2*cos(theta_6(i))+r7*cos(theta_7(i));
        r4*h4(i)^2*sin(theta_4(i))-r6*h6p(i)*cos(theta_6(i))+r6*h6(i)^2*sin(theta_6(i))+r7*sin(theta_7(i))];
    X(i,:) = AA_\BB;
    h4p(i) = X(i,1);
    hr5p(i) = X(i,2);
    
    
          %% spring_length and damper velocity
    rd(i) = sqrt((PS(i,1)-PH(1,1))^2+(PS(i,2)-PH(1,2))^2);
    spring_length(i) = (PO2(1,2)-PG(1,2)) -r5(i);
    
    
    theta_d(i) = atan2(PH(1,2)-PS(i,2),PH(1,1)-PS(i,1));
    
    AA_ = [cos(theta_d(i)) -rd(i)*sin(theta_d(i));sin(theta_d(i)) rd(i)*cos(theta_d(i))];
    BB = [r2*sin(theta_2(i))+h3(i)*r3p*sin(theta_3(i) + theta_3p);-r2*cos(theta_2(i))-h3(i)*r3p*cos(theta_3(i) + theta_3p)];
    
    X(i,:) = AA_\BB;
    
    fd(i) = X(i,1);
    hd(i) = X(i,2);
    
    rdv(i,2) = fd(i)*omega_2(i,3);
    rdv(i,3) = 0;
    if dead_center > spring_length(i)
        dead_center = spring_length(i);
        
    end
   
    
     %% kinematic coefficient
    f2x(i) = -L2*sin(theta_2(i));
    f2y(i) = L2*cos(theta_2(i));
    
    f2xp(i) = -L2*cos(theta_2(i));
    f2yp(i) = -L2*sin(theta_2(i));
    
    f3x(i) = -r2*sin(theta_2(i))-L3*h3(i)*sin(theta_3(i)+phi3);
    f3y(i) = r2*cos(theta_2(i))+L3*h3(i)*cos(theta_3(i)+phi3);
    
    f3xp(i) = -r2*cos(theta_2(i))-L3*h3p(i)*sin(theta_3(i)+phi3)-L3*(h3(i)^2)*cos(theta_3(i)+phi3);
    f3yp(i) = -r2*sin(theta_2(i))+L3*h3p(i)*cos(theta_3(i)+phi3)-L3*(h3(i)^2)*sin(theta_3(i)+phi3);
    
    f4x(i) = -r2*sin(theta_2(i))-h3(i)*r3*sin(theta_3(i))-L4*h4(i)*sin(theta_4(i)+phi4);
    f4y(i) = r2*cos(theta_2(i))+h3(i)*r3*cos(theta_3(i))+L4*h4(i)*cos(theta_4(i)+phi4);
    
    f4xp(i) = -r2*cos(theta_2(i))-(h3(i)^2)*r3*cos(theta_3(i))-h3p(i)*r3*sin(theta_3(i))-(h4(i)^2)*L4*cos(theta_4(i))-h4p(i)*L4*sin(theta_4(i));
    f4yp(i) = -r2*sin(theta_2(i))-(h3(i)^2)*r3*sin(theta_3(i))+h3p(i)*r3*cos(theta_3(i))-(h4(i)^2)*L4*sin(theta_4(i))+h4p(i)*L4*cos(theta_4(i));
    
    
    f5y(i) = -hr5(i);        %f5y , positive when go up
    f5yp(i)= -hr5p(i);

    f6x(i) = r7*sin(theta_7(i))-h6(i)*L6*sin(theta_6(i));
    f6y(i) = -r7*cos(theta_7(i))+h6(i)*L6*cos(theta_6(i));
    
    f6xp(i) = -r7*cos(theta_7(i))-(h6(i)^2)*L6*cos(theta_6(i))-h6p(i)*L6*sin(theta_6(i));
    f6yp(i) = -r7*sin(theta_7(i))-(h6(i)^2)*L6*sin(theta_6(i))+h6p(i)*L6*cos(theta_6(i));
    
    f7x(i) = L7*sin(theta_7(i));
    f7y(i) = -L7*cos(theta_7(i));
    
    f7xp(i) = -L7*cos(theta_7(i));
    f7yp(i) = -L7*sin(theta_7(i));
    
    %% Link center mass displacement
    

    
    %% Displacement vector
    for j = 1:3
        vector_O2_A(i,j) = PA(i,j)-PO2(1,j);
        vector_A_B(i,j) = PB(i,j)-PA(i,j);
        vector_B_C(i,j) = PC(i,j)-PB(i,j);
        vector_O7_D(i,j) = PD(i,j)-PO7(1,j);
        vector_D_B(i,j) = PB(i,j)-PD(i,j);
        
    end
    
   
    %% Linkage mass-center displacement

    PG_r2(i,1) = PO2(1,1) + L2*cos(theta_2(i)+phi2);
    PG_r2(i,2) = PO2(1,2) + L2*sin(theta_2(i)+phi2);
    PG_r2(i,3) = 0;
    
    PG_r3(i,1) = PA(i,1) + L3*cos(theta_3(i)+phi3);
    PG_r3(i,2) = PA(i,2) + L3*sin(theta_3(i)+phi3);
    PG_r3(i,3) = 0;
    
    PG_r4(i,1) = PB(i,1) + L4*cos(theta_4(i)+phi4);
    PG_r4(i,2) = PB(i,2) + L4*sin(theta_4(i)+phi4);
    PG_r4(i,3) = 0;
    
    PG_r5(i,1) = PB(i,1) + r4*cos(theta_4(i));
    PG_r5(i,2) = PB(i,2) + r4*sin(theta_4(i));
    PG_r5(i,3) = 0; 
    
    PG_r6(i,1) = PD(i,1) + L6*cos(theta_6(i)+phi6);
    PG_r6(i,2) = PD(i,2) + L6*sin(theta_6(i)+phi6);
    PG_r6(i,3) = 0;
    
    PG_r7(i,1) = PO7(1,1) + L7*cos(theta_7(i)+phi7);
    PG_r7(i,2) = PO7(1,2) + L7*sin(theta_7(i)+phi7);
    PG_r7(i,3) = 0;
     
     %% Displacement center-mass vector
    
    for j = 1:3
        vector_O2_Gr2(i,j) = PG_r2(i,j)-PO2(1,j);
        vector_A_Gr3(i,j) = PG_r3(i,j)-PA(i,j);
        vector_B_Gr4(i,j) = PG_r4(i,j)-PB(i,j);
        vector_O7_Gr7(i,j) = PG_r7(i,j)-PO7(1,j);
        vector_D_Gr6(i,j) = PG_r6(i,j)-PD(i,j);
        
    end
    
    %% Linkage velocity

    VO2 = [0 0 0];
    VO7 = [0 0 0];
    VA(i,:) = cross(omega_2(i,:),vector_O2_A(i,:));
    VB(i,:) = cross(omega_3(i,:),vector_A_B(i,:)) + VA(i,:);
    VD(i,:) = cross(omega_7(i,:),vector_O7_D(i,:)); 
    



    %% Link center mass velocity
    Vr2(i,:) = cross(omega_2(i,:),vector_O2_Gr2(i,:));
    Vr3(i,:) = cross(omega_3(i,:),vector_A_Gr3(i,:))+VA(i,:);
    Vr4(i,:) = cross(omega_4(i,:),vector_B_Gr4(i,:))+VB(i,:);
    Vr6(i,:) = cross(omega_6(i,:),vector_D_Gr6(i,:))+VD(i,:);
    Vr7(i,:) = cross(omega_7(i,:),vector_O7_Gr7(i,:));
    
    %% Force
    Fs(i) = stiffness*(spring_length(i)-spring_initial);     
    Fc(i) = rdv(i,2)*damping;
    
    if spring_length(i)-dead_center <= 0.005 && r5v(i,2) < 0
        Fw(i) = 89000;
%         Fw(i) = 0;
    else
        Fw(i) = 0;
    end
    
     %% alpha_2 alpha_7   
   sigmaA(i) = m2*(f2x(i)^2+f2y(i)^2) + I2*1 + m3*(f3x(i)^2+f3y(i)^2) + I3*h3(i)^2 + m4*(f4x(i)^2+f4y(i)^2)+I4*h4(i)^2 ...
        + m6*(f6x(i)^2+f6y(i)^2)+ I6*h6(i)^2 + m7*(f7x(i)^2+f7y(i)^2) + I7*(-1)^2 + m5*(f5y(i)^2);
    sigmaB = m2*(f2x(i)*f2xp(i)+f2y(i)*f2yp(i)) + m3*(f3x(i)*f3xp(i)+f3y(i)*f3yp(i))+I3*h3(i)*h3p(i) + m4*(f4x(i)*f4xp(i)+f4y(i)*f4yp(i))+I4*h4(i)*h4p(i) ...
        + m6*(f6x(i)*f6xp(i)+f6y(i)*f6yp(i))+I6*h6(i)*h6p(i) + m7*(f7x(i)*f7xp(i)+f7y(i)*f7yp(i)) + m5*(f5y(i)*f5yp(i));
    sigmaUg = -(m2*f2y(i)+m3*f3y(i)+m4*f4y(i)+m5*f5y(i)+m6*f6y(i)+m7*f7y(i))*g;
    Uspring = stiffness*(spring_initial-spring_length(i))*f5y(i) ;
    Udamping = -damping*fd(i)^2;
    UFw = Fw(i)*f5y(i);
    
    Motor_rpm = omega_2(i,3)*gr*sr*60/(2*pi);
    
    if Motor_rpm >= 0 && Motor_rpm <=1500
        Torque(i) = a*(27000-7*Motor_rpm)*gr*sr;
        
    elseif Motor_rpm >= 1500 && Motor_rpm <=1800
        Torque(i) = a*(16500-55*(Motor_rpm-1500))*gr*sr;
        
    else Torque(i) = 0;

    end
    
    alpha_2(i,3) = (Torque(i) - sigmaB*omega_2(i,3)^2 + sigmaUg + Uspring + Udamping*omega_2(i,3) + UFw)/(sigmaA(i));
    alpha_7(i,3) = -alpha_2(i,3);
    
    alpha_3(i,3) = h3p(i)*omega_2(i,3)^2 + h3(i)*alpha_2(i,3);
    alpha_6(i,3) = h6p(i)*omega_2(i,3)^2 + h6(i)*alpha_2(i,3);
    
    alpha_4(i,3) = h4p(i)*omega_2(i,3)^2 + h4(i)*alpha_2(i,3);
    r5a(i,2) = -(hr5p(i)*omega_2(i,3)^2 + hr5(i)*alpha_2(i,3));
    
    

      %%  Linkage acceleration
    
    AA(i,:) = cross(omega_2(i,:),VA(i,:)) + cross(alpha_2(i,:),vector_O2_A(i,:));
    AB(i,:) = AA(i,:) + cross(omega_3(i,:),cross(omega_3(i,:),vector_A_B(i,:))) + cross(alpha_3(i,:),vector_A_B(i,:));
    AC(i,:) = r5a(i,:);
    AD(i,:) = cross(omega_7(i,:),VD(i,:)) + cross(alpha_7(i,:),vector_O7_D(i,:));
    
    %% Linkage mass acceleration
   
    Ar2(i,:) = cross(omega_2(i,:),Vr2(i,:)) +  cross(alpha_2(i,:),vector_O2_Gr2(i,:));
    Ar3(i,:) = AA(i,:) + cross(omega_3(i,:),cross(omega_3(i,:),vector_A_Gr3(i,:))) + cross(alpha_3(i,:),vector_A_Gr3(i,:));
    Ar4(i,:) = AB(i,:) + cross(omega_4(i,:),cross(omega_4(i,:),vector_B_Gr4(i,:))) + cross(alpha_4(i,:),vector_B_Gr4(i,:));
    Ar6(i,:) = AD(i,:) + cross(omega_6(i,:),cross(omega_6(i,:),vector_D_Gr6(i,:))) + cross(alpha_6(i,:),vector_D_Gr6(i,:));
    Ar7(i,:) = cross(omega_7(i,:),Vr7(i,:)) +  cross(alpha_7(i,:),vector_O7_Gr7(i,:));
   
end
    
    LDCP = [0,PG(2)+dead_center,0];
    
    for j = 1:3
        vector_O2_LDCP(j) = PO2(j)-LDCP(j);
        vector_O7_LDCP(j) = PO7(j)-LDCP(j);
        vector_H_LDCP(j) = PH(j)-LDCP(j);
    end
  
    %% 18*18 matrix
    for i= 1:1:length(t)
         
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
   
    Fg2x = cos(theta_1+90*dtr-pressure_angle);
    Fg2y = sin(theta_1+90*dtr-pressure_angle);
    Fg7x = -cos(theta_1+90*dtr-pressure_angle);
    Fg7y = -sin(theta_1+90*dtr-pressure_angle);
    
    g1(i) = Fg2x*(r1/2*sin(theta_1+pi)+L2*sin(theta_2(i)))-Fg2y*(r1/2*cos(theta_1+pi)+L2*cos(theta_2(i)));   % WHY
    g2(i) = -Fg7x*(r1/2*sin(theta_1+pi)-L7*sin(theta_7(i)))-Fg2y*(r1/2*cos(theta_1+pi)-L7*cos(theta_7(i)));  % WHY
    
 
 %%force   1       2       3       4       5       6       7       8       9       10      11      12      13      14      15      16      17      18          Force
 force = [ -1      0       1       0       0       Fg2x    0       0       0       0       0       0       0       0       0       0       0       0     ;     %F12x     1   Link2x
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
           0       0       0       0       0       0       0       0       0       0       t_      -tp     0       0       0       0       0       0     ;     %N15p     12  Link5z
           0       0       0       0       0       0       0       0       0       0       0       0       -1      0       -1      0       0       0     ;     %F64x     13  Link6x
           0       0       0       0       0       0       0       0       0       0       0       0       0       1       0       1       0       0     ;     %F64y     14  Link6y
           0       0       0       0       0       0       0       0       0       0       0       0       B6y(i)  B6x(i)  D6y(i)  D6x(i)  0       0     ;     %F67x     15  Link6z
           0       0       0       0       0       Fg7x    0       0       0       0       0       0       0       0       1       0       -1      0     ;     %F67y     16  Link7x
           0       0       0       0       0       Fg7y    0       0       0       0       0       0       0       0       0       -1      0       1     ;     %F17x     17  Link7y
           0       0       0       0       0       g2(i)   0       0       0       0       0       0       0       0       -D7y(i) -D7x(i) O7y(i)  O7x(i)];    %F17y     18  Link7z
       
       
       maG2x(i) = Ar2(i,1)*m2;
       maG2y(i) = Ar2(i,2)*m2+m2*g;
       moment2(i) = I2*alpha_2(i,3);
       maG3x(i) = Ar3(i,1)*m3 - Fc(i)*cos(theta_d(i));
       maG3y(i) = Ar3(i,2)*m3 - Fc(i)*sin(theta_d(i))+m3*g;
       moment3(i) = I3*alpha_3(i,3)+Fc(i)*cos(theta_d(i))*S3y(i)-Fc(i)*sin(theta_d(i))*S3x(i);
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
       moment7(i) = I7*alpha_7(i,3);
       
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
       
       Force_matrix = inv(force)*acceleration_matrix;
       
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
       
       F12m(i,:) = cross(vector_O2_LDCP(1,:),[F12x(i),-F12y(i),0]);
       F17m(i,:) = cross(vector_O7_LDCP(1,:),[F17x(i),-F17y(i),0]);
       Fcm(i,:) = cross(vector_H_LDCP(1,:),[-Fc(i)*cos(theta_d(i)),-Fc(i)*sin(theta_d(i)),0]);
    
       Shaking_x(i) = F12x(i) + F17x(i) - Fc(i)*cos(theta_d(i)) - N15(i) - N15p(i);
       Shaking_y(i) = -F12y(i) - F17y(i) - Fc(i)*sin(theta_d(i)) - stiffness*(spring_initial-spring_length(i)) - Fw(i); 
       Shaking_result(i) = sqrt(Shaking_x(i)^2+Shaking_y(i)^2);
       Shaking_M(i) = -T(i) + F12m(i,3) + F17m(i,3) + Fcm(i,3);
       
       
       
       
%        F12(i)=sqrt(F12x(i)^2+F12y(i)^2);
%        F32(i)=sqrt(F32x(i)^2+F32y(i)^2);
%        F43(i)=sqrt(F43x(i)^2+F43y(i)^2);
%        F64(i)=sqrt(F64x(i)^2+F64y(i)^2);
%        F54(i)=sqrt(F54x(i)^2+F54y(i)^2);
% 
% 
%        F17(i)=sqrt(F17x(i)^2+F17y(i)^2);
%        F67(i)=sqrt(F67x(i)^2+F67y(i)^2);
%       
      
       %take time = 3-5 Torque data
       if t(i) > 4 && t(i) < 6
           temp_data_1(i) = Torque(i);
           
       end
       %take time = 8-9 Torque data
       if t(i) > 12 && t(i) < 14
           temp_data_2(i) = Torque(i);
           
       end
       
       
       

%     energy_Velocity(i) = m2*(Ar2(i,1)*Vr2(i,1)+Ar2(i,2)*Vr2(i,2)) + m3*(Ar3(i,1)*Vr3(i,1)+Ar3(i,2)*Vr3(i,2))+m4*(Ar4(i,1)*Vr4(i,1)+Ar4(i,2)*Vr4(i,2))...
%                         +m5*(r5a(i,2)*r5v(i,2)) + m6*(Ar6(i,1)*Vr6(i,1)+Ar6(i,2)*Vr6(i,2)) + m7*(Ar7(i,1)*Vr7(i,1)+Ar7(i,2)*Vr7(i,2));
%
%     energy_Rotation(i) = I2*omega_2(i,3)*alpha_2(i,3) + I3*omega_3(i,3)*alpha_3(i,3) + I4*omega_4(i,3)*alpha_4(i,3) + I6*omega_6(i,3)*alpha_6(i,3) +  I7*omega_7(i,3)*alpha_7(i,3);
%
%     energy_gravity(i) = -(m2*g*Vr2(i,2)+m3*g*Vr3(i,2)+m4*g*Vr4(i,2)+m5*g*r5v(i,2)+m6*g*Vr6(i,2)+m7*g*Vr7(i,2));    %with minus,because gravity direction must be equal to velocity's
%     energy_dmaping(i) = Fc(i)*(rdv(i,2));
%     energy_spring(i) = stiffness*(spring_initial - spring_length(i))*r5v(i,2);
%     energy_Fw (i)= Fw(i)*r5v(i,2);
%     T_power(i) = (energy_Velocity(i) + energy_Rotation(i) - energy_gravity(i) - energy_dmaping(i) - energy_spring(i) - energy_Fw(i) )/omega_2(i,3);
    end
    
    
    
       
    %find peak
    peak_1 = findpeaks(temp_data_1);
    %find peak
    peak_1 = findpeaks(peak_1);
    
    
    
    for  i = 4000:1:6000
        for j = 1:1:length(peak_1)
            
            if peak_1(j) == Torque(i)
                temp_t_1(j) = t(i);
                
            end
        end
    end
    
    
    period_1 = (max(temp_t_1)-min(temp_t_1))/(length(temp_t_1)-1);
    period_1 = round(period_1,4);
    
   
    
    %find peak
    peak_2 = findpeaks(temp_data_2);
    %find peak
    peak_2 = findpeaks(peak_2);


 
 for  i = 12000:1:14000
     for j = 1:1:length(peak_2)

         if peak_2(j) == Torque(i)
             temp_t_2(j) = t(i);
             
         end
     end
 end
    
 period_2 = (max(temp_t_2)-min(temp_t_2))/(length(temp_t_2)-1);
 period_2 = round(period_2,4);
%T_power not same



%% Excel

OP_1 = round(min(temp_t_1),3)*1000;
OP_2 = round((min(temp_t_1)+period_1),3)*1000;

OP_3 = round(min(temp_t_2),3)*1000;
OP_4 = round((min(temp_t_2)+period_2),3)*1000;


for  i = uint64(OP_1):1:uint64(OP_2)
    omega_2_st(i-OP_1+1) = omega_2(i,3);
    Torque_st(i-OP_1+1) = Torque(i);
    Shaking_M_st(i-OP_1+1) = Shaking_M(i);
    Shaking_result_st(i-OP_1+1) = Shaking_result(i);
end
OP_3 = round(OP_3,3);


for  j = uint64(OP_3):1:uint64(OP_4)
    omega_2_nd(j-OP_3+1) = omega_2(j,3);
    Torque_nd(j-OP_3+1) = Torque(j);
    Shaking_M_nd(j-OP_3+1) = Shaking_M(j);
    Shaking_result_nd(j-OP_3+1) = Shaking_result(j);
end





  for  i = uint64(min(temp_t_1)*1000):1:uint64((min(temp_t_1)+period_2)*1000)
       omega_2_st(i-3999) = omega_2(i,3);
       Torque_st(i-3999) = Torque(i);
       Shaking_M_st(i-3999) = Shaking_M(i);
       Shaking_result_st(i-3999) = Shaking_result(i);
 end

 for  i = uint64(min(temp_t_2)*1000):1:uint64((min(temp_t_2)+period_2)*1000)
       omega_2_nd(i-11999) = omega_2(i,3);
       Torque_nd(i-11999) = Torque(i);
       Shaking_M_nd(i-11999) = Shaking_M(i);
       Shaking_result_nd(i-11999) = Shaking_result(i);
 end
 
fluctuation_omega_2_st = max(omega_2_st)-min(omega_2_st);
fluctuation_omega_2_nd = max(omega_2_nd)-min(omega_2_nd);
average_omega_2_st = mean(omega_2_st);
average_omega_2_nd = mean(omega_2_nd);
coefficient_omega_2_st = fluctuation_omega_2_st/average_omega_2_st;
coefficient_omega_2_nd = fluctuation_omega_2_nd/average_omega_2_nd;
%
fluctuation_Torque_st = max(Torque_st)-min(Torque_st);
fluctuation_Torque_nd = max(Torque_nd)-min(Torque_nd);
average_Torque_st = mean(Torque_st);
average_Torque_nd = mean(Torque_nd);
coefficient_Torque_st = fluctuation_Torque_st/average_Torque_st;
coefficient_Torque_nd = fluctuation_Torque_nd/average_Torque_nd;

%
fluctuation_Shaking_result_st = max(Shaking_result_st)-min(Shaking_result_st);
fluctuation_Shaking_result_nd = max(Shaking_result_nd)-min(Shaking_result_nd);
average_Shaking_result_st = mean(Shaking_result_st);
average_Shaking_result_nd = mean(Shaking_result_nd);
coefficient_Shaking_result_st = fluctuation_Shaking_result_st/average_Shaking_result_st;
coefficient_Shaking_result_nd = fluctuation_Shaking_result_nd/average_Shaking_result_nd;
%
fluctuation_Shaking_M_st = max(Shaking_M_st)-min(Shaking_M_st);
fluctuation_Shaking_M_nd = max(Shaking_M_nd)-min(Shaking_M_nd);
average_Shaking_M_st = mean(Shaking_M_st);
average_Shaking_M_nd = mean(Shaking_M_nd);
coefficient_Shaking_M_st = fluctuation_Shaking_M_st/average_Shaking_M_st;
coefficient_Shaking_M_nd = fluctuation_Shaking_M_nd/average_Shaking_M_nd;

%                 1 steady state                 2 steady state
data_matrix = [   max(omega_2_st)                max(omega_2_nd)              ; 
                  min(omega_2_st)                min(omega_2_nd)              ; 
                  fluctuation_omega_2_st         fluctuation_omega_2_nd       ;
                  average_omega_2_st             average_omega_2_nd           ;
                  coefficient_omega_2_st         coefficient_omega_2_nd       ;
                  
                  max(Torque_st)                 max(Torque_nd)               ; 
                  min(Torque_st)                 min(Torque_nd)               ; 
                  fluctuation_Torque_st          fluctuation_Torque_nd        ;
                  average_Torque_st              average_Torque_nd            ;
                  coefficient_Torque_st          coefficient_Torque_nd        ;
                   
                  max(Shaking_result_st)         max(Shaking_result_nd)       ;
                  min(Shaking_result_st)         min(Shaking_result_nd)       ; 
                  fluctuation_Shaking_result_st  fluctuation_Shaking_result_nd;
                  average_Shaking_result_st      average_Shaking_result_nd    ;
                  coefficient_Shaking_result_st  coefficient_Shaking_result_nd;
                  
                  max(Shaking_M_st)              max(Shaking_M_nd)            ; 
                  min(Shaking_M_st)              min(Shaking_M_nd)            ; 
                  fluctuation_Shaking_M_st       fluctuation_Shaking_M_nd     ;
                  average_Shaking_M_st           average_Shaking_M_nd         ;
                  coefficient_Shaking_M_st       coefficient_Shaking_M_nd     ;

];

%% plot time
% % inertia
plot_number = 0;
% 
% %omega 2
% plot_number = plot_number + 1; set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
% figure(plot_number)
% plot(t,omega_2(:,3));
% figure_chienli2(300,260,0.23,0.25,0.7,0.6)
% title('\omega_{2} v.s Time');
% xlabel('Time (s)'); 
% ylabel('\omega_{2} (rad/s)');
% 
% %Torque
plot_number = plot_number + 1; set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
figure(plot_number)
figure_chienli2(300,260,0.23,0.25,0.7,0.6)
plot(t,Torque);
title('1^{st} s.s. Input Torque');
xlabel('Time (s)'); 
ylabel('Torque (Nm)');
xlim([min(temp_t_1) , min(temp_t_1) + period_1]);
% 
% plot_number = plot_number + 1; set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
% figure(plot_number)
% figure_chienli2(300,260,0.23,0.25,0.7,0.6)
% plot(t,Torque);
% title('2^{nd} s.s. Input Torque');
% xlabel('Time (s)'); 
% ylabel('Torque (Nm)');
% xlim([min(temp_t_2) , min(temp_t_2) + period_2]);
% %F12
% plot_number = plot_number + 1; set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
% figure(plot_number)
% figure_chienli2(300,260,0.23,0.25,0.7,0.6)
% plot(t,F12x,t,F12y,'--');
% title('1^{st} s.s. F_{12}  ');
% xlabel('Time (s)'); 
% ylabel('Force (N)');
% xlim([min(temp_t_1) , min(temp_t_1) + period_1]);
% zzz = legend('F_{12x}','F_{12y}');
% set(zzz,'box','off','FontSize',10,'FontName','Times New Roman','location','best')
% 
% plot_number = plot_number + 1; set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
% figure(plot_number)
% figure_chienli2(300,260,0.23,0.25,0.7,0.6)
% plot(t,F12x,t,F12y,'--');
% title('2^{nd} s.s. F_{12}  ');
% xlabel('Time (s)'); 
% ylabel('Force (N)');
% xlim([min(temp_t_2) , min(temp_t_2) + period_2]);
% zzz = legend('F_{12x}','F_{12y}');
% set(zzz,'box','off','FontSize',10,'FontName','Times New Roman','location','best')
% %F32
% plot_number = plot_number + 1; set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
% figure(plot_number)
% figure_chienli2(300,260,0.23,0.25,0.7,0.6)
% plot(t,F32x,t,F32y,'--');
% title('1^{st} s.s. F_{32}  ');
% xlabel('Time (s)'); 
% ylabel('Force (N)');
% xlim([min(temp_t_1) , min(temp_t_1) + period_1]);
% zzz = legend('F_{32x}','F_{32y}');
% set(zzz,'box','off','FontSize',10,'FontName','Times New Roman','location','best')
% 
% plot_number = plot_number + 1; set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
% figure(plot_number)
% figure_chienli2(300,260,0.23,0.25,0.7,0.6)
% plot(t,F32x,t,F32y,'--');
% title('2^{nd} s.s. F_{32}  ');
% xlabel('Time (s)'); 
% ylabel('Force (N)');
% xlim([min(temp_t_2) , min(temp_t_2) + period_2]);
% zzz = legend('F_{32x}','F_{32y}');
% set(zzz,'box','off','FontSize',10,'FontName','Times New Roman','location','best')
% %F43
% plot_number = plot_number + 1; set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
% figure(plot_number)
% figure_chienli2(300,260,0.23,0.25,0.7,0.6)
% plot(t,F43x,t,F43y,'--');
% title('1^{st} s.s. F_{43}  ');
% xlabel('Time (s)'); 
% ylabel('Force (N)');
% xlim([min(temp_t_1) , min(temp_t_1) + period_1]);
% zzz = legend('F_{43x}','F_{43y}');
% set(zzz,'box','off','FontSize',10,'FontName','Times New Roman','location','best')
% 
% plot_number = plot_number + 1; set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
% figure(plot_number)
% figure_chienli2(300,260,0.23,0.25,0.7,0.6)
% plot(t,F43x,t,F43y,'--');
% title('2^{nd} s.s. F_{43}  ');
% xlabel('Time (s)'); 
% ylabel('Force (N)');
% xlim([min(temp_t_2) , min(temp_t_2) + period_2]);
% zzz = legend('F_{43x}','F_{43y}');
% set(zzz,'box','off','FontSize',10,'FontName','Times New Roman','location','best')
% 
% %F54
% plot_number = plot_number + 1; set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
% figure(plot_number)
% figure_chienli2(300,260,0.23,0.25,0.7,0.6)
% plot(t,F54x,t,F54y,'--');
% title('1^{st} s.s. F_{54}  ');
% xlabel('Time (s)'); 
% ylabel('Force (N)');
% xlim([min(temp_t_1) , min(temp_t_1) + period_1]);
% zzz = legend('F_{54x}','F_{54y}');
% set(zzz,'box','off','FontSize',10,'FontName','Times New Roman','location','best')
% 
% plot_number = plot_number + 1; set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
% figure(plot_number)
% figure_chienli2(300,260,0.23,0.25,0.7,0.6)
% plot(t,F54x,t,F54y,'--');
% title('2^{nd} s.s. F_{54}  ');
% xlabel('Time (s)'); 
% ylabel('Force (N)');
% xlim([min(temp_t_2) , min(temp_t_2) + period_2]);
% zzz = legend('F_{54x}','F_{54y}');
% set(zzz,'box','off','FontSize',10,'FontName','Times New Roman','location','best')
% 
% %F64
% plot_number = plot_number + 1; set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
% figure(plot_number)
% figure_chienli2(300,260,0.23,0.25,0.7,0.6)
% plot(t,F64x,t,F64y,'--');
% title('1^{st} s.s. F_{64}  ');
% xlabel('Time (s)'); 
% ylabel('Force (N)');
% xlim([min(temp_t_1) , min(temp_t_1) + period_1]);
% zzz = legend('F_{64x}','F_{64y}');
% set(zzz,'box','off','FontSize',10,'FontName','Times New Roman','location','best')
% 
% plot_number = plot_number + 1; set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
% figure(plot_number)
% figure_chienli2(300,260,0.23,0.25,0.7,0.6)
% plot(t,F64x,t,F64y,'--');
% title('2^{nd} s.s. F_{64}  ');
% xlabel('Time (s)'); 
% ylabel('Force (N)');
% xlim([min(temp_t_2) , min(temp_t_2) + period_2]);
% zzz = legend('F_{64x}','F_{64y}');
% set(zzz,'box','off','FontSize',10,'FontName','Times New Roman','location','best')
% %F67
% plot_number = plot_number + 1; set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
% figure(plot_number)
% figure_chienli2(300,260,0.23,0.25,0.7,0.6)
% plot(t,F67x,t,F67y,'--');
% title('1^{st} s.s. F_{67}  ');
% xlabel('Time (s)'); 
% ylabel('Force (N)');
% xlim([min(temp_t_1) , min(temp_t_1) + period_1]);
% zzz = legend('F_{67x}','F_{67y}');
% set(zzz,'box','off','FontSize',10,'FontName','Times New Roman','location','best')
% 
% plot_number = plot_number + 1; set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
% figure(plot_number)
% figure_chienli2(300,260,0.23,0.25,0.7,0.6)
% plot(t,F67x,t,F67y,'--');
% title('2^{nd} s.s. F_{67}  ');
% xlabel('Time (s)'); 
% ylabel('Force (N)');
% xlim([min(temp_t_2) , min(temp_t_2) + period_2]);
% zzz = legend('F_{67x}','F_{67y}');
% set(zzz,'box','off','FontSize',10,'FontName','Times New Roman','location','best')
% %F17
% plot_number = plot_number + 1; set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
% figure(plot_number)
% figure_chienli2(300,260,0.23,0.25,0.7,0.6)
% plot(t,F17x,t,F17y,'--');
% title('1^{st} s.s. F_{17}  ');
% xlabel('Time (s)'); 
% ylabel('Force (N)');
% xlim([min(temp_t_1) , min(temp_t_1) + period_1]);
% zzz = legend('F_{17x}','F_{17y}');
% set(zzz,'box','off','FontSize',10,'FontName','Times New Roman','location','best')
% 
% plot_number = plot_number + 1; set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
% figure(plot_number)
% figure_chienli2(300,260,0.23,0.25,0.7,0.6)
% plot(t,F17x,t,F17y,'--');
% title('2^{nd} s.s. F_{17}  ');
% xlabel('Time (s)'); 
% ylabel('Force (N)');
% xlim([min(temp_t_2) , min(temp_t_2) + period_2]);
% zzz = legend('F_{17x}','F_{17y}');
% set(zzz,'box','off','FontSize',10,'FontName','Times New Roman','location','best')
% %Shaking Moment
% plot_number = plot_number + 1; set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
% figure(plot_number)
% figure_chienli2(300,260,0.23,0.25,0.7,0.6)
% plot(t,Shaking_M);
% % title('  1^{st} s.s. Shaking Moment');
% title('  1^{st} s.s. Shaking_M');
% xlabel('Time (s)'); 
% ylabel('Moment (Nm)');
% xlim([min(temp_t_1) , min(temp_t_1) + period_1]);
% 
% plot_number = plot_number + 1; set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
% figure(plot_number)
% figure_chienli2(300,260,0.23,0.25,0.7,0.6)
% plot(t,Shaking_M);
% % title('  2^{nd} s.s. Shaking Moment');
% title('  2^{nd} s.s. Shaking_M');
% xlabel('Time (s)'); 
% ylabel('Moment (Nm)');
% xlim([min(temp_t_2) , min(temp_t_2) + period_2]);
% set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
% % 
% for i = plot_number:-1:1
%     
%     figure(i)
%     
% end
% % 
% plot_number = plot_number + 1; set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
% figure(plot_number)
% title('  1^{st} s.s. S.F.');
% xlabel('x Component (N)'); 
% ylabel('y Component (N)');
% figure_chienli2(300,260,0.23,0.25,0.7,0.6);
% set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
% 
% hold all
% for k = min(temp_t_1)*1000:1:(min(temp_t_1) + period_1)*1000
%         plot([0,Shaking_x(k)],[0,Shaking_y(k)],'-');
% end
% plot(Shaking_x(min(temp_t_1)*1000:1:(min(temp_t_1) + period_1)*1000),Shaking_y(min(temp_t_1)*1000:1:(min(temp_t_1) + period_1)*1000),'o');
% 
% 
% 
% plot_number = plot_number + 1; set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
% figure(plot_number)
% title('  2^{nt} s.s. S.F.');
% xlabel('x Component (N)'); 
% ylabel('y Component (N)');
% figure_chienli2(300,260,0.23,0.25,0.7,0.6);
% set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
% 
% hold all
% for k = uint64(min(temp_t_2)*1000):1:uint64((min(temp_t_2) + period_2)*1000)
%         plot([0,Shaking_x(k)],[0,Shaking_y(k)],'-');
% end
% 
% 
% plot(Shaking_x( uint64(min(temp_t_2)*1000):1:   uint64((min(temp_t_2) + period_2)*1000)) ,Shaking_y(  uint64(min(temp_t_2)*1000)  :1: uint64((min(temp_t_2) + period_2)*1000)),'o');
% 
% hold all
% 
% for i = plot_number:-1:1
%     
%     figure(i)
%     
% end

% for k =  min(temp_t_2)*1000:1: (min(temp_t_2) + period_2)*1000
%         plot([0,Shaking_x(k)],[0,Shaking_y(k)],'-');
% end
% plot(Shaking_x(min(temp_t_2)*1000:1:(min(temp_t_2) + period_2)*1000),Shaking_y(min(temp_t_2)*1000:1:(min(temp_t_2) + period_2)*1000),'o');
% 
% 
% file = load('C:\Users\user\Desktop\omega2.txt');
% file = load('C:\Users\user\Desktop\omega2.txt');
% data_1 = file(:,1);
% data_2 = file(:,2);
% 
% 
% figure_chienli2(700,450,0.23,0.25,0.7,0.6)
% 
% plot(t,omega_2(:,3));
% hold on
% plot(data_1,data_2,'--r','Linewidth',2);
% 
% title('\omega_{2} v.s Time');
% xlabel('Time (s)'); %bf 粗體 /it 斜體 /rm 恢復正常
% ylabel('\omega_{2} (rad/s)');
% zzz = legend(' Matlab',' Recurdyn');
% set(zzz,'box','off','FontSize',20,'FontName','Times New Roman','location','best')
% 
% set(gca,'Fontsize',20,'FontName','Times New Roman','FontWeight','Bold');
% 
% 
% 
% title('Torque v.s Time');
% xlabel('Time (s)'); %bf 粗體 /it 斜體 /rm 恢復正常
% ylabel('Torque (Nm)');
% zzz = legend(' Matlab',' Recurdyn');
% set(zzz,'box','off','FontSize',20,'FontName','Times New Roman','location','best')
% 
% set(gca,'Fontsize',20,'FontName','Times New Roman','FontWeight','Bold');
% 





plot_number = plot_number + 1; set(gca,'Fontsize',12,'FontName','Times New Roman','FontWeight','Bold');
figure(plot_number)
plot(t,omega_2(:,3));
figure_chienli2(300,260,0.23,0.25,0.7,0.6)
title('\omega_{2} v.s Time');
xlabel('Time (s)'); 
ylabel('\omega_{2} (rad/s)');





