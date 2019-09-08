function dy = overwatch(t,y)

dtr = pi/180;

global a m2 m6 m7 phi2 phi6 phi7 I2 I6 I7 L2 L6 L7


stiffness = 37000; damping = 11000;

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

% m2 = 53;
m3 = 7.1;
m4 = 10.59;
m5 = 47.93;
% m6 = 50.03;
% m7 = 60.07;

% L2 = 0.0247;
L3 = 0.1287;
L4 = 0.1467;
% L6 = 0.0749;
% L7 = 0.0277;

% phi2 = 0;
phi3 = 13*dtr;

phi4 = 0;
phi5 = 0;
% phi6 = 0;
% phi7 = 0;
pressure_angle = 20*dtr;

g = 9.81;


% I2 = 0.651;
I3 = 0.047;
I4 = 0.887;
% I6 = 0.313;
% I7 = 0.955;
IJ = 0.9;%kg-m*m
IF = 0;%kg-m*m
spring_initial = 0.27;  %m

theta_1 = (360-107.47)*dtr;
theta_2i = 0;
theta_3p = 27*dtr;
theta_7i = 287.45*dtr;
theta_5 = 270*dtr;

PG = [0 -0.75 0];
PH = [0.22 0.44 0];

PO2 = [0 0 0];
PO7 = [PO2(1,1) + r1*cos(theta_1) ,PO2(1,2) + r1*sin(theta_1) , 0];

dead_center = 0.139013027945341;



theta_2 = rem(y(1),2*pi);
omega_2 = [0,0,y(2)];
omega_7 = [0,0,-y(2)];


%% displacement

theta_7 = rem(theta_7i-(theta_2-0),2*pi);
A = 2*(r1*r6*cos(theta_1)-r2*r6*cos(theta_2)+r6*r7*cos(theta_7));
B = 2*(r1*r6*sin(theta_1)-r2*r6*sin(theta_2)+r6*r7*sin(theta_7));
C = r3^2-r1^2-r2^2-r6^2-r7^2+2*r1*r2*cos(theta_1-theta_2)-2*r1*r7*cos(theta_1-theta_7)+2*r2*r7*cos(theta_2-theta_7);

if  (2*atan2(B+((B*B+A*A-C*C)^0.5),A+C) + 360*dtr)/dtr < 360
    theta_6 = 2*atan2(B+((B*B+A*A-C*C)^0.5),A+C ) + 360*dtr;
else
    theta_6 = 2*atan2(B+((B*B+A*A-C*C)^0.5),A+C) -360*dtr;
end

theta_3 = -(acos((r1*cos(theta_1)-r2*cos(theta_2)+r6*cos(theta_6)+r7*cos(theta_7))/r3))+360*dtr;

D = 2*(-r1*cos(theta_1-theta_5)+r5p*cos(theta_5)-r6*cos(theta_5-theta_6)-r7*cos(theta_5-theta_7));
E = r1^2-r4^2+r5p^2+r6^2+r7^2-2*r1*r5p*cos(theta_1)+2*r1*r6*cos(theta_1-theta_6)+2*r1*r7*cos(theta_1-theta_7)...
    +2*r6*r7*cos(theta_6-theta_7)-2*r5p*r6*cos(theta_6)-2*r5p*r7*cos(theta_7);

r5 =(-D+sqrt(D^2-4*E))/2;

theta_4 = acos((-r1*cos(theta_1)+r5*cos(theta_5)+r5p-r6*cos(theta_6)-r7*cos(theta_7))/r4);

if 0 < theta_4 && theta_4 < 90*dtr
    theta_4 = -theta_4+360*dtr;
else
    theta_4 = 360*dtr-theta_4;
end




%% h

AA_ = [r3*sin(theta_3) -r6*sin(theta_6);-r3*cos(theta_3) r6*cos(theta_6)];
BB = [-r2*sin(theta_2)-r7*sin(theta_7);r2*cos(theta_2)+r7*cos(theta_7)];
X(1,:) = AA_\BB;
h3 = X(1);
h6 = X(2);
omega_3(3) = h3*omega_2(3);
omega_6(3) = h6*omega_2(3);


AA_ = [-r4*sin(theta_4) 0;r4*cos(theta_4) 1];
BB = [r6*h6*sin(theta_6)-r7*sin(theta_7);-r6*h6*cos(theta_6)+r7*cos(theta_7)];
X(1,:) = AA_\BB;
h4 = X(1);
hr5 = X(2);   %positive,when go down
omega_4(3) = h4*omega_2(3);
r5v(2) = -hr5*omega_2(3);


%% hp

AA_ = [r3*sin(theta_3) -r6*sin(theta_6);-r3*cos(theta_3) r6*cos(theta_6)];
BB = [h6^2*r6*cos(theta_6)+r7*cos(theta_7)-r2*cos(theta_2)-h3^2*r3*cos(theta_3);
    h6^2*r6*sin(theta_6)+r7*sin(theta_7)-r2*sin(theta_2)-h3^2*r3*sin(theta_3)];
X(1,:) = AA_\BB;
h3p = X(1);
h6p = X(2);
% alpha_3(3) = h3p*omega_2(3)^2;
% alpha_6(3) = h6p*omega_2(3)^2;

AA_ = [-r4*sin(theta_4) 0;r4*cos(theta_4) 1];
BB = [r4*h4^2*cos(theta_4)+r6*h6p*sin(theta_6)+r6*h6^2*cos(theta_6)+r7*cos(theta_7);
    r4*h4^2*sin(theta_4)-r6*h6p*cos(theta_6)+r6*h6^2*sin(theta_6)+r7*sin(theta_7)];
X(1,:) = AA_\BB;
h4p = X(1);
hr5p = X(2);

% alpha_4(3) = h4p*omega_2(3)^2;
% r5a(2) = -hr5p*omega_2(3)^2;


%% Link center mass displacement


    PA(1) = PO2(1) + r2*cos(theta_2);
    PA(2) = PO2(2) + r2*sin(theta_2);
    PA(3) = 0;
    
    PB(1) = PA(1) + r3*cos(theta_3);
    PB(2) = PA(2) + r3*sin(theta_3);
    PB(3) = 0;
    
    PD(1) = PO7(1) + r7*cos(theta_7);
    PD(2) = PO7(2) + r7*sin(theta_7);
    PD(3) = 0;
    
    PR5(1) = PO2(1);
    PR5(2) = PO2(2) - r5;
    PR5(3) = 0;
    
    PC(1) = PB(1) + r4*cos(theta_4);
    PC(2) = PB(2) + r4*sin(theta_4);
    PC(3) = 0;
    
    PS(1) = PA(1) + r3p*cos(theta_3 + theta_3p);
    PS(2) = PA(2) + r3p*sin(theta_3 + theta_3p);
    PS(3) = 0;





%% spring_length and damper velocity
rd = sqrt((PS(1,1)-PH(1,1))^2+(PS(1,2)-PH(1,2))^2);
spring_length = (PO2(1,2)-PG(1,2)) -r5;

theta_d = atan2(PH(1,2)-PS(2),PH(1,1)-PS(1));

AA_ = [cos(theta_d) -rd*sin(theta_d);sin(theta_d) rd*cos(theta_d)];
BB = [r2*sin(theta_2)+h3*r3p*sin(theta_3 + theta_3p);-r2*cos(theta_2)-h3*r3p*cos(theta_3 + theta_3p)];

X(1,:) = AA_\BB;

fd = X(1);
hd = X(2);


%% kinematic coefficient
f2x = -L2*sin(theta_2 + phi2);
f2y = L2*cos(theta_2  + phi2);

f2xp = -L2*cos(theta_2 + phi2);
f2yp = -L2*sin(theta_2 + phi2);

f3x = -r2*sin(theta_2)-L3*h3*sin(theta_3+phi3);
f3y = r2*cos(theta_2)+L3*h3*cos(theta_3+phi3);

f3xp = -r2*cos(theta_2)-L3*h3p*sin(theta_3+phi3)-L3*(h3^2)*cos(theta_3+phi3);
f3yp = -r2*sin(theta_2)+L3*h3p*cos(theta_3+phi3)-L3*(h3^2)*sin(theta_3+phi3);

f4x = -r2*sin(theta_2)-h3*r3*sin(theta_3)-L4*h4*sin(theta_4+phi4);
f4y = r2*cos(theta_2)+h3*r3*cos(theta_3)+L4*h4*cos(theta_4+phi4);

f4xp = -r2*cos(theta_2)-(h3^2)*r3*cos(theta_3)-h3p*r3*sin(theta_3)-(h4^2)*L4*cos(theta_4)-h4p*L4*sin(theta_4);
f4yp = -r2*sin(theta_2)-(h3^2)*r3*sin(theta_3)+h3p*r3*cos(theta_3)-(h4^2)*L4*sin(theta_4)+h4p*L4*cos(theta_4);

f5y = -hr5;        %f5y , positive when go up
f5yp= -hr5p;

f6x = r7*sin(theta_7)-h6*L6*sin(theta_6 + phi6);
f6y = -r7*cos(theta_7)+h6*L6*cos(theta_6 + phi6);

f6xp = -r7*cos(theta_7)-(h6^2)*L6*cos(theta_6 + phi6)-h6p*L6*sin(theta_6 + phi6);
f6yp = -r7*sin(theta_7)-(h6^2)*L6*sin(theta_6 + phi6)+h6p*L6*cos(theta_6 + phi6);

f7x = L7*sin(theta_7 + phi7);
f7y = -L7*cos(theta_7 + phi7);

f7xp = -L7*cos(theta_7 + phi7);
f7yp = -L7*sin(theta_7 + phi7);


LDCP = [0,PG(2)+dead_center,0];
% 
% for j = 1:3
%     vector_O2_LDCP(j) = LDCP(j)-PO2(j);
%     vector_O7_LDCP(j) = LDCP(j)-PO7(j);
%     vector_H_LDCP(j) = LDCP(j)-PH(j);
% end

time = 8;

if t > time    
    
    if spring_length - dead_center <= 0.005 && r5v(2) < 0
        Fw = 89000;
% Fw = 0;
    else
        Fw = 0;
    end
    
else
    Fw = 0;
end


sigmaA = m2*(f2x^2+f2y^2) + I2*1 + m3*(f3x^2+f3y^2) + I3*h3^2 + m4*(f4x^2+f4y^2)+I4*h4^2 ...
         + m6*(f6x^2+f6y^2)+ I6*h6^2 + m7*(f7x^2+f7y^2) + I7*(-1)^2 + m5*(f5y^2);
sigmaB = m2*(f2x*f2xp+f2y*f2yp) + m3*(f3x*f3xp+f3y*f3yp)+I3*h3*h3p + m4*(f4x*f4xp+f4y*f4yp)+I4*h4*h4p ...
         + m6*(f6x*f6xp+f6y*f6yp)+I6*h6*h6p + m7*(f7x*f7xp+f7y*f7yp) + m5*(f5y*f5yp);
sigmaUg = -(m2*f2y+m3*f3y+m4*f4y+m5*f5y+m6*f6y+m7*f7y)*g;
Uspring = +stiffness*(spring_initial-spring_length)*f5y ;
Udamping = -damping*fd^2; 
UFw = Fw*f5y;



Motor_rpm = omega_2(3)*gr*sr*60/(2*pi);

if Motor_rpm >= 0 && Motor_rpm <=1500
    Torque = a*(27000-7*Motor_rpm)*gr*sr;

elseif Motor_rpm >= 1500 && Motor_rpm <=1800
    Torque = a*(16500-55*(Motor_rpm-1500))*gr*sr;
    
else Torque = 0;


end


dy(1) = y(2);
dy(2) = (Torque - sigmaB*y(2)^2 + sigmaUg + Uspring + Udamping*y(2) + UFw)/(sigmaA);
dy = [dy(1); dy(2)];