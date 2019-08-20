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

I2 = 0.651;
I3 = 0.047;
I4 = 0.887;
I6 = 0.313;
I7 = 0.955;
%%
theta_1 = -107.47*dtr;
theta_2i = 0;
theta_3p = 27*dtr;
theta_7i = 287.45*dtr;
theta_5 = 270*dtr;

theta_2 = zeros(361,1);
theta_7 = zeros(361,1);


theta_6 = zeros(361,1);
theta_3 = zeros(361,1);
r5 = zeros(361,1);
theta_4 = zeros(361,1);
theta_d = zeros(361,1);

h3 = zeros(361,1);
hr5 = zeros(361,1);
h4 = zeros(361,1);
h6 = zeros(361,1);

h3p = zeros(361,1);
hr5p = zeros(361,1);
h4p = zeros(361,1);
h6p = zeros(361,1);


%% Displacement
interval = 1*dtr;
A = 2*(r1*r6*cos(theta_1)-r2*r6*cos(theta_2)+r6*r7*cos(theta_7));
B = 2*(r1*r6*sin(theta_1)-r2*r6*sin(theta_2)+r6*r7*sin(theta_7));
C = r3^2-r1^2-r2^2-r6^2-r7^2+2*r1*r2*cos(theta_1-theta_2)-2*r1*r7*cos(theta_1-theta_7)+2*r2*r7*cos(theta_2-theta_7);

for i = 1:361
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


%% Angular velocity
omega = [10 10 10];
for k = 1:1

omega_2 = [0 0 omega(k)];
omega_7 = [0 0 -omega(k)];
omega_6 = zeros(361,3);
omega_3 = zeros(361,3);
omega_4 = zeros(361,3);
r5v =  zeros(361,3);
X = zeros(361,2);

for i = 1:361
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

alpha_6 = zeros(361,3);
alpha_3 = zeros(361,3);
alpha_4 = zeros(361,3);
r5a = zeros(361,3);

for i = 1:361
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

PA = zeros(361,3);
PB = zeros(361,3);
PD = zeros(361,3);
PR5 = zeros(361,3);
PC = zeros(361,3);
PS = zeros(361,3);    %the joint between damper and r3` 


for i = 1:361
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

vector_O2_A = zeros(361,3); 
vector_A_B = zeros(361,3); 
vector_B_C = zeros(361,3); 
vector_O7_D = zeros(361,3); 
vector_D_B = zeros(361,3); 

for i = 1:361
    for j = 1:2
vector_O2_A(i,j) = PA(i,j)-PO2(1,j);
vector_A_B(i,j) = PB(i,j)-PA(i,j);
vector_B_C(i,j) = PC(i,j)-PB(i,j);
vector_O7_D(i,j) = PD(i,j)-PO7(1,j);
vector_D_B(i,j) = PB(i,j)-PD(i,j);

    end
end

%% Linkage mass-center displacement

PG_r2 = zeros(361,3);
PG_r3 = zeros(361,3);
PG_r4 = zeros(361,3);
PG_r5 = zeros(361,3);
PG_r6 = zeros(361,3);
PG_r7 = zeros(361,3);


vector_O2_Gr2 = zeros(361,3);
vector_A_Gr3 = zeros(361,3);
vector_B_Gr4 = zeros(361,3);
vector_O7_Gr7 = zeros(361,3);
vector_D_Gr6 = zeros(361,3);
vector_B_C = zeros(361,3);

for i = 1 :361
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
VA = zeros(361,3);
VB = zeros(361,3);
VD = zeros(361,3);


for i = 1:361

VA(i,:) = cross(omega_2(1,:),vector_O2_A(i,:));
VB(i,:) = cross(omega_3(i,:),vector_A_B(i,:)) + VA(i,:);
VD(i,:) = cross(omega_7(1,:),vector_O7_D(i,:));

end

%%  Linkage mass velocity
Vr2 = zeros(361,3);
Vr3 = zeros(361,3);
Vr4 = zeros(361,3);
Vr6 = zeros(361,3);
Vr7 = zeros(361,3);


for i = 1:361
Vr2(i,:) = cross(omega_2(1,:),vector_O2_Gr2(i,:));
Vr3(i,:) = cross(omega_3(i,:),vector_A_Gr3(i,:))+VA(i,:);
Vr4(i,:) = cross(omega_4(i,:),vector_B_Gr4(i,:))+VB(i,:);
Vr6(i,:) = cross(omega_6(i,:),vector_D_Gr6(i,:))+VD(i,:);
Vr7(i,:) = cross(omega_7(1,:),vector_O7_Gr7(i,:));
end

%%  Linkage acceleration

AA = zeros(361,3);
AB = zeros(361,3);
AC = zeros(361,3);
AD = zeros(361,3);

for i = 1:361
AA(i,:) = cross(omega_2(1,:),VA(i,:));
AB(i,:) = AA(i,:) + cross(omega_3(i,:),cross(omega_3(i,:),vector_A_B(i,:))) + cross(alpha_3(i,:),vector_A_B(i,:));
AC(i,:) = r5a(i,:);
AD(i,:) = cross(omega_7(1,:),VD(i,:));
end

%% Linkage mass acceleration
Ar2 = zeros(361,3);
Ar3 = zeros(361,3);
Ar4 = zeros(361,3);
Ar6 = zeros(361,3);
Ar7 = zeros(361,3);

for i = 1:361
Ar2(i,:) = cross(omega_2(1,:),Vr2(i,:));
Ar3(i,:) = AA(i,:) + cross(omega_3(i,:),cross(omega_3(i,:),vector_A_Gr3(i,:))) + cross(alpha_3(i,:),vector_A_Gr3(i,:));

Ar4(i,:) = AB(i,:) + cross(omega_4(i,:),cross(omega_4(i,:),vector_B_Gr4(i,:))) + cross(alpha_4(i,:),vector_B_Gr4(i,:));

Ar6(i,:) = AD(i,:) + cross(omega_6(i,:),cross(omega_6(i,:),vector_D_Gr6(i,:))) + cross(alpha_6(i,:),vector_D_Gr6(i,:));

Ar7(i,:) = cross(omega_7(1,:),Vr7(i,:));

end



%% spring_length and damper velocity
spring_length = zeros(361,1);
rd =  zeros(361,1);
rdv = zeros(361,3);

for i = 1:361
    rd(i) = sqrt((PS(i,1)-PH(1,1))^2+(PS(i,2)-PH(1,2))^2);
    spring_length(i) = (PO2(1,2)-PG(1,2)) -r5(i);
    theta_d(i) = atan2(PH(1,2)-PS(i,2),PH(1,1)-PS(i,1));
    
    AA = [cos(theta_d(i)) -rd(i)*sin(theta_d(i));sin(theta_d(i)) rd(i)*cos(theta_d(i))];
    BB = [omega_2(1,3)*r2*sin(theta_2(i))+r3p*omega_3(i,3)*sin(theta_3(i) + theta_3p);-omega_2(1,3)*r2*cos(theta_2(i))-r3p*omega_3(i,3)*cos(theta_3(i) + theta_3p)];
    
    X(i,:) = AA\BB;
    
    
    rdv(i,2) = X(i,1);

    
end

% i = 1:1:361;
% x = theta_2(i)/dtr;
% y = rdv(i,2);
% 
% 
% plot(x,y,'Linewidth',2);
% 
% hold on;
% title('Damper Velocity vs \theta_2','FontSize',20,'FontName','Times New Roman');
% xlabel('\theta_2 (\circ)','FontSize',20,'FontName','Times New Roman'); %bf 粗體 /it 斜體 /rm 恢復正常
% ylabel('Velocity (m/s)','FontSize',20,'FontName','Times New Roman');
% set(gca,'Fontsize',10,'FontName','Times New Roman');
% set(gcf,'Position',[900 150 700 700]);


%% load file


% file = load('E:\機器動力學\作業一\Report_plot\spring_length.txt');
% data_1 = file(:,1);
% data_2 = file(:,2);
% 
% figure_chienli2(300,270,0.23,0.25,0.7,0.6)
%  
% i = 1:1:360;
% 
% x = theta_2(i)/dtr;
% y = spring_length(i);
% xlim([0,360]);
% 
% title('{Spring Length} vs \theta_2');
% xlabel('\theta_2 (\circ)'); %bf 粗體 /it 斜體 /rm 恢復正常
% ylabel('Length (mm)');
% set(gca,'Fontsize',12,'FontName','Times New Roman','Xtick',60*(0:6),'FontWeight','Bold');
% % set(gcf,'Position',[900 150 700 700]);
% 
% 
% 
% plot(x,y*1000,'Linewidth',2);
% hold on
% 
% plot(data_1,data_2+0.150192030536832*1000,'--r','Linewidth',2);
% zzz = legend(' Matlab',' Recurdyn');
% set(zzz,'box','off','FontSize',12,'FontName','Times New Roman','location','best')
% 


% title('\itv_{5} vs \theta_2');
% xlabel('\theta_2 (\circ)');
% ylabel('\omega_4 (rad/s)');









%% Mass-center track
% 
% 
% 
% i = 1:1:361;
% x = PG_r7(i,1);
% y = PG_r7(i,2);
% 
% figure
% plot(1000*x,1000*y);
% 
% 
% title('Mass Center - Link2','FontSize',20,'FontName','Times New Roman');
% xlabel('x (mm)','FontSize',20,'FontName','Times New Roman'); %bf 粗體 /it 斜體 /rm 恢復正常
% ylabel('y (mm)','FontSize',20,'FontName','Times New Roman');
% %set(gca,'Fontsize',10,'FontName','Times New Roman');
% set(gcf,'Position',[900 150 700 700]);

%% Mass-center track x ,y 
 
% i = 1:1:361;
% x = theta_2(i)/dtr;
% y1 = 1000*PG_r5(i,1);
% y2 = 1000*PG_r5(i,2);
% plot(x,y1,x,y2,'Linewidth',2);
% 
% title('Mass Center - Link5','FontSize',20,'FontName','Times New Roman');
% xlabel('\theta_2(\circ)','FontSize',24,'FontName','Times New Roman'); %bf 粗體 /it 斜體 /rm 恢復正常
% ylabel('Position (mm)','FontSize',24,'FontName','Times New Roman');
% zzz = legend(' x position',' y position');
% set(zzz,'box','off','FontSize',16,'FontName','Times New Roman','location','northeast')
% 
% 
% 

%% h v.s theta_2

% i = 1:1:360;
% x = theta_2(i)/dtr;
% y = hr5(i);
% 
% figure
% plot(x,y);
% 
% xlim([0,360]);
% 
% title('f_{5} vs \theta_2','FontSize',20,'FontName','Times New Roman');
% xlabel('\theta_2(\circ)','FontSize',20,'FontName','Times New Roman'); %bf 粗體 /it 斜體 /rm 恢復正常
% ylabel('f_{5}','FontSize',20,'FontName','Times New Roman');
% set(gca,'Fontsize',20,'FontName','Times New Roman');
% %set(gcf,'Position',[900 150 700 700]);

%% hp v.s theta_2
%  
% i = 1:1:361;
% x = theta_2(i)/dtr;
% y = hr5p(i);
% 
% 
% plot(x,y,'Linewidth',2);
% 
% title('f_{5} vs \theta_2','FontName','Times New Roman');
% xlabel('\theta_2 (\circ)','FontName','Times New Roman');
% ylabel('f_{5}','FontName','Times New Roman');
% 

%% omega v.s theta_2
 
% i = 1:1:360;
% x = theta_2(i)/dtr;
% y = omega_3(i,3);
% 
% 
% 
% 
% plot(x,y,'Linewidth',2);
% hold on
% title('\itv_{5} vs \theta_2');
% xlabel('\theta_2 (\circ)');
% ylabel('\omega_4 (rad/s)');


%% r5v v.s theta_2
%  
% i = 1:1:360;
% x = theta_2(i)/dtr;
% y = rd(i);
% 
% figure
% plot(x,y);
% 
% title('\omega_3 vs \theta_2');
% xlabel('\theta_2 (\circ)');
% ylabel('\omega_3 (rad/s)');
% 


%% Mass-center velocity v.s theta_2

% i = 1:1:360;
% x = theta_2(i)/dtr;
% y = sqrt(Vr6(i,1).^2+Vr6(i,2).^2);
% 
% hold on;
% plot(x,y,'Linewidth',2);
% 
% title('\itV_{g6} vs \theta_2');
% xlabel('\theta_2 (\circ)');
% ylabel('\itV_{g6} (m/s)');

%% Alpha_X v.s theta_2

% i = 1:1:361;
% x = theta_2(i)/dtr;
% y = r5a(i,2);
% 
% plot(x,y,'Linewidth',2);
% 
% hold on;
% title('a_{5} vs \theta_2');
% xlabel('\theta_2 (angle)');
% ylabel('a_{5} (m/s^2)');
% 


%% Mass-center acceleration v.s theta_2
% 
% i = 1:1:361;
% x = theta_2(i)/dtr;
% y = sqrt(Ar4(i,1).^2+Ar4(i,2).^2);
% 
% 
% plot(x,y,'Linewidth',2);
% 
% hold on;
% 
% title('\itA_{g4} vs \theta_2');
% xlabel('\theta_2 (\circ)');
% ylabel('\itA_{g4} (m/s^2)');





    %box on
    %set(gca,'Fontsize',24,'FontName','Times New Roman');
    %zzz = legend('positive','negative','all');
    %set(zzz,'box','off','FontWeight','bold','FontSize',16,'FontName','Times New Roman','orientation','horizontal','location','northwest')
    %xlim([0 360]);
    %ylim([22 33]);
    %set(gcf,'Position',[900 150 680 680*0.8]);
    %set(gca,'Position',[0.18 0.18 0.7 0.7]);
    %xlabel('Cam Angle(\circ)','FontSize',24,'FontName','Times New Roman'); %bf 粗體 /it 斜體 /rm 恢復正常
    %ylabel('s(mm)','FontSize',24,'FontName','Times New Roman');
   

    box on
%     set(gca,'Fontsize',24,'FontName','Times New Roman');
    %zzz = legend('positive','negative','all');
    %set(zzz,'box','off','FontWeight','bold','FontSize',16,'FontName','Times New Roman','orientation','horizontal','location','northwest')
    xlim([0 360]);
  
    %set(gcf,'Position',[900 150 680 680*0.8]);
    %set(gca,'Position',[0.18 0.18 0.7 0.7]);
    %xlabel('\theta_2(\circ)','FontSize',24,'FontName','Times New Roman'); %bf 粗體 /it 斜體 /rm 恢復正常
    %ylabel('\itv_{5} (m/s)','FontSize',24,'FontName','Times New Roman');
    %zzz = legend(' \omega = 10 rad/s',' \omega = 60 rad/s',' \omega = 120 rad/s');
    %set(zzz,'box','off','FontSize',16,'FontName','Times New Roman','location','northeast')

end


%% Graphic simulation

for i = 1:361
   clf
   
   plot(PO7(1,1),PO7(1,2),'s','MarkerSize',8); hold on %o2
   plot(PO2(1,1),PO2(1,2),'o','MarkerSize',8); %B

  
   plot([PO7(1,1),PO2(1,1)],[PO7(1,2),PO2(1,2)],'LineWidth',2); %R2
   
   plot(PA(i,1),PA(i,2),'o','MarkerSize',8); %C
    
   plot([PO2(1,1),PA(i,1)],[PO2(1,2),PA(i,2)],'LineWidth',2); %R3
   
   plot(PS(i,1),PS(i,2),'o','MarkerSize',8); %C
   
   plot([PA(i,1),PS(i,1)],[PA(i,2),PS(i,2)],'LineWidth',2); %R3
   
   plot(PB(i,1),PB(i,2),'o','MarkerSize',8); %D
   plot([PA(i,1),PB(i,1)],[PA(i,2),PB(i,2)],'LineWidth',2); %R3'
   
   plot(PD(i,1),PD(i,2),'o','MarkerSize',8); %E
   plot([PO7(1,1),PD(i,1)],[PO7(1,2),PD(i,2)],'LineWidth',2); %R5
   
   plot([PD(i,1),PB(i,1)],[PD(i,2),PB(i,2)],'LineWidth',2); %R5
   
   plot(PR5(i,1),PR5(i,2),'o','MarkerSize',8); %F

   
   plot(PC(i,1),PC(i,2),'s','MarkerSize',8); %o4
   plot([PR5(i,1),PC(i,1)],[PR5(i,2),PC(i,2)],'LineWidth',2); %R3'
   
   plot([PB(i,1),PC(i,1)],[PB(i,2),PC(i,2)],'LineWidth',2); %R3'
  
   plot(PG_r2(i,1),PG_r2(i,2),'o','MarkerSize',8); 
   plot(PG_r3(i,1),PG_r3(i,2),'o','MarkerSize',8); 
   plot(PG_r4(i,1),PG_r4(i,2),'o','MarkerSize',8); 
   plot(PG_r5(i,1),PG_r5(i,2),'o','MarkerSize',8); 
   plot(PG_r6(i,1),PG_r6(i,2),'o','MarkerSize',8); 
   plot(PG_r7(i,1),PG_r7(i,2),'o','MarkerSize',8);
   
   
   hold off
   axis equal
   xlim([-0.2 0.2]); ylim([-0.7 0.06]); 
   grid on
   pause(0.001)
end


