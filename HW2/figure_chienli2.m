function [] = figure_chienli2(horz,vert,l,d,r,u)
figure
ss = get(0,'screensize'); %The screen size
width = ss(3);
height = ss(4);
% horz = 335; %for two figures in A4 page
% horz = 230; %for three figures in A4 page
% horz = 700; %for two figures in A4 page
set(gcf,'Position',[width-horz-5, height-(vert+80), horz, vert])
% set(gcf,'Position',[width-horz-5, height-(vert+80), horz, vert])
% set(gca,'Position',[(l/horz),(d/vert),(horz-r)/horz,(vert-u)/vert])
set(gca,'Position',[l,d,r,u])
set(gca,'FontWeight','bold','FontSize',10,'FontName','Times New Roman')
hold all
box on
fontname= 'Times New Roman';
set(0,'defaultaxesfontname',fontname);
set(0,'defaulttextfontname',fontname);
fontsize= 10;
set(0,'defaultaxesfontsize',fontsize);
set(0,'defaulttextfontsize',fontsize);
fontweight='bold';
set(0,'defaultaxesfontweight',fontweight);
set(0,'defaulttextfontweight',fontweight);
% xlabel('Time (sec)')
% ylabel('D_y (mm)')
% title('Yaw motor motion')
% legend('experiment','trajectory')
% set(gca,'XTick',[0:0.5:4])
% set(gca,'XTickLabel',['0';' ';'1';' ';'2';' ';'3';' ';'4'])
% xlim([time_part(1),time_part(end)])
% ylim([min([data1_part;data2_part])-10,max([data1_part;data2_part])+10])
% xtick=roundn((time_part(end)-time_part(1))/4,-1);
% xticki=(fix(time_part(1)/xtick)+1)*xtick;
% xticke=fix(time_part(end)/xtick)*xtick;
% set(gca,'XTick',[xticki:xtick:xticke])
end

