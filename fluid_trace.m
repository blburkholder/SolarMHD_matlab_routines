%   for i = 5
%   prepare3(i)
%    mentrop_prepare_fluelt(i)
%   end
l = 0;
  load(['mhd',num2str(l),'.mat']); load(['flt_mentrop',num2str(l),'.mat']);
% wy = [-flipud(y(4:end));y];
% wbx = [-fliplr(flipud(bx(3:end,:,:)));bx(2:end,:,:)];
% wby = [-fliplr(flipud(by(3:end,:,:)));by(2:end,:,:)];
% wbz = [fliplr(flipud(bz(3:end,:,:)));bz(2:end,:,:)];
% wsx = [-fliplr(flipud(sx(3:end,:,2)));sx(2:end,:,2)];
% wsy = [-fliplr(flipud(sy(3:end,:,2)));sy(2:end,:,2)];
% wrho = [fliplr(flipud(rho(3:end,:,2)));rho(2:end,:,2)];
vx = sx./rho;
vy = sy./rho;

%use line symmetric boundaries to get every footpoint in 1st quadrant
[xend1_t0,yend1_t0] = fix_coords(x,y,xfin1,yfin1);
zend1_t0 = zfin1;
[xend2_t0,yend2_t0] = fix_coords(x,y,xfin2,yfin2);
zend2_t0 = zfin2;
%keep starting points of field lines at top boundary
x_fluelt_t0 = xpo;
y_fluelt_t0 = ypo;
z_fluelt_t0 = zpo;

xend_open_tf = zeros(size(xfin1));
yend_open_tf = zeros(size(xfin1));
xend_conv_tf = zeros(size(xfin1));
yend_conv_tf = zeros(size(xfin1));
xend_t0 = zeros(size(xfin1));
yend_t0 = zeros(size(xfin1));
time1 = time;
  
  l=6;
  load(['mhd',num2str(l),'.mat']); load(['flt_mentrop',num2str(l),'.mat']);

  %take initial location of both ends of the field line (whether on p-sphere or not)
  %and push with photospheric convection  
  [xend1_conv_tf,yend1_conv_tf] = velocity_perturbation(x,y,xend1_t0,yend1_t0,time-time1);
  [xend2_conv_tf,yend2_conv_tf] = velocity_perturbation(x,y,xend2_t0,yend2_t0,time-time1);
  %do this because its fucking necessary
  xend1_conv_tf = xend1_conv_tf';
  xend2_conv_tf = xend2_conv_tf';
  yend1_conv_tf = yend1_conv_tf';
  yend2_conv_tf = yend2_conv_tf';
  %get location of actual photospheric footpoints at end of simulation
  [xend1_open_tf,yend1_open_tf] = fix_coords(x,y,xfin1,yfin1);
  zend1_tf = zfin1;
  [xend2_open_tf,yend2_open_tf] = fix_coords(x,y,xfin2,yfin2);
  zend2_tf = zfin2;

  x_fluelt_tf = xpo;
  y_fluelt_tf = ypo;
  z_fluelt_tf = zpo;

    for i = 1:length(xend_t0(:))    
        if zend1_t0(i) < 1 && zend1_tf(i) < 1
            xend_t0(i) = xend1_t0(i); yend_t0(i) = yend1_t0(i);
            xend_open_tf(i) = xend1_open_tf(i); yend_open_tf(i) = yend1_open_tf(i);
            xend_conv_tf(i) = xend1_conv_tf(i); yend_conv_tf(i) = yend1_conv_tf(i); 
        elseif zend2_t0(i) < 1 && zend2_tf(i) < 1
            xend_t0(i) = xend2_t0(i); yend_t0(i) = yend2_t0(i);
            xend_open_tf(i) = xend2_open_tf(i); yend_open_tf(i) = yend2_open_tf(i);
            xend_conv_tf(i) = xend2_conv_tf(i); yend_conv_tf(i) = yend2_conv_tf(i); 
        %elseif zend1_t0(i) < 1 && zend2_t0(i) < 1
        else
%             if zend2_tf(i) < 1
%                 xend_t0(i) = xend2_t0(i); yend_t0(i) = yend2_t0(i);
%                 xend_open_tf(i) = xend2_open_tf(i); yend_open_tf(i) = yend2_open_tf(i);
%                 xend_conv_tf(i) = xend2_conv_tf(i); yend_conv_tf(i) = yend2_conv_tf(i); 
%             elseif zend1_tf(i) < 1
%                 xend_t0(i) = xend1_t0(i); yend_t0(i) = yend1_t0(i);
%                 xend_open_tf(i) = xend1_open_tf(i); yend_open_tf(i) = yend1_open_tf(i);
%                 xend_conv_tf(i) = xend1_conv_tf(i); yend_conv_tf(i) = yend1_conv_tf(i); 
%             end
%         else
            xend_t0(i) = xend_t0(i-1); yend_t0(i) = yend_t0(i-1);
            xend_open_tf(i) = xend_open_tf(i-1); yend_open_tf(i) = yend_open_tf(i-1);
            xend_conv_tf(i) = xend_conv_tf(i-1); yend_conv_tf(i) = yend_conv_tf(i-1); 
%              xend_t0(i) = nan; yend_t0(i) = nan;
%              xend_open_tf(i) = nan; yend_open_tf(i) = nan;
%              xend_conv_tf(i) = nan; yend_conv_tf(i) = nan; 
            %xend_t0(i) = 0; yend_t0(i) = 0;
            %xend_open_tf(i) = 0; yend_open_tf(i) = 0;
            %xend_conv_tf(i) = 0; yend_conv_tf(i) = 0; 
        end
    end

% rect_flueltst0_cy = [-fliplr(flipud(y_fluelt_t0(2:end,:)));y_fluelt_t0];
% rect_flueltst0_cx = [x(end-1)-fliplr(flipud(x_fluelt_t0(2:end,:)));x_fluelt_t0];
% 
% rect_flueltstf_cy = [-fliplr(flipud(y_fluelt_tf(2:end,:)));y_fluelt_tf];
% rect_flueltstf_cx = [x(end-1)-fliplr(flipud(x_fluelt_tf(2:end,:)));x_fluelt_tf];
% 
% rect_xend_t0 = [fliplr(flipud(x(end-1) - xend_t0(2:end,:)));xend_t0];
% rect_yend_t0 = [fliplr(flipud(-yend_t0(2:end,:)));yend_t0];
% rect_xend_tf = [fliplr(flipud(x(end-1) - xend_open_tf(2:end,:)));xend_open_tf];
% rect_yend_tf = [fliplr(flipud(-yend_open_tf(2:end,:)));yend_open_tf];
% rect_xend_tflu = [fliplr(flipud(x(end-1) - xend_conv_tf(2:end,:)));xend_conv_tf];
% rect_yend_tflu = [fliplr(flipud(-yend_conv_tf(2:end,:)));yend_conv_tf];

tubetag1 = zeros(size(y_fluelt_t0));
tubetag2 = zeros(size(y_fluelt_t0));

ti_tube1 = (xend_t0 > 20 & xend_t0 < 60 & yend_t0 > 40 & yend_t0 < 70) | xend_t0 == 0;
ti_tube2 = yend_t0 > 70;
ti_tube3 =  (yend_t0 > 0 & yend_t0 < 34) | xend_t0 > 80; 

tf_tube1 = (xend_open_tf > 20 & xend_open_tf < 60 & yend_open_tf > 40 & yend_open_tf < 70) | xend_open_tf == 0;
tf_tube2 = yend_open_tf > 70;
tf_tube3 = (yend_t0 > 0 & yend_t0 < 34) | xend_open_tf > 80;

tubetag1(ti_tube1) = 1;
tubetag1(ti_tube2) = 2;
tubetag1(ti_tube3) = 3;

tubetag2(tf_tube1) = 1;
tubetag2(tf_tube2) = 2;
tubetag2(tf_tube3) = 3;

ti_xtube1 = x_fluelt_t0(ti_tube1); ti_ytube1 = y_fluelt_t0(ti_tube1);
ti_xtube2 = x_fluelt_t0(ti_tube2); ti_ytube2 = y_fluelt_t0(ti_tube2);
ti_xtube3 = x_fluelt_t0(ti_tube3); ti_ytube3 = y_fluelt_t0(ti_tube3);

tf_xtube1 = x_fluelt_tf(tf_tube1); tf_ytube1 = y_fluelt_tf(tf_tube1);
tf_xtube2 = x_fluelt_tf(tf_tube2); tf_ytube2 = y_fluelt_tf(tf_tube2);
tf_xtube3 = x_fluelt_tf(tf_tube3); tf_ytube3 = y_fluelt_tf(tf_tube3);

figure; hold on
scatter(xend_t0(ti_tube1),yend_t0(ti_tube1),'r.')
scatter(xend_t0(ti_tube2),yend_t0(ti_tube2),'g.')
scatter(xend_t0(ti_tube3),yend_t0(ti_tube3),'b.')
%scatter(xend_t0(:),yend_t0(:),'r.')
daspect([1 1 1])
axis([x(1) x(end) y(1) y(end)])
xlabel('x')
ylabel('y')
title('open flux at photosphere')

figure; hold on
scatter(ti_xtube1,ti_ytube1,'r.')
scatter(ti_xtube2,ti_ytube2,'g.')
scatter(ti_xtube3,ti_ytube3,'b.')
daspect([1 1 1])
axis([x(1) x(end) y(1) y(end)])
xlabel('x')
ylabel('y')
title('open flux tubes initial equilibrium')

figure; hold on
scatter(tf_xtube1,tf_ytube1,'r.')
scatter(tf_xtube2,tf_ytube2,'g.')
scatter(tf_xtube3,tf_ytube3,'b.')
flipx = x_fluelt_tf(tubetag1 ~= tubetag2 & tubetag1 > 0 & y_fluelt_tf > 50);
flipy = y_fluelt_tf(tubetag1 ~= tubetag2 & tubetag1 > 0 & y_fluelt_tf > 50);
scatter(flipx,flipy,'y.','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
daspect([1 1 1])
axis([x(1) x(end) y(1) y(end)])
xlabel('x')
ylabel('y')
title('reconnection between flux tubes')

figure; hold on
scatter(tf_xtube1,tf_ytube1,'r.')
scatter(tf_xtube2,tf_ytube2,'g.')
scatter(tf_xtube3,tf_ytube3,'b.')
daspect([1 1 1])
axis([x(1) x(end) y(1) y(end)])
xlabel('x')
ylabel('y')
title(['open flux tubes t =',num2str(time)])

bzall = interp3(x,y,z,bz,xpo,ypo,60*ones(size(xpo)));
bzflip = interp3(x,y,z,bz,flipx,flipy,60*ones(size(flipx)));
[net,tot] = flux(bzall,x(2)-x(1),y(2)-y(1))
sum(abs(bzflip*(x(2)-x(1))*(y(2)-y(1))))/tot

% load('mhd0.mat');
% figure; hold on
% pcolor(x,y,sqrt(vx(:,:,1).^2+vy(:,:,1).^2))
% colormap(gray)
% colorbar
% shading interp
% scatter(xend_t0(ti_tube1),yend_t0(ti_tube1),'r.')
% scatter(xend_t0(ti_tube2),yend_t0(ti_tube2),'g.')
% scatter(xend_t0(ti_tube3),yend_t0(ti_tube3),'b.')
% [xxx,yyy] = meshgrid(x,y);
% h = streamslice(xxx,yyy,vx(:,:,1),vy(:,:,1));
% set(h,'Color','k');
% reso = 5000;
% h = streamline(stream3(x,y,z,-bx,-by,-bz,ti_xtube3(1:reso:end),ti_ytube3(1:reso:end),60*ones(size(ti_ytube3(1:reso:end)))));
% set(h,'Color','b')
% h = streamline(stream3(x,y,z,bx,by,bz,ti_xtube2(1:reso:end),ti_ytube2(1:reso:end),60*ones(size(ti_ytube2(1:reso:end)))));
% set(h,'Color','g')
% h = streamline(stream3(x,y,z,bx,by,bz,ti_xtube1(1:reso:end),ti_ytube1(1:reso:end),60*ones(size(ti_ytube1(1:reso:end)))));
% set(h,'Color','r')
% scatter3(ti_xtube1,ti_ytube1,61*ones(size(ti_xtube1)),'r.','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75)
% scatter3(ti_xtube2,ti_ytube2,61*ones(size(ti_xtube2)),'g.','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75)
% scatter3(ti_xtube3,ti_ytube3,61*ones(size(ti_xtube3)),'b.','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75)
% axis([0 x(end) 0 y(end)])
% daspect([1 1 1])
% xlabel('x'); ylabel('y'); zlabel('z')
% view([45 30])
%this adds some field lines from the very thin layer of open fux
% [xxx,yyy] = meshgrid(xpo,ypo);
% qsl = zfin1 > 50 & yyy > 7.5 & yyy < 17.5;
% xx = xxx(qsl); yy = yyy(qsl);
% streamline(x,y,z,bx,by,bz,xx(1:35:end),yy(1:35:end),zeros(size(xx(1:35:end))));

% bx_top = interp3(x,y,z,bx,cxt,cyt,czt);
% by_top = interp3(x,y,z,by,cxt,cyt,czt);
% bz_top = interp3(x,y,z,bz,cxt,cyt,czt);
% open = abs(bz_top) > abs(by_top) & abs(bz_top) > abs(bx_top);

load('mhd6.mat'); figure
cond1 = yend_t0 > 55 & tubetag1 ~= tubetag2 & xend_t0 < 60 & xend_t0 > 20;
scatter(xend_t0(yend_t0 > 70 & cond1),yend_t0(cond1 & yend_t0 > 70),'g.');
hold on
scatter(xend_t0(cond1 & yend_t0 < 70),yend_t0(cond1 & yend_t0 < 70),'r.');
scatter(xend_conv_tf(cond1 & yend_conv_tf > 70),yend_conv_tf(cond1 & yend_conv_tf > 70),'g.');
scatter(xend_conv_tf(cond1 & yend_conv_tf < 70),yend_conv_tf(cond1 & yend_conv_tf < 70),'r.');
scatter(xend_open_tf(cond1),yend_open_tf(cond1),'k.');
reso = 200;
flines_x = x_fluelt_tf(cond1 & yend_t0 > 70); flines_y = y_fluelt_tf(cond1 & yend_t0 > 70);
flines_z = 60*ones(size(flines_y));
h1 = streamline(stream3(x,y,z,bx,by,bz,flines_x(1:reso:end),flines_y(1:reso:end),flines_z(1:reso:end)));
set(h1,'Color','g')
gdtx = x_fluelt_tf(cond1 & yend_open_tf > 70); gdty = y_fluelt_tf(cond1 & yend_open_tf > 70);
gdtz = 61*ones(size(x_fluelt_tf(cond1 & yend_open_tf > 70)));
s = scatter3(gdtx,gdty,gdtz,'k','filled'); s.MarkerFaceAlpha = 0.5;
flines_x = x_fluelt_tf(cond1 & yend_t0 < 70); flines_y = y_fluelt_tf(cond1 & yend_t0 < 70);
flines_z = 60*ones(size(flines_y));
h1 = streamline(stream3(x,y,z,bx,by,bz,flines_x(1:reso:end),flines_y(1:reso:end),flines_z(1:reso:end)));
set(h1,'Color','r')
gdtx = x_fluelt_tf(cond1 & yend_open_tf < 70); gdty = y_fluelt_tf(cond1 & yend_open_tf < 70);
gdtz = 61*ones(size(x_fluelt_tf(cond1 & yend_open_tf < 70)));
s = scatter3(gdtx,gdty,gdtz,'k','filled'); s.MarkerFaceAlpha = 0.5; 
axis([0 x(end) 50 90]); daspect([1 1 1]); view(3)

load('mhd0.mat'); figure
cond1 = yend_t0 > 55 & tubetag1 ~= tubetag2 & xend_t0 < 60 & xend_t0 > 20;
scatter(xend_t0(yend_t0 > 70 & cond1),yend_t0(cond1 & yend_t0 > 70),'g.');
hold on
scatter(xend_t0(cond1 & yend_t0 < 70),yend_t0(cond1 & yend_t0 < 70),'r.');
scatter(xend_conv_tf(cond1 & yend_conv_tf > 70),yend_conv_tf(cond1 & yend_conv_tf > 70),'g.');
scatter(xend_conv_tf(cond1 & yend_conv_tf < 70),yend_conv_tf(cond1 & yend_conv_tf < 70),'r.');
reso = 200;
flines_x = x_fluelt_t0(cond1 & yend_t0 > 70); flines_y = y_fluelt_t0(cond1 & yend_t0 > 70);
flines_z = 60*ones(size(flines_y));
h1 = streamline(stream3(x,y,z,bx,by,bz,flines_x(1:reso:end),flines_y(1:reso:end),flines_z(1:reso:end)));
set(h1,'Color','g')
gdtx = x_fluelt_t0(cond1 & yend_open_tf > 70); gdty = y_fluelt_t0(cond1 & yend_open_tf > 70);
gdtz = 61*ones(size(y_fluelt_t0(cond1 & yend_open_tf > 70)));
s = scatter3(gdtx,gdty,gdtz,'k','filled'); s.MarkerFaceAlpha = 0.5; 
flines_x = x_fluelt_t0(cond1 & yend_t0 < 70);flines_y = y_fluelt_t0(cond1 & yend_t0 < 70);
flines_z = 60*ones(size(flines_y));
h1 = streamline(stream3(x,y,z,bx,by,bz,flines_x(1:reso:end),flines_y(1:reso:end),flines_z(1:reso:end)));
set(h1,'Color','r')
gdtx = x_fluelt_t0(cond1 & yend_open_tf < 70); gdty = y_fluelt_t0(cond1 & yend_open_tf < 70);
gdtz = 61*ones(size(y_fluelt_t0(cond1 & yend_open_tf < 70)));
s = scatter3(gdtx,gdty,gdtz,'k','filled'); s.MarkerFaceAlpha = 0.5; 
axis([0 x(end) 50 90]); daspect([1 1 1]); view(3)
 
% load('mhd10.mat');
% figure; hold on
% reso = 50;
% flines_x = rect_flueltst0_cx(cond1 & rect_yend_t0 > 70);
% flines_y = rect_flueltst0_cy(cond1 & rect_yend_t0 > 70);
% flines_z = 60*ones(size(flines_y));
% h1 = streamline(stream3(x,y,z,bx,by,bz,flines_x(1:reso:end),flines_y(1:reso:end),flines_z(1:reso:end)));
% set(h1,'Color','g')
% flines_x = rect_flueltst0_cx(cond1 & rect_yend_t0 < 70);
% flines_y = rect_flueltst0_cy(cond1 & rect_yend_t0 < 70);
% flines_z = 60*ones(size(flines_y));
% h1 = streamline(stream3(x,y,z,bx,by,bz,flines_x(1:reso:end),flines_y(1:reso:end),flines_z(1:reso:end)));
% set(h1,'Color','r')
% axis([0 x(end) 50 90 0 40])
% daspect([1 1 1])
% plot3(mmm,ppp,nnn,'k','LineWidth',3)

