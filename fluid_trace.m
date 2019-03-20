% for i = [0,10]
%  prepare3(i)
%  mentrop_prepare3(i)
% end

load('mhd0.mat'); load('flt_mentrop0.mat');
wy = [-flipud(y(4:end));y];
wbx = [-fliplr(flipud(bx(3:end,:,:)));bx(2:end,:,:)];
wby = [-fliplr(flipud(by(3:end,:,:)));by(2:end,:,:)];
wbz = [fliplr(flipud(bz(3:end,:,:)));bz(2:end,:,:)];
wsx = [-fliplr(flipud(sx(3:end,:,2)));sx(2:end,:,2)];
wsy = [-fliplr(flipud(sy(3:end,:,2)));sy(2:end,:,2)];
wrho = [fliplr(flipud(rho(3:end,:,2)));rho(2:end,:,2)];
vx = wsx./wrho;
vy = wsy./wrho;

%use line symmetric boundaries to get every footpoint in 1st quadrant
[xend1_t0,yend1_t0] = fix_coords(x,y,xfin1,yfin1);
zend1_t0 = zfin1;
[xend2_t0,yend2_t0] = fix_coords(x,y,xfin2,yfin2);
zend2_t0 = zfin2;
%keep starting points of field lines at top boundary
x_fluelt_t0 = xpo;
y_fluelt_t0 = ypo;
z_fluelt_t0 = zpo;

xendiff = zeros(size(xfin1));
yendiff = zeros(size(xfin1));
xend_open_tf = zeros(size(xfin1));
yend_open_tf = zeros(size(xfin1));
xend_conv_tf = zeros(size(xfin1));
yend_conv_tf = zeros(size(xfin1));
xend_t0 = zeros(size(xfin1));
yend_t0 = zeros(size(xfin1));
  
  l=10;
  load(['mhd',num2str(l),'.mat']); load(['flt_mentrop',num2str(l),'.mat']);

  %take initial location of both ends of the field line (whether on p-sphere or not)
  %and push with photospheric convection  
  [xend1_conv_tf,yend1_conv_tf] = velocity_perturbation(x,y,xend1_t0,yend1_t0,time);
  [xend2_conv_tf,yend2_conv_tf] = velocity_perturbation(x,y,xend2_t0,yend2_t0,time);
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

  for i = 1:length(xendiff(:))
    %if end1 was on the photosphere for all times
    if (zend1_t0(i) < 1 && zend1_tf(i) < 1)
      xendiff(i) = xend1_conv_tf(i) - xend1_open_tf(i); yendiff(i) = yend1_conv_tf(i) - yend1_open_tf(i);
      xend_t0(i) = xend1_t0(i); yend_t0(i) = yend1_t0(i);
      xend_open_tf(i) = xend1_open_tf(i); yend_open_tf(i) = yend1_open_tf(i);
      xend_conv_tf(i) = xend1_conv_tf(i); yend_conv_tf(i) = yend1_conv_tf(i);
    %if end2 was on the photosphere for all times
    elseif (zend2_t0(i) < 1 && zend2_tf(i) < 1)
      xendiff(i) = xend2_conv_tf(i) - xend2_open_tf(i); yendiff(i) = yend2_conv_tf(i) - yend2_open_tf(i);
      xend_t0(i) = xend2_t0(i); yend_t0(i) = yend2_t0(i);
      xend_open_tf(i) = xend2_open_tf(i); yend_open_tf(i) = yend2_open_tf(i);
      xend_conv_tf(i) = xend2_conv_tf(i); yend_conv_tf(i) = yend2_conv_tf(i);
    else
      xendiff(i) = NaN;      yendiff(i) = NaN;
      xend_t0(i) = NaN;      yend_t0(i) = NaN;
      xend_open_tf(i) = NaN;      yend_open_tf(i) = NaN;
      xend_conv_tf(i) = NaN;      yend_conv_tf(i) = NaN;    
    end
   % i/length(xendiff(:))
  end

rect_flueltst0_cy = [-fliplr(flipud(y_fluelt_t0(2:end,:)));y_fluelt_t0];
rect_flueltst0_cx = [x(end-1)-fliplr(flipud(x_fluelt_t0(2:end,:)));x_fluelt_t0];

rect_flueltstf_cy = [-fliplr(flipud(y_fluelt_tf(2:end,:)));y_fluelt_tf];
rect_flueltstf_cx = [x(end-1)-fliplr(flipud(x_fluelt_tf(2:end,:)));x_fluelt_tf];

rect_xend_t0 = [fliplr(flipud(x(end-1) - xend_t0(2:end,:)));xend_t0];
rect_yend_t0 = [fliplr(flipud(-yend_t0(2:end,:)));yend_t0];
rect_xend_tf = [fliplr(flipud(x(end-1) - xend_open_tf(2:end,:)));xend_open_tf];
rect_yend_tf = [fliplr(flipud(-yend_open_tf(2:end,:)));yend_open_tf];
rect_xend_tflu = [fliplr(flipud(x(end-1) - xend_conv_tf(2:end,:)));xend_conv_tf];
rect_yend_tflu = [fliplr(flipud(-yend_conv_tf(2:end,:)));yend_conv_tf];

tubetag1 = zeros(size(rect_flueltst0_cy));
tubetag2 = zeros(size(rect_flueltst0_cy));

ti_tube1 = rect_xend_t0 > 20 & rect_xend_t0 < 60 & ((rect_yend_t0 > 40 & rect_yend_t0 < 70) | (rect_yend_t0 < -40 & rect_yend_t0 > -70));
ti_tube2 = rect_yend_t0 > 70 | rect_yend_t0 < -70;
ti_tube3 = (rect_yend_t0 > 12 & rect_yend_t0 < 34) | (rect_yend_t0 < -12 & rect_yend_t0 > -34);
ti_tube4 = rect_yend_t0 < 12 & rect_yend_t0 > -12;

tf_tube1 = rect_xend_tf > 20 & rect_xend_tf < 60 & ((rect_yend_tf > 40 & rect_yend_tf < 70) | (rect_yend_tf < -40 & rect_yend_tf > -70));
tf_tube2 = rect_yend_tf > 70 | rect_yend_tf < -70;
tf_tube3 = (rect_yend_tf > 12 & rect_yend_tf < 34) | (rect_yend_tf < -12 & rect_yend_tf > -34);
tf_tube4 = rect_yend_tf < 12 & rect_yend_tf > -12;

tubetag1(ti_tube1) = 1;
tubetag1(ti_tube2) = 2;
tubetag1(ti_tube3) = 3;
tubetag1(ti_tube4) = 4;

tubetag2(tf_tube1) = 1;
tubetag2(tf_tube2) = 2;
tubetag2(tf_tube3) = 3;
tubetag2(tf_tube4) = 4;

ti_xtube1 = rect_flueltst0_cx(ti_tube1); ti_ytube1 = rect_flueltst0_cy(ti_tube1);
ti_xtube2 = rect_flueltst0_cx(ti_tube2); ti_ytube2 = rect_flueltst0_cy(ti_tube2);
ti_xtube3 = rect_flueltst0_cx(ti_tube3); ti_ytube3 = rect_flueltst0_cy(ti_tube3);
ti_xtube4 = rect_flueltst0_cx(ti_tube4); ti_ytube4 = rect_flueltst0_cy(ti_tube4);

tf_xtube1 = rect_flueltstf_cx(tf_tube1); tf_ytube1 = rect_flueltstf_cy(tf_tube1);
tf_xtube2 = rect_flueltstf_cx(tf_tube2); tf_ytube2 = rect_flueltstf_cy(tf_tube2);
tf_xtube3 = rect_flueltstf_cx(tf_tube3); tf_ytube3 = rect_flueltstf_cy(tf_tube3);
tf_xtube4 = rect_flueltstf_cx(tf_tube4); tf_ytube4 = rect_flueltstf_cy(tf_tube4);

% figure; hold on
% scatter(rect_xend_t0(ti_tube1),rect_yend_t0(ti_tube1),'r.')
% scatter(rect_xend_t0(ti_tube2),rect_yend_t0(ti_tube2),'g.')
% scatter(rect_xend_t0(ti_tube3),rect_yend_t0(ti_tube3),'b.')
% scatter(rect_xend_t0(ti_tube4),rect_yend_t0(ti_tube4),'y.')
% daspect([1 1 1])
% axis([x(1) x(end) -y(end) y(end)])
% xlabel('x')
% ylabel('y')
% title('open flux at photosphere')
% 
% figure; hold on
% scatter(ti_xtube1,ti_ytube1,'r.')
% scatter(ti_xtube2,ti_ytube2,'g.')
% scatter(ti_xtube3,ti_ytube3,'b.')
% scatter(ti_xtube4,ti_ytube4,'y.')
% daspect([1 1 1])
% axis([x(1) x(end) -y(end) y(end)])
% xlabel('x')
% ylabel('y')
% title('open flux tubes initial equilibrium')
% 
% figure; hold on
% scatter(tf_xtube1,tf_ytube1,'r.')
% scatter(tf_xtube2,tf_ytube2,'g.')
% scatter(tf_xtube3,tf_ytube3,'b.')
% scatter(tf_xtube4,tf_ytube4,'y.')
% flipx = rect_flueltstf_cx(tubetag1 ~= tubetag2);
% flipy = rect_flueltstf_cy(tubetag1 ~= tubetag2);
% scatter(flipx,flipy,'k.','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1)
% daspect([1 1 1])
% axis([x(1) x(end) -y(end) y(end)])
% xlabel('x')
% ylabel('y')
% title('reconnection between flux tubes convected coordinates')
% 
% figure; hold on
% scatter(tf_xtube1,tf_ytube1,'r.')
% scatter(tf_xtube2,tf_ytube2,'g.')
% scatter(tf_xtube3,tf_ytube3,'b.')
% scatter(tf_xtube4,tf_ytube4,'y.')
% daspect([1 1 1])
% axis([x(1) x(end) -y(end) y(end)])
% xlabel('x')
% ylabel('y')
% title('open flux tubes t = 146')
% 
% figure; hold on
% scatter(ti_xtube1,ti_ytube1,'r.')
% scatter(ti_xtube2,ti_ytube2,'g.')
% scatter(ti_xtube3,ti_ytube3,'b.')
% scatter(ti_xtube4,ti_ytube4,'y.')
% flipx = rect_flueltst0_cx(tubetag1 ~= tubetag2);
% flipy = rect_flueltst0_cy(tubetag1 ~= tubetag2);
% scatter(flipx,flipy,'w.','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
% daspect([1 1 1])
% axis([x(1) x(end) -y(end) y(end)])
% xlabel('x')
% ylabel('y')
% title('reconnection between flux tubes initial coordinates')

% figure; hold on
% pcolor(x,wy,sqrt(vx.^2+vy.^2))
% colormap(gray)
% colorbar
% shading interp
% scatter(rect_xend_t0(ti_tube1),rect_yend_t0(ti_tube1),'r.')
% scatter(rect_xend_t0(ti_tube2),rect_yend_t0(ti_tube2),'g.')
% scatter(rect_xend_t0(ti_tube3),rect_yend_t0(ti_tube3),'b.')
% scatter(rect_xend_t0(ti_tube4),rect_yend_t0(ti_tube4),'y.')
% [xxx,yyy] = meshgrid(x,wy);
% h = streamslice(xxx,yyy,vx,vy);
% set(h,'Color','k');
% reso = 500;
% h = streamline(stream3(x,wy,z,-wbx,-wby,-wbz,ti_xtube4(1:reso:end),ti_ytube4(1:reso:end),60*ones(size(ti_ytube4(1:reso:end)))));
% set(h,'Color','y')
% reso = 2000;
% h = streamline(stream3(x,wy,z,-wbx,-wby,-wbz,ti_xtube3(1:reso:end),ti_ytube3(1:reso:end),60*ones(size(ti_ytube3(1:reso:end)))));
% set(h,'Color','b')
% h = streamline(stream3(x,wy,z,wbx,wby,wbz,ti_xtube2(1:reso:end),ti_ytube2(1:reso:end),60*ones(size(ti_ytube2(1:reso:end)))));
% set(h,'Color','g')
% h = streamline(stream3(x,wy,z,wbx,wby,wbz,ti_xtube1(1:reso:end),ti_ytube1(1:reso:end),60*ones(size(ti_ytube1(1:reso:end)))));
% set(h,'Color','r')
% scatter3(ti_xtube1,ti_ytube1,61*ones(size(ti_xtube1)),'r.','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75)
% scatter3(ti_xtube2,ti_ytube2,61*ones(size(ti_xtube2)),'g.','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75)
% scatter3(ti_xtube3,ti_ytube3,61*ones(size(ti_xtube3)),'b.','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75)
% scatter3(ti_xtube4,ti_ytube4,61*ones(size(ti_xtube4)),'y.','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75)
% axis([0 92.8 -92.8 92.8])
% daspect([1 1 1])
% xlabel('x'); ylabel('y'); zlabel('z')
% view(3)

% bx_top = interp3(x,y,z,bx,cxt,cyt,czt);
% by_top = interp3(x,y,z,by,cxt,cyt,czt);
% bz_top = interp3(x,y,z,bz,cxt,cyt,czt);
% open = abs(bz_top) > abs(by_top) & abs(bz_top) > abs(bx_top);

load('mhd10.mat');
figure
cond1 = rect_yend_t0 > 55 & tubetag1 ~= tubetag2 & rect_xend_t0 < 60 & rect_xend_t0 > 20;
scatter(rect_xend_t0(rect_yend_t0 > 70 & cond1),rect_yend_t0(cond1 & rect_yend_t0 > 70),'g.');
hold on
scatter(rect_xend_t0(cond1 & rect_yend_t0 < 70),rect_yend_t0(cond1 & rect_yend_t0 < 70),'r.');
scatter(rect_xend_tflu(cond1 & rect_yend_tflu > 70),rect_yend_tflu(cond1 & rect_yend_tflu > 70),'g.');
scatter(rect_xend_tflu(cond1 & rect_yend_tflu < 70),rect_yend_tflu(cond1 & rect_yend_tflu < 70),'r.');
scatter(rect_xend_tf(cond1),rect_yend_tf(cond1),'k.');
%scatter(rect_xend_tf(cond1 & rect_yend_t0 > 70),rect_yend_tf(cond1 & rect_yend_t0 > 70),'g.');
%scatter(rect_xend_tf(cond1 & rect_yend_t0 < 70),rect_yend_tf(cond1 & rect_yend_t0 < 70),'r.');
reso = 50;
flines_x = rect_flueltstf_cx(cond1 & rect_yend_t0 > 70);
flines_y = rect_flueltstf_cy(cond1 & rect_yend_t0 > 70);
flines_z = 60*ones(size(flines_y));
h1 = streamline(stream3(x,y,z,bx,by,bz,flines_x(1:reso:end),flines_y(1:reso:end),flines_z(1:reso:end)));
set(h1,'Color','g')
scatter3(rect_flueltstf_cx(cond1 & rect_yend_tflu > 70),rect_flueltstf_cy(cond1 & rect_yend_tflu > 70),...
    60*ones(size(rect_flueltstf_cy(cond1 & rect_yend_tflu > 70))),'k.')
flines_x = rect_flueltstf_cx(cond1 & rect_yend_t0 < 70);
flines_y = rect_flueltstf_cy(cond1 & rect_yend_t0 < 70);
flines_z = 60*ones(size(flines_y));
h1 = streamline(stream3(x,y,z,bx,by,bz,flines_x(1:reso:end),flines_y(1:reso:end),flines_z(1:reso:end)));
set(h1,'Color','r')
scatter3(rect_flueltstf_cx(cond1 & rect_yend_tflu < 70),rect_flueltstf_cy(cond1 & rect_yend_tflu < 70),...
    60*ones(size(rect_flueltstf_cy(cond1 & rect_yend_tflu < 70))),'k.')
axis([0 x(end) 50 90])
daspect([1 1 1])
view(3)

load('mhd0.mat');
figure
cond1 = rect_yend_t0 > 55 & tubetag1 ~= tubetag2 & rect_xend_t0 < 60 & rect_xend_t0 > 20;
scatter(rect_xend_t0(rect_yend_t0 > 70 & cond1),rect_yend_t0(cond1 & rect_yend_t0 > 70),'g.');
hold on
scatter(rect_xend_t0(cond1 & rect_yend_t0 < 70),rect_yend_t0(cond1 & rect_yend_t0 < 70),'r.');
scatter(rect_xend_tflu(cond1 & rect_yend_tflu > 70),rect_yend_tflu(cond1 & rect_yend_tflu > 70),'g.');
scatter(rect_xend_tflu(cond1 & rect_yend_tflu < 70),rect_yend_tflu(cond1 & rect_yend_tflu < 70),'r.');
reso = 50;
flines_x = rect_flueltst0_cx(cond1 & rect_yend_t0 > 70);
flines_y = rect_flueltst0_cy(cond1 & rect_yend_t0 > 70);
flines_z = 60*ones(size(flines_y));
h1 = streamline(stream3(x,y,z,bx,by,bz,flines_x(1:reso:end),flines_y(1:reso:end),flines_z(1:reso:end)));
set(h1,'Color','g')
scatter3(rect_flueltst0_cx(cond1 & rect_yend_tflu > 70),rect_flueltst0_cy(cond1 & rect_yend_tflu > 70),...
    60*ones(size(rect_flueltst0_cy(cond1 & rect_yend_tflu > 70))),'k.')
flines_x = rect_flueltst0_cx(cond1 & rect_yend_t0 < 70);
flines_y = rect_flueltst0_cy(cond1 & rect_yend_t0 < 70);
flines_z = 60*ones(size(flines_y));
h1 = streamline(stream3(x,y,z,bx,by,bz,flines_x(1:reso:end),flines_y(1:reso:end),flines_z(1:reso:end)));
set(h1,'Color','r')
scatter3(rect_flueltst0_cx(cond1 & rect_yend_tflu < 70),rect_flueltst0_cy(cond1 & rect_yend_tflu < 70),...
    60*ones(size(rect_flueltst0_cy(cond1 & rect_yend_tflu < 70))),'k.')
axis([0 x(end) 50 90])
daspect([1 1 1])
view(3)
 
load('mhd10.mat');
figure; hold on
reso = 50;
flines_x = rect_flueltst0_cx(cond1 & rect_yend_t0 > 70);
flines_y = rect_flueltst0_cy(cond1 & rect_yend_t0 > 70);
flines_z = 60*ones(size(flines_y));
h1 = streamline(stream3(x,y,z,bx,by,bz,flines_x(1:reso:end),flines_y(1:reso:end),flines_z(1:reso:end)));
set(h1,'Color','g')
flines_x = rect_flueltst0_cx(cond1 & rect_yend_t0 < 70);
flines_y = rect_flueltst0_cy(cond1 & rect_yend_t0 < 70);
flines_z = 60*ones(size(flines_y));
h1 = streamline(stream3(x,y,z,bx,by,bz,flines_x(1:reso:end),flines_y(1:reso:end),flines_z(1:reso:end)));
set(h1,'Color','r')
axis([0 x(end) 50 90 0 40])
daspect([1 1 1])
plot3(mmm,ppp,nnn,'k','LineWidth',3)

