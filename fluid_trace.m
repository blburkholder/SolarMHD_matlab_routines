% for i = [0,10]
%  prepare3(i)
%  mentrop_prepare3(i)
% end

load('mhd0.mat'); load('flt_mentrop0.mat');

zfin11 = zfin1;
zfin22 = zfin2;
xfin11 = xfin1;
xfin22 = xfin2;
yfin11 = yfin1;
yfin22 = yfin2;
[ffx1,ffy1] = fix_coords(x,y,xfin1,yfin1);
[ffx2,ffy2] = fix_coords(x,y,xfin2,yfin2);
cxt = xpo;
cyt = ypo;
czt = zpo;
xendiff = zeros(size(ffx1));
yendiff = zeros(size(ffx1));
xend1 = zeros(size(ffx1));
yend1 = zeros(size(ffx1));
xend2 = zeros(size(ffx1));
yend2 = zeros(size(ffx1));
xstart = zeros(size(ffx1));
ystart = zeros(size(ffx1));

for l=10
  load(['mhd',num2str(l),'.mat']); load(['flt_mentrop',num2str(l),'.mat']);

  [xf1,yf1] = velocity_perturbation(x,y,ffx1,ffy1,time);
  [xf2,yf2] = velocity_perturbation(x,y,ffx2,ffy2,time);
  [fx1,fy1] = fix_coords(x,y,xfin1,yfin1);
  [fx2,fy2] = fix_coords(x,y,xfin2,yfin2);

  xf1 = xf1';
  xf2 = xf2';
  yf1 = yf1';
  yf2 = yf2';

  cxt1 = xpo;
  cyt1 = ypo;
  czt1 = zpo;

  for i = 1:length(xendiff(:))
    if (zfin1(i) < 1 && zfin11(i) < 1) || (zfin2(i) < 1 && zfin22(i) > 1 && zfin1(i) < 1)
      xendiff(i) = xf1(i) - fx1(i);
      yendiff(i) = yf1(i) - fy1(i);
      xstart(i) = ffx1(i);
      ystart(i) = ffy1(i);
      xend1(i) = fx1(i);
      yend1(i) = fy1(i);
      xend2(i) = xf1(i);
      yend2(i) = yf1(i);
    elseif (zfin2(i) < 1 && zfin22(i) < 1) || (zfin1(i) < 1 && zfin11(i) > 1 && zfin2(i) < 1)
      xendiff(i) = xf2(i) - fx2(i);
      yendiff(i) = yf2(i) - fy2(i);
      xstart(i) = ffx2(i);
      ystart(i) = ffy2(i);
      xend1(i) = fx2(i);
      yend1(i) = fy2(i);
      xend2(i) = xf2(i);
      yend2(i) = yf2(i);
    else
      xendiff(i) = NaN;
      yendiff(i) = NaN;
      xstart(i) = NaN;
      ystart(i) = NaN;
      xend1(i) = NaN;
      yend1(i) = NaN;
      xend2(i) = NaN;
      yend2(i) = NaN;    
    end
    i/length(xendiff(:))
  end

  figure
  dist = sqrt(xendiff.^2+yendiff.^2);
  scatter(cxt1(:),cyt1(:),[],dist(:),'.');
  colorbar
  title('footpoint dist from ideal evolution')
  daspect([1 1 1])
end

bx_top = interp3(x,y,z,bx,cxt,cyt,czt);
by_top = interp3(x,y,z,by,cxt,cyt,czt);
bz_top = interp3(x,y,z,bz,cxt,cyt,czt);
open = abs(bz_top) > abs(by_top) & abs(bz_top) > abs(bx_top);

figure
scatter(cxt(:),cyt(:),[],dist(:),'.')
colorbar
daspect([1 1 1])
title('footpoint distance from ideal evolution')

figure
thresh = 15;
scatter(cxt1(open),cyt1(open),[],dist(open),'.')
colorbar
daspect([1 1 1])
title('footpoint distance from ideal evolution')

% figure
% scatter(xstart(dist > thresh),ystart(dist > thresh),'.');
% hold on
% scatter(xend1(dist > thresh),yend1(dist > thresh),[],sqrt(xend2(dist > thresh).^2+yend2(dist > thresh).^2),'.');
% reso = 10;
% [xx,yy] = meshgrid(x(1:reso:end),y(1:reso:end));
% h = quiver(xx,yy,sx(1:reso:end,1:reso:end,1)./rho(1:reso:end,1:reso:end,1),...
%     sy(1:reso:end,1:reso:end,1)./rho(1:reso:end,1:reso:end,1));
% set(h,'Color','k')
% axis([0 x(end) 0 y(end)])
% daspect([1 1 1])
% title('photospheric footpoint open flux')
% colormap(jet)
% 
% figure
% scatter(xstart(dist > thresh),ystart(dist > thresh),'.');
% hold on
% scatter(xend2(dist > thresh),yend2(dist > thresh),[],sqrt(xend2(dist > thresh).^2+yend2(dist > thresh).^2),'.');
% [xx,yy] = meshgrid(x(1:reso:end),y(1:reso:end));
% h = quiver(xx,yy,sx(1:reso:end,1:reso:end,1)./rho(1:reso:end,1:reso:end,1),...
%     sy(1:reso:end,1:reso:end,1)./rho(1:reso:end,1:reso:end,1));
% set(h,'Color','k')
% axis([0 x(end) 0 y(end)])
% daspect([1 1 1])
% title('photospheric footpoint fluid trace')
% colormap(jet)

figure
scatter(xstart(:),ystart(:),'k.')
hold on
scatter(xend1(:),yend1(:),[],sqrt((xend2(:)-xend1(:)).^2+(yend1(:) - yend2(:)).^2),'.')
reso = 10;
[xx,yy] = meshgrid(x(1:reso:end),y(1:reso:end));
h = quiver(xx,yy,sx(1:reso:end,1:reso:end,1)./rho(1:reso:end,1:reso:end,1),...
    sy(1:reso:end,1:reso:end,1)./rho(1:reso:end,1:reso:end,1));
set(h,'Color','k')
axis([0 x(end) 0 y(end)])
daspect([1 1 1])
title('photospheric footpoint open flux')
colormap(jet)

figure
scatter(xstart(:),ystart(:),'k.')
hold on
scatter(xend2(:),yend2(:),'g.')
[xx,yy] = meshgrid(x(1:reso:end),y(1:reso:end));
h = quiver(xx,yy,sx(1:reso:end,1:reso:end,1)./rho(1:reso:end,1:reso:end,1),...
    sy(1:reso:end,1:reso:end,1)./rho(1:reso:end,1:reso:end,1));
set(h,'Color','k')
axis([0 x(end) 0 y(end)])
daspect([1 1 1])
title('photospheric footpoint fluid trace')
colormap(jet)

load('mhd10.mat');
cond1 = cxt1 > 6 & cxt1 < 77 & cyt1 > 55 & dist > thresh;
figure
scatter(xstart(cond1 & ystart > 70),ystart(cond1 & ystart > 70),'b.');
hold on
scatter(xstart(cond1 & ystart < 70),ystart(cond1 & ystart < 70),'m.');
scatter(xend2(cond1 & ystart > 70),yend2(cond1 & ystart > 70),'b.');
scatter(xend2(cond1 & ystart < 70),yend2(cond1 & ystart< 70),'m.');
scatter(xend1(cond1),yend1(cond1),'g.');
reso = 10;
flines_x = cxt1(cond1 & yend1 > 70);
flines_y = cyt1(cond1 & yend1 > 70);
flines_z = czt1(cond1 & yend1 > 70);
h1 = streamline(stream3(x,y,z,bx,by,bz,flines_x(1:reso:end),flines_y(1:reso:end),flines_z(1:reso:end)));
set(h1,'Color','r')
scatter3(cxt1(cond1 & yend1 > 70),cyt1(cond1 & yend1 > 70),czt1(cond1 & yend1 > 70),'k.')
flines_x = cxt1(cond1 & yend1 < 70);
flines_y = cyt1(cond1 & yend1 < 70);
flines_z = czt1(cond1 & yend1 < 70);
h1 = streamline(stream3(x,y,z,bx,by,bz,flines_x(1:reso:end),flines_y(1:reso:end),flines_z(1:reso:end)));
scatter3(cxt1(cond1 & yend1 < 70),cyt1(cond1 & yend1 < 70),czt1(cond1 & yend1 < 70),'k.')
set(h1,'Color','c')
axis([0 x(end) 50 90])
daspect([1 1 1])
view(3)

load('mhd0.mat');
figure
scatter(xstart(cond1 & ystart > 70),ystart(cond1 & ystart > 70),'b.');
hold on
scatter(xstart(cond1 & ystart < 70),ystart(cond1 & ystart < 70),'m.');
scatter(xend2(cond1 & ystart > 70),yend2(cond1 & ystart > 70),'b.');
scatter(xend2(cond1 & ystart < 70),yend2(cond1 & ystart< 70),'m.');
flines_x = cxt(cond1 & yend1 > 70);
flines_y = cyt(cond1 & yend1 > 70);
flines_z = czt(cond1 & yend1 > 70);
h1 = streamline(stream3(x,y,z,bx,by,bz,flines_x(1:reso:end),flines_y(1:reso:end),flines_z(1:reso:end)));
set(h1,'Color','r')
scatter3(cxt(cond1 & yend1 > 70),cyt(cond1 & yend1 > 70),czt(cond1 & yend1 > 70),'k.')
flines_x = cxt(cond1 & yend1 < 70);
flines_y = cyt(cond1 & yend1 < 70);
flines_z = czt(cond1 & yend1 < 70);
h1 = streamline(stream3(x,y,z,bx,by,bz,flines_x(1:reso:end),flines_y(1:reso:end),flines_z(1:reso:end)));
set(h1,'Color','c')
scatter3(cxt(cond1 & yend1 < 70),cyt(cond1 & yend1 < 70),czt(cond1 & yend1 < 70),'k.')
axis([0 x(end) 50 90])
daspect([1 1 1])
view(3)

load('mhd1.mat');
figure
scatter(xstart(cond1 & ystart > 70),ystart(cond1 & ystart > 70),'b.');
hold on
scatter(xstart(cond1 & ystart < 70),ystart(cond1 & ystart < 70),'m.');
reso = 5;
flines_x = cxt(cond1 & yend1 > 70);
flines_y = cyt(cond1 & yend1 > 70);
flines_z = czt(cond1 & yend1 > 70);
h1 = streamline(stream3(x,y,z,bx,by,bz,flines_x(1:reso:end),flines_y(1:reso:end),flines_z(1:reso:end)));
set(h1,'Color','r')
flines_x = cxt(cond1 & yend1 < 70);
flines_y = cyt(cond1 & yend1 < 70);
flines_z = czt(cond1 & yend1 < 70);
h1 = streamline(stream3(x,y,z,bx,by,bz,flines_x(1:reso:end),flines_y(1:reso:end),flines_z(1:reso:end)));
set(h1,'Color','c')
axis([0 x(end) 50 90 0 20])
daspect([1 1 1])
plot3(mmm,ppp,nnn,'k','LineWidth',3)

