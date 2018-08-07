%for i = 0:10
  %prepare3(i)
  %mentrop_prepare3(i)
%end

load('mhd0.mat'); load('flt_mentrop0.mat');

zfin11 = zfin1;
zfin22 = zfin2;
[ffx1,ffy1] = fix_coords(x,y,xfin1,yfin1);
[ffx2,ffy2] = fix_coords(x,y,xfin2,yfin2);
cxt = xpo;
cyt = ypo;
czt = zpo;
xendiff = zeros(size(ffx1));
yendiff = zeros(size(ffx1));

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
    if (zfin1(i) < 1 && zfin11(i) < 1)
      xendiff(i) = xf1(i) - fx1(i);
      yendiff(i) = yf1(i) - fy1(i);
    elseif (zfin2(i) < 1 && zfin22(i) < 1)
      xendiff(i) = xf2(i) - fx2(i);
      yendiff(i) = yf2(i) - fy2(i);
    elseif (zfin1(i) < 1 && zfin11(i) > 1 && zfin2(i) < 1)
      xendiff(i) = xf2(i) - fx2(i);
      yendiff(i) = yf2(i) - fy2(i);
    elseif (zfin2(i) < 1 && zfin22(i) > 1 && zfin1(i) < 1)
      xendiff(i) = xf1(i) - fx1(i);
      yendiff(i) = yf1(i) - fy1(i);
    else
      xendiff(i) = NaN;
      yendiff(i) = NaN;
    end
  end

%   figure
%   dist = sqrt(xendiff(:).^2+yendiff(:).^2);
%   scatter(cxt1(:),cyt1(:),[],dist(:),'.');
%   colorbar
%   title('footpoint dist from ideal evolution')
%   daspect([1 1 1])

%   figure
%   scatter(fx1(:),fy1(:),[],'r.')
%   axis([x(1) x(end) y(1) y(end)])
%   figure
%   scatter(xf(:),yf(:),[],'b.')
%   axis([x(1) x(end) y(1) y(end)])
  
end
figure
scatter(cxt(:),cyt(:),[],sqrt(xendiff(:).^2+yendiff(:).^2),'.')
colorbar
daspect([1 1 1])
title('footpoint distance from ideal evolution')

% figure
% scatter(fx(:),fy(:),'.');
% hold on
% k = scatter(fx1(:),fy1(:),'.');
% reso = 10;
% [xx,yy] = meshgrid(x(1:reso:end),y(1:reso:end));
% h = quiver(xx,yy,sx(1:reso:end,1:reso:end,1)./rho(1:reso:end,1:reso:end,1),...
%     sy(1:reso:end,1:reso:end,1)./rho(1:reso:end,1:reso:end,1));
% set(h,'Color','k')
% axis([0 x(end) 0 y(end)])
% daspect([1 1 1])
% title('photospheric footpoint open flux')
% 
% figure
% scatter(fx(:),fy(:),'.');
% hold on
% k = scatter(xf(:),yf(:),'.');
% reso = 10;
% [xx,yy] = meshgrid(x(1:reso:end),y(1:reso:end));
% h = quiver(xx,yy,sx(1:reso:end,1:reso:end,1)./rho(1:reso:end,1:reso:end,1),...
%     sy(1:reso:end,1:reso:end,1)./rho(1:reso:end,1:reso:end,1));
% set(h,'Color','k')
% axis([0 x(end) 0 y(end)])
% daspect([1 1 1])
% title('photospheric footpoint fluid trace')
% 
