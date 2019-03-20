  %load('mentrop10.mat')

%   reso = 1;
%   %ftepar(ftepar > 0.02) = 0.02;
%   %ftepar(ftepar < - 0.01) = -0.01;
%   wftepar = [fliplr(flipud(ftepar(3:end,:)));ftepar(2:end,:)];
%   wwy = [fliplr(flipud(-ypo(3:end,:)));ypo(2:end,:)];
%   [wwwx,wwwy] = meshgrid(xpo,wwy);
%   figure
%   %scatter(wwwx(1:10:end),wwwy(1:10:end),[1],wftepar(1:10:end))
%   pcolor(wwwx(1:reso:end,1:reso:end),wwwy(1:reso:end,1:reso:end),wftepar(1:reso:end,1:reso:end))
%   %colormap(jet)
%   shading interp
%   colorbar
%   daspect([1 1 1])
%   xlabel('X')
%   ylabel('Y')
%   title('integrated E_{||}')

yellow_blue = ftepar(1:2500,1:5000);
figure
pcolor(yellow_blue)
shading interp
colorbar
figure
straight_yb = yellow_blue(:);
histogram(straight_yb(straight_yb > 0.02))


