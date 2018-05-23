% for i = [0,10]
%   prepare3(i)
%   mentrop_prepare2(i)
%   flprepare2(i)
% end

load('flt0.mat')
load('mhd0.mat')
load('mentrop0.mat')

%s = 1-bottom 2-middle 3-top
s = 3;

slice = squeeze(fepo(:,:,s,:));
coor_x = slice(:,:,1);
coor_y = slice(:,:,2);
coor_z = slice(:,:,3);

fx = interp2(xpo,ypo,xfin1,coor_x',coor_y','nearest');
fy = interp2(xpo,ypo,yfin1,coor_x',coor_y','nearest');

for i = 1:length(fx)
    for j = 1:length(fy)
        if fx(i,j) > x(end)
            fx(i,j) = x(end) - mod(fx(i,j),x(end));
            fy(i,j) = y(end) - fy(i,j);
        end
        if fx(i,j) < x(1)
            fx(i,j) = abs(fx(i,j));
            fy(i,j) = y(end) - fy(i,j);
        end
        if fy(i,j) > y(end)
            fy(i,j) = y(end) - mod(fy(i,j),y(end));
            fx(i,j) = x(end) - fx(i,j);
        end
        if fy(i,j) < y(1)
            fy(i,j) = abs(fy(i,j));
            fx(i,j) = x(end) - fx(i,j);
        end
    end
end

bzend = interp2(x,y,bz(:,:,2),fx,fy);

%for l = 6
l=10;

load(['flt',num2str(l),'.mat'])
load(['mhd',num2str(l),'.mat'])
load(['mentrop',num2str(l),'.mat'])

fx1 = interp2(xpo,ypo,xfin1,coor_x1',coor_y1','nearest');
fy1 = interp2(xpo,ypo,yfin1,coor_x1',coor_y1','nearest');
ftjparl = ftjpar;
ftjparl(abs(ftjparl) > 16) = 16*sign(ftjparl(abs(ftjparl) > 16));
c1 = interp2(xpo,ypo,ftjparl,coor_x1',coor_y1','nearest');

for i = 1:length(fx1)
    for j = 1:length(fy1)
        if fx1(i,j) > x(end)
            fx1(i,j) = x(end) - mod(fx1(i,j),x(end));
            fy1(i,j) = y(end) - fy1(i,j);
        end
        if fx1(i,j) < x(1)
            fx1(i,j) = abs(fx1(i,j));
            fy1(i,j) = y(end) - fy1(i,j);
        end
        if fy1(i,j) > y(end)
            fy1(i,j) = y(end) - mod(fy1(i,j),y(end));
            fx1(i,j) = x(end) - fx1(i,j);
        end
        if fy1(i,j) < y(1)
            fy1(i,j) = abs(fy1(i,j));
            fx1(i,j) = x(end) - fx1(i,j);
        end
    end
end

bzend1 = interp2(x,y,bz(:,:,2),fx1,fy1);

slice1 = squeeze(fepo(:,:,s,:));
coor_x1 = slice1(:,:,1);
coor_y1 = slice1(:,:,2);
coor_z1 = slice1(:,:,3);

cxt1 = coor_x1';
cyt1 = coor_y1';
bzendiff = bzend - bzend1;

figure
h = scatter(cxt1(:),cyt1(:),[],bzendiff(:),'.');
%pcolor(cxt,cyt,bzendiff)
%shading interp
colorbar
title('closed top')

figure
surface(cxt,cyt,coor_z1-coor_z)
shading interp
colorbar
title('closed top')

figure
h = scatter(fx(:),fy(:),'.');
hold on
k = scatter(fx1(:),fy1(:),'.');
title('closed top')

figure
scatter(cxt1(:),cyt1(:),[],c1(:),'.')
colorbar

bzendiff(abs(c1) < 7) = 0;
bzendiff(abs(bzendiff)> 10) = 10*sign(bzendiff(abs(bzendiff) > 10));
figure
scatter(cxt1(:),cyt1(:),[],bzendiff(:),'.')
colorbar





