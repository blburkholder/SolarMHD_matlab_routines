
b = sqrt((bx.^2+by.^2+bz.^2));
Va = b./sqrt((rho));

[six,siy,siz] = meshgrid(xpo,ypo,zp1);
xxfin = zeros(size(xfin1));
yyfin = zeros(size(yfin1));


n = 2;
for i = 1:length(xfin1)
    for j = 1:length(xfin1)
        if xfin1(i,j) > x(end-1)
            xxfin(i,j) = x(end-1) - mod(xfin1(i,j),x(end-1));
        elseif xfin1(i,j) < x(end-1)
            xxfin(i,j) = abs(xfin1(i,j));
        end

        if yfin1(i,j) > y(end-1)
            yyfin(i,j) = y(end-1) - mod(yfin1(i,j),y(end-1));
        elseif yfin1(i,j) < y(end-1)
            yyfin(i,j) = abs(yfin1(i,j));
        end
    end
end

        

bbz = interp3(y,x,z,bz,yyfin,xxfin,z(n)*ones(size(zfin1)));
bby = interp3(y,x,z,by,yyfin,xxfin,z(n)*ones(size(zfin1)));
bbx = interp3(y,x,z,bx,yyfin,xxfin,z(n)*ones(size(zfin1)));
VVa = interp3(y,x,z,Va,yyfin,xxfin,z(n)*ones(size(zfin1)));
vx = interp3(y,x,z,sx./rho,yyfin,xxfin,z(n)*ones(size(zfin1)));
vy = interp3(y,x,z,sy./rho,yyfin,xxfin,z(n)*ones(size(zfin1)));

% bbz = interp3(y,x,z,bz,yfin1,xfin1,z(n)*ones(size(zfin1)));
% bby = interp3(y,x,z,by,yfin1,xfin1,z(n)*ones(size(zfin1)));
% bbx = interp3(y,x,z,bx,yfin1,xfin1,z(n)*ones(size(zfin1)));
% VVa = interp3(y,x,z,Va,yfin1,xfin1,z(n)*ones(size(zfin1)));
% vx = interp3(y,x,z,sx./rho,yfin1,xfin1,z(n)*ones(size(zfin1)));
% vy = interp3(y,x,z,sy./rho,yfin1,xfin1,z(n)*ones(size(zfin1)));

% bbz = interp3(y,x,z,bz,mod(abs(yfin1),y(end-1)),mod(abs(xfin1),x(end-1)),z(n)*ones(size(zfin1)));
% bby = interp3(y,x,z,by,mod(abs(yfin1),y(end-1)),mod(abs(xfin1),x(end-1)),z(n)*ones(size(zfin1)));
% bbx = interp3(y,x,z,bx,mod(abs(yfin1),y(end-1)),mod(abs(xfin1),x(end-1)),z(n)*ones(size(zfin1)));
% VVa = interp3(y,x,z,Va,mod(abs(yfin1),y(end-1)),mod(abs(xfin1),x(end-1)),z(n)*ones(size(zfin1)));
% vx = interp3(y,x,z,sx./rho,mod(abs(yfin1),y(end-1)),mod(abs(xfin1),x(end-1)),z(n)*ones(size(zfin1)));
% vy = interp3(y,x,z,sy./rho,mod(abs(yfin1),y(end-1)),mod(abs(xfin1),x(end-1)),z(n)*ones(size(zfin1)));

dbx = vx.*(bbx.^2 + bby.^2 + bbz.^2)./VVa;
dby = vy.*(bbx.^2 + bby.^2 + bbz.^2)./VVa;


% figure
% pcolor(ypo,xpo,dbx)
% shading interp
% colorbar
% daspect([1 1 1])
% 
% figure
% pcolor(ypo,xpo,dby)
% shading interp
% colorbar
% daspect([1 1 1])

figure
reso = 50;
quiver(ypo(1:reso:end),ypo(1:reso:end),dby(1:res:end,1:res:end),dbx(1:res:end,1:res:end))
%streamslice(ypo(1:reso:end),ypo(1:reso:end),dby(1:res:end,1:res:end),dbx(1:res:end,1:res:end))
