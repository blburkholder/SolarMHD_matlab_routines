% zzz = 149;
% 
% %first load an mhd.mat and a mentrop.mat
% 
% %perpendicular velocity
% pp = sqrt(sx.^2+sy.^2)./rho;
% figure ; pcolor(y,x,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1])
% title(strcat('z = ',num2str(z(zzz)),' v_\perp'))
% 
% %field aligned velocity
% pp = sz./rho;
% figure ; pcolor(y,x,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1])   
% title(strcat('z = ',num2str(z(zzz)),' v_{||}'))
% 
% %parallel current
% [jx,jy,jz,pp] = get_j(nx,ny,nz,bx,by,bz,b0x,b0y,b0z,res,difx,dify,difz,2);
% figure ; pcolor(y,x,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1])  
% title(strcat('z = ',num2str(z(zzz)),' J_{||}'))
% 
% %perpendicular current
% [jx,jy,jz,pp] = get_j(nx,ny,nz,bx,by,bz,b0x,b0y,b0z,res,difx,dify,difz,3);
% figure ; pcolor(y,x,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1])  
% title(strcat('z = ',num2str(z(zzz)),' J_\perp'))
% 
% %specific entropy
% pp = p./(rho).^(5/3);
% figure ; pcolor(y,x,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1]) 
% title(strcat('z = ',num2str(z(zzz)),' specific entropy'))
% 
% %amount that angle wrt z has changed
% bz_n0 = b0z./(sqrt(b0x.^2+b0y.^2+b0z.^2));
% a0 = acos(bz_n0);
% bz_n = bz./(sqrt(bx.^2+by.^2+bz.^2));
% pp = a0 - acos(bz_n);
% figure ; pcolor(y,x,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1]) 
% title(strcat('z = ',num2str(z(zzz)),' misalignment'))
% 
% %all 3 components of field
% pp = bz;
% figure ; pcolor(y,x,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1])
% hold on
% ree = 10;
% quiver(y(1:ree:end),x(1:ree:end),by(1:ree:end,1:ree:end,zzz),bx(1:ree:end,1:ree:end,zzz),'r') 
% title(strcat('z = ',num2str(z(zzz)),' magnetic field'))
% 
% %plasma beta
% pp = p./(bx.^2+by.^2+bz.^2);
% figure ; pcolor(y,x,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1])
% title(strcat('z = ',num2str(z(zzz)),' beta'))
% 
% %density
% pp = rho;
% figure ; pcolor(y,x,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1])
% title(strcat('z = ',num2str(z(zzz)),' density'))

%magnetic perturbation for top boundary
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
vz = interp3(y,x,z,sz./rho,yyfin,xxfin,z(n)*ones(size(zfin1)));
dbx = vx.*(bbx.^2 + bby.^2 + bbz.^2)./VVa;
dby = vy.*(bbx.^2 + bby.^2 + bbz.^2)./VVa;
dbz = vz.*(bbx.^2 + bby.^2 + bbz.^2)./VVa;

figure
reso = 20;
pcolor(ypo,xpo,dbz) ; shading interp ; colorbar
hold on
quiver(ypo(1:reso:end),xpo(1:reso:end),dby(1:reso:end,1:reso:end),dbx(1:reso:end,1:reso:end))
title(strcat('z = ',num2str(z(zzz)),' magnetic perturbation')) ; daspect([1 1 1])

