%  prepare3(10)
  %mentrop_prepare2(15)
  %load('mhd15.mat')
  %load('mentrop15.mat')

%first load an mhd.mat and a mentrop.mat from same time in simulation

% for zzz = [90,149]
% %nz 253
%for zzz = [90,149]
% [jx,jy,jz] = get_j(nx,ny,nz,bx,by,bz,b0x,b0y,b0z,difx,dify,difz);
% [e_par,j_par,j_perp,ex,ey,ez] = get_j2(nx,ny,nz,res,jx,jy,jz,bx,by,bz,sx./rho,sy./rho,sz./rho);
%[dbx,dby,dbz,poynt] = magn_perturb(bx,by,bz,rho,xpo,ypo,xfin1,yfin1,zfin1,x,y,z,sx,sy,sz,ex,ey,ez);
ree = 40;
%ftvol(ftvol > 400) = 400;
for zzz = 149
% 
    %%%perpendicular velocity
    %%%|v x B|/|B|
    %pp = sqrt((sy.*bz - sz.*by).^2 + (sz.*bx - sx.*bz).^2 + (sx.*by - sy.*bx).^2)./(rho.*sqrt(bx.^2+by.^2+bz.^2));
    pp = (sx.^2+sy.^2+sz.^2)./rho;
    figure ; pcolor(x,y,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1])
    title(strcat('z = ',num2str(z(zzz)),' v_\perp'))

    %%%field aligned velocity
    %%%(v . B)/|B|
%     pp = (sx.*bx + sy.*by + sz.*bz)./(rho.*sqrt(bx.^2+by.^2+bz.^2));
%     figure ; pcolor(x,y,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1])   
%     title(strcat('z = ',num2str(z(zzz)),' v_{||}'))
% 
    %%%parallel current
    %pp = j_par;
    pp = log10(jx.^2+jy.^2+jz.^2);
    figure ; pcolor(x,y,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1])  
    title(strcat('z = ',num2str(z(zzz)),' J_{||}'))
% 
%     %%%perpendicular current
%     pp = j_perp;
%     figure ; pcolor(x,y,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1])  
%     title(strcat('z = ',num2str(z(zzz)),' J_\perp'))
% 
%     specific entropy
%     pp = p./(rho).^(5/3);
%     figure ; pcolor(x,y,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1]) 
%     title(strcat('z = ',num2str(z(zzz)),' specific entropy'))
% 
%     amount that angle wrt z has changed
%     bz_n0 = b0z./(sqrt(b0x.^2+b0y.^2+b0z.^2));
%     a0 = acos(bz_n0);
%     bz_n = bz./(sqrt(bx.^2+by.^2+bz.^2));
%     pp = a0 - acos(bz_n);
%     figure ; pcolor(x,y,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1]) 
%     title(strcat('z = ',num2str(z(zzz)),' misalignment'))
% 
%     all 3 components of field
%     pp = bz;
%     figure ; pcolor(x,y,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1])
%     hold on
%     ree = 10;
%     quiver(x(1:ree:end),y(1:ree:end),bx(1:ree:end,1:ree:end,zzz),by(1:ree:end,1:ree:end,zzz),'r') 
%     title(strcat('z = ',num2str(z(zzz)),' magnetic field'))
% 
%     plasma beta
%     pp = p./(bx.^2+by.^2+bz.^2);
%     figure ; pcolor(x,y,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1])
%     title(strcat('z = ',num2str(z(zzz)),' beta'))
% 
%     density
%     pp = rho;
%     figure ; pcolor(x,y,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1])
%     title(strcat('z = ',num2str(z(zzz)),' density'))
end

%magnetic perturbation for conjugate footpoint top boundary
% figure ; pcolor(xpo,ypo,dbz) ; shading interp ; colorbar ; hold on ; colormap('HSV')
% quiver(xpo(1:ree:end),ypo(1:ree:end),dbx(1:ree:end,1:ree:end),dby(1:ree:end,1:ree:end),'color',[0 0 0])
% title(strcat('z = ',num2str(zp1),' magnetic perturbation conjugate footpoint')) ; daspect([1 1 1])
% 
% %field aligned E cross B_perturbation for conjugate footpoint top boundary
% figure ; pcolor(xpo,ypo,poynt) ; shading interp ; colorbar ; hold on
% title(strcat('z = ',num2str(zp1),' field aligned E x B conjugate footpoint')) ; daspect([1 1 1])

%ave_all_fucking_figs ; close all ; clear all