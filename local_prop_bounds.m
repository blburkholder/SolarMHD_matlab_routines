prepare3(6)
mentrop_prepare2(6)
load('mhd6.mat')
load('mentrop6.mat')

%first load an mhd.mat and a mentrop.mat from same time in simulation
for zzz = [90,149]
    [jx,jy,jz] = get_j(nx,ny,nz,bx,by,bz,b0x,b0y,b0z,difx,dify,difz);
    [e_par,j_par,j_perp,ex,ey,ez] = get_j2(nx,ny,nz,res,jx,jy,jz,bx,by,bz,sx,sy,sz,rho);

    %perpendicular velocity
    pp = sqrt(sx.^2+sy.^2)./rho;
    figure ; pcolor(y,x,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1])
    title(strcat('z = ',num2str(z(zzz)),' v_\perp'))

    %field aligned velocity
    pp = sz./rho;
    figure ; pcolor(y,x,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1])   
    title(strcat('z = ',num2str(z(zzz)),' v_{||}'))

    %parallel current
    pp = j_par;
    figure ; pcolor(y,x,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1])  
    title(strcat('z = ',num2str(z(zzz)),' J_{||}'))

    %perpendicular current
    pp = j_perp;
    figure ; pcolor(y,x,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1])  
    title(strcat('z = ',num2str(z(zzz)),' J_\perp'))

    %specific entropy
    pp = p./(rho).^(5/3);
    figure ; pcolor(y,x,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1]) 
    title(strcat('z = ',num2str(z(zzz)),' specific entropy'))

    %amount that angle wrt z has changed
    bz_n0 = b0z./(sqrt(b0x.^2+b0y.^2+b0z.^2));
    a0 = acos(bz_n0);
    bz_n = bz./(sqrt(bx.^2+by.^2+bz.^2));
    pp = a0 - acos(bz_n);
    figure ; pcolor(y,x,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1]) 
    title(strcat('z = ',num2str(z(zzz)),' misalignment'))

    %all 3 components of field
    pp = bz;
    figure ; pcolor(y,x,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1])
    hold on
    ree = 10;
    quiver(y(1:ree:end),x(1:ree:end),by(1:ree:end,1:ree:end,zzz),bx(1:ree:end,1:ree:end,zzz),'r') 
    title(strcat('z = ',num2str(z(zzz)),' magnetic field'))

    %plasma beta
    pp = p./(bx.^2+by.^2+bz.^2);
    figure ; pcolor(y,x,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1])
    title(strcat('z = ',num2str(z(zzz)),' beta'))

    %density
    pp = rho;
    figure ; pcolor(y,x,pp(:,:,zzz)) ; shading interp ; colorbar ; daspect([1 1 1])
    title(strcat('z = ',num2str(z(zzz)),' density'))
end

%magnetic perturbation for conjugate footpoint top boundary
[dbx,dby,dbz,poynt] = magn_perturb(bx,by,bz,rho,xpo,ypo,zp1,xfin1,yfin1,zfin1,x,y,z,sx,sy,sz,ex,ey,ez);
ree = 20;
figure ; pcolor(ypo,xpo,dbz) ; shading interp ; colorbar ; hold on
quiver(ypo(1:ree:end),xpo(1:ree:end),dby(1:ree:end,1:ree:end),dbx(1:ree:end,1:ree:end))
title(strcat('z = ',num2str(zp1),' magnetic perturbation conjugate footpoint')) ; daspect([1 1 1])

%field aligned E cross B_perturbation for conjugate footpoint top boundary
figure ; pcolor(ypo,xpo,poynt) ; shading interp ; colorbar ; hold on
title(strcat('z = ',num2str(zp1),' field aligned E x B conjugate footpoint')) ; daspect([1 1 1])

save_all_fucking_figs ; close all ; clear all