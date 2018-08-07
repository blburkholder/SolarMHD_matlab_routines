
xxfin = zeros(size(xfin1));
yyfin = zeros(size(yfin1));

for i = 1:length(xfin1)
    for j = 1:length(xfin1)
        if xfin1(i,j) > x(end-1)
            xxfin(i,j) = x(end-1) - mod(xfin1(i,j),x(end-1));
            yyfin(i,j) = y(end-1) - yfin1(i,j);
        elseif xfin1(i,j) < x(2)
            xxfin(i,j) = abs(xfin1(i,j));
            yyfin(i,j) = y(end-1) - yfin1(i,j);
        elseif yfin1(i,j) > y(end-1)
            yyfin(i,j) = y(end-1) - mod(yfin1(i,j),y(end-1));
            xxfin(i,j) = x(end-1) - xfin1(i,j);
        elseif yfin1(i,j) < y(2)
            yyfin(i,j) = abs(yfin1(i,j));
            xxfin(i,j) = x(end-1) - xfin1(i,j);
        else
            xxfin(i,j) = xfin1(i,j);
            yyfin(i,j) = yfin1(i,j);
        end
    end
end

    %ree = 2;
    %[ssx,ssy] = meshgrid(x(2:ree:end-1),y(2:ree:end-1));
    %flines = perio_stream3(x,y,z,bx,by,bz,ssx,ssy,z(end-1)*ones(size(ssx)));
    %flines = stream3(x,y,z,-bx,-by,-bz,ssx,ssy,z(end-1)*ones(size(ssx)));


%figure
%streamline(flines)
% 
%     xxfin = zeros(length(ssx)^2,1);
%     yyfin = zeros(length(ssx)^2,1);
%     zzfin = zeros(length(ssx)^2,1);

%     for i = 1:length(flines)
%         xxfin(i) = flines{i}(end,2);
%         yyfin(i) = flines{i}(end,1);
%         zzfin(i) = flines{i}(end,3);
%     end
% % 
% [jx,jy,jz] = get_j(nx,ny,nz,bx,by,bz,b0x,b0y,b0z,difx,dify,difz);
% [e_par,j_par,j_perp,ex,ey,ez] = get_j2(nx,ny,nz,res,jx,jy,jz,bx,by,bz,sx./rho,sy./rho,sz./rho);
% [dbx,dby,dbz,poynt] = magn_perturb(bx,by,bz,rho,xpo,ypo,xfin1,yfin1,zfin1,x,y,z,sx,sy,sz,ex,ey,ez);

%bbx = interp3(x,y,bx(:,:,2),xxfin,yyfin);
%bby = interp3(x,y,by(:,:,2),xxfin,yyfin);

vvx = interp3(x,y,z,sx./rho,xxfin,yyfin,zfin1);
vvy = interp3(x,y,z,sy./rho,xxfin,yyfin,zfin1);
vv = vvx.^2+vvy.^2;
%vv(vv > 0.01) = 0.01;

bbx = interp3(x,y,z,bx,xxfin,yyfin,zfin1);
bby = interp3(x,y,z,by,xxfin,yyfin,zfin1);
bbz = interp3(x,y,z,bz,xxfin,yyfin,zfin1);

figure
pcolor(xpo,ypo,bbz)
%pcolor(ssx,ssy,reshape(vvx.^2+vvy.^2,[length(ssx),length(ssx)]))
hold on
h = quiver(xpo(1:40:end,1:40:end),ypo(1:40:end,1:40:end),vvx(1:40:end,1:40:end),vvy(1:40:end,1:40:end));
%h = quiver(ssx,ssy,vvx,vvy);
set(h,'Color','k')
shading interp
daspect([1 1 1])
colorbar
legend({'$B_z$','$\overrightarrow{v}$'},'Interpreter','latex')
title('Photosphere conditions on open flux')

[jx,jy,jz] = get_j(nx,ny,nz,bx,by,bz,b0x,b0y,b0z,difx,dify,difz);
[e_par,j_par,j_perp,ex,ey,ez] = get_j2(nx,ny,nz,res,jx,jy,jz,bx,by,bz,sx./rho,sy./rho,sz./rho);

figure
pcolor(x,y,j_par(:,:,end-1))
shading interp
hold on
[xx,yy] = meshgrid(x,y);
h = quiver(xx(1:10:end,1:10:end),yy(1:10:end,1:10:end),sx(1:10:end,1:10:end,end-1)./rho(1:10:end,1:10:end,end-1),sy(1:10:end,1:10:end,end-1)./rho(1:10:end,1:10:end,end-1));
set(h,'Color','k')
colorbar
title('Local conditions on open flux')
legend({'$j_{||}$','$\overrightarrow{v}$'},'Interpreter','latex')
daspect([1 1 1])

figure
pcolor(xpo,ypo,dbx.^2+dby.^2)
hold on
h = quiver(xpo(1:40:end,1:40:end),ypo(1:40:end,1:40:end),dbx(1:40:end,1:40:end),dby(1:40:end,1:40:end));
set(h,'Color','k')
shading interp
daspect([1 1 1])
legend({'$\delta B_x^2 + \delta B_y^2$','$\delta \overrightarrow{B}_\perp$'},'Interpreter','latex')
title('magnetic perturbation at photosphere')

figure
pcolor(xpo,ypo,sqrt(bbx.^2+bby.^2+bzfin1.^2))
%pcolor(ssx,ssy,reshape(bbz,[size(ssx)]))
hold on
h1 = quiver(xpo(1:40:end,1:40:end),ypo(1:40:end,1:40:end),abs(bzfin1(1:40:end,1:40:end)),abs(bzfin1(1:40:end,1:40:end)));
set(h1,'Color','k')
h2 = quiver(xpo(1:40:end,1:40:end),ypo(1:40:end,1:40:end),bbx(1:40:end,1:40:end),bby(1:40:end,1:40:end));
set(h2,'Color','r')
shading interp
daspect([1 1 1])
legend({'$|B|$','$|B_z|$','$\overrightarrow{B}_\perp$'},'Interpreter','latex')
title('photospheric field')

figure
pcolor(xpo,ypo,poynt)
shading interp
daspect([1 1 1])
title('(E\times \delta B)\cdot B')

