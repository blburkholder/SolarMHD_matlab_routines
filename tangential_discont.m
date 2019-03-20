%define 2 vertical lines
%mhd10 3.9
tx1 = 34;
ty1 = 79;
tx2 = 34;
ty2 = 71;

tx11 = 39;
ty11 = 79;
tx22 = 39;
ty22 = 71;

% tx1 = 34;
% ty1 = 86;
% tx2 = 34;
% ty2 = 75;
% 
% tx11 = 39;
% ty11 = 86;
% tx22 = 39;
% ty22 = 75;

% tx1 = 67;
% ty1 = 30;
% tx2 = 67;
% ty2 = 16;
% 
% tx11 = 69;
% ty11 = 30;
% tx22 = 69;
% ty22 = 16;

%%%%%%%calculate current density and fac density for simulation volume
%  [jx,jy,jz] = get_j(nx,ny,nz,bx,by,bz,b0x,b0y,b0z,difx,dify,difz);
%  [e_par,j_par,j_perp,ex,ey,ez] = get_j2(nx,ny,nz,res,jx,jy,jz,bx,by,bz,...
%     sx./rho,sy./rho,sz./rho);
% j_par = j_par.^2+j_perp.^2;

%find maximum current density along vertical lines
%nnz = 75;
%nnz = 55;
nnz = 25;
txx1 = tx1*ones(length(ty2:0.01:ty1),1);
tyy1 = ty2:0.01:ty1;
txx2 = tx11*ones(length(ty22:0.01:ty11),1);
tyy2 = ty22:0.01:ty11;
tjj1 = interp2(x,y,j_par(:,:,nnz),txx1,tyy1');
tjj2 = interp2(x,y,j_par(:,:,nnz),txx2,tyy2');
tmax_max1 = abs(tjj1) == max(abs(tjj1));
tx_max1 = txx1(tmax_max1);
ty_max1 = tyy1(tmax_max1);
tj_max1 = tjj1(tmax_max1);
tmax_max2 = abs(tjj2) == max(abs(tjj2));
tx_max2 = txx2(tmax_max2);
ty_max2 = tyy2(tmax_max2);
tj_max2 = tjj2(tmax_max2);

%find line connecting maximum current along 2 vertical lines to define a
%normal to the boundary in the 2D x-y plane 
%mmm = (tx_max1-3):0.01:(tx_max2+2);
mmm = tx_max1:0.01:tx_max2;
mslope1 = (ty_max1 - ty_max2)/(tx_max1 - tx_max2); 
lll = mslope1*((mmm) - tx_max1) + ty_max1;
ppp = (-1/mslope1)*((mmm) - tx_max1) + ty1;
%ppp = (-1/mslope1)*((mmm)- (tx_max2+tx_max1)/2) + (ty_max1+ty_max2)/2; 

diffx = mmm(end)-mmm(1);
diffy = ppp(end)-ppp(1);
magg = sqrt(diffx^2+diffy^2);
dirx = diffx/magg;
diry = diffy/magg;

%make plot showing direction in 2D plane at given height
[xxx,yyy] = meshgrid(x,y);

figure
%pcolor(x,y(190:end),j_par(190:end,:,nnz))
pcolor(x,y,j_par(:,:,nnz))
shading interp
daspect([1 1 1])
hold on 
%h = streamslice(x,y,bx(:,:,nnz),by(:,:,nnz),3);
%set(h,'Color','k')
 plot(txx1,tyy1,'k')
 plot(tx_max1,ty_max1,'rp')
 plot(txx2,tyy2,'k')
 plot(tx_max2,ty_max2,'rp')
%line connecting the maxima
 plot(mmm,lll,'g')
 plot(mmm,ppp,'g')
title(['time=',num2str(time),' z=',num2str(z(nnz))])
colorbar

% figure
% pcolor(x,y,bx(:,:,nnz)); hold on
% plot(mmm,ppp,'k')
% shading interp
% colormap(jet)
% daspect([1 1 1])
% title('B_x')
% colorbar
% 
% figure
% pcolor(x,y,by(:,:,nnz)); hold on
% plot(mmm,ppp,'k')
% shading interp
% colormap(jet)
% daspect([1 1 1])
% colorbar
% title('B_y')
% 
% figure
% pcolor(x,y,bz(:,:,nnz)); hold on
% plot(mmm,ppp,'k')
% shading interp
% colormap(jet)
% daspect([1 1 1])
% title('B_z')
% colorbar

%find normal to boundary in 3d geometrically?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_y0 = ty2*mslope1 + tx_max1;
x_yf = -(y(end-1)-ty2)*mslope1+tx_max1;

mm = x_yf:0.01:x_y0;
pp = (-1/mslope1)*((mm) - tx_max1) + ty2;

xx = zeros(length(mm),153);
yy = zeros(length(mm),153);
zz = zeros(length(mm),153);
for i = 1:153
    xx(:,i) = mm;
    yy(:,i) = pp;
    zz(:,i) = z(i)*ones(size(mm));
end
%plot(mm,pp)
ang_slice = interp3(x,y,z,j_par,xx,yy,zz);
bzang_slice = interp3(x,y,z,bz,xx,yy,zz);

%define 2 horizontal lines
y1 = 72;
z1 = z(nnz)+1.5;
y2 = 77;
z2 = z(nnz)+1.5;

y11 = 72;
z11 = z(nnz)-1.5;
y22 = 77;
z22 = z(nnz)-1.5;

% y1 = 15;
% z1 = z(nnz)+2.5;
% y2 = 20;
% z2 = z(nnz)+2.5;
% 
% y11 = 15;
% z11 = z(nnz)-2.5;
% y22 = 20;
% z22 = z(nnz)-2.5;

zz1 = z1*ones(length(y1:0.01:y2),1);
yy1 = y1:0.01:y2;
zz2 = z11*ones(length(y11:0.01:y22),1);
yy2 = y11:0.01:y22;
jj1 = interp2(pp,z,ang_slice',yy1',zz1);
jj2 = interp2(pp,z,ang_slice',yy2',zz2);
max_max1 = jj1 == max(jj1);
z_max1 = zz1(max_max1);
y_max1 = yy1(max_max1);
j_max1 = jj1(max_max1);
max_max2 = jj2 == max(jj2);
z_max2 = zz2(max_max2);
y_max2 = yy2(max_max2);
j_max2 = jj2(max_max2);

nnn = z_max2:(abs(z_max2 - z_max1)/(length(mmm)-1)):z_max1;
mslope = (y_max1 - y_max2)/(z_max1 - z_max2); 
iii = mslope*((z_max2:0.1:z_max1) - z_max1) + y_max1;
qqq = (-1/mslope)*((nnn) - z_max1) + y1;

figure
pcolor(pp,z(20:end),ang_slice(:,20:end)')
%pcolor(pp,z(60:end),ang_slice(:,60:end)')
shading interp
daspect([1 1 1])
hold on
reso = 10;
qz = z(1:65);
qu = 0*bzang_slice;
qv = bzang_slice;
%h = quiver(pp(1:2*reso:end),qz(1:reso:end),qu(1:2*reso:end,1:reso:end)',qv(1:2*reso:end,1:reso:end)');
%set(h,'Color','k')
plot(yy1,zz1,'k')
plot(y_max1,z_max1,'rp')
plot(yy2,zz2,'k')
plot(y_max2,z_max2,'rp')
%line connecting the maxima
%plot(iii,10:0.1:15,'g','LineWidth',0.1)
plot(iii,z_max2:0.1:z_max1,'g')
plot(qqq,nnn,'g','LineWidth',2)

risetot = sqrt((mmm(1)-mmm(end))^2+(ppp(1)-ppp(end))^2)*mslope;
z_one = z(nnz) - risetot/2;
z_two = z(nnz) + risetot/2;
nnn = z_one:-(z_one-z_two)/(length(mmm)-1):z_two;

diffz = nnn(end)-nnn(1);
diffyy = qqq(end)-qqq(1);
maggg = sqrt(diffz^2+diffyy^2);
dirz = diffz/maggg;
frac = diffyy/maggg;
dirx = frac*dirx;
diry = frac*diry;
dirdir_g = [dirx,diry,dirz];
%dirdir_g = [dirx,diry,0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Methods from dissertationdraft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%interpolated quantities along 3d line

binterped_x = interp3(x,y,z,bx,mmm,ppp,nnn);
binterped_y = interp3(x,y,z,by,mmm,ppp,nnn);
binterped_z = interp3(x,y,z,bz,mmm,ppp,nnn);
jinterped_x = interp3(x,y,z,jx,mmm,ppp,nnn);
jinterped_y = interp3(x,y,z,jy,mmm,ppp,nnn);
jinterped_z = interp3(x,y,z,jz,mmm,ppp,nnn);
jpar_interped = interp3(x,y,z,j_par,mmm,ppp,nnn);
jperp_interped = interp3(x,y,z,j_perp,mmm,ppp,nnn);
vinterped_x = interp3(x,y,z,sx./rho,mmm,ppp,nnn);
vinterped_y = interp3(x,y,z,sy./rho,mmm,ppp,nnn);
vinterped_z = interp3(x,y,z,sz./rho,mmm,ppp,nnn);
rho_interped = interp3(x,y,z,rho,mmm,ppp,nnn);
p_interped = interp3(x,y,z,p,mmm,ppp,nnn);
epar_interped = interp3(x,y,z,e_par,mmm,ppp,nnn);

bis = sqrt(binterped_x.^2+binterped_y.^2+binterped_z.^2);
jis = sqrt(jinterped_x.^2+jinterped_y.^2+jinterped_z.^2);
vis = sqrt(vinterped_x.^2+vinterped_y.^2+vinterped_z.^2);

jdb = (jinterped_x.*binterped_x+jinterped_y.*binterped_y+...
    jinterped_z.*binterped_z)./(bis.*jis);

%index max j
imjdbs = find(abs(jpar_interped) == max(abs(jpar_interped)));
%%grid cells away from "middle" to perform averages
bl = 150;
%bl = 300;

u_jdb = imjdbs+bl;
d_jdb = imjdbs-bl;

bxs = binterped_x(d_jdb:u_jdb);
bys = binterped_y(d_jdb:u_jdb);
bzs = binterped_z(d_jdb:u_jdb);
vxs = vinterped_x(d_jdb:u_jdb);
vys = vinterped_y(d_jdb:u_jdb);
vzs = vinterped_z(d_jdb:u_jdb);

[rb1,rb2,rv1,rv2] = riemann_problemator(bxs(1),bxs(end),bys(1),bys(end),bzs(1),bzs(end),vxs(1),vys(1),vzs(1),vxs(end),vys(end),vzs(end))

% bxs(1)
% bxs(end)
% bys(1)
% bys(end)
% bzs(1)
% bzs(end)
% vxs(1)
% vys(1)
% vzs(1)
% vxs(end)
% vys(end)
% vzs(end)

%do the hard work
% V_DHT = get_DHT_frame(bxs,bys,bzs,vxs,vys,vzs);
% vx = -sx./rho + V_DHT(1);
% vy = -sy./rho + V_DHT(2);
% vz = -sz./rho + V_DHT(3);
% [e_par,j_par,j_perp,ex,ey,ez] = get_j2(nx,ny,nz,res,jx,jy,jz,bx,by,bz,vx,vy,vz);
% edht_interped_x = interp3(x,y,z,ex,mmm,ppp,nnn);
% edht_interped_y = interp3(x,y,z,ey,mmm,ppp,nnn);
% edht_interped_z = interp3(x,y,z,ez,mmm,ppp,nnn);
% exs = edht_interped_x(d_jdb:u_jdb);
% eys = edht_interped_y(d_jdb:u_jdb);
% ezs = edht_interped_z(d_jdb:u_jdb);
% [dir_cp_B,dir_cp_E] = normal_dir_cp(binterped_x,binterped_y,binterped_z,...
%     edht_interped_x,edht_interped_y,edht_interped_z,imjdbs,u_jdb,d_jdb);
% [dir_MVA_B,dir_MVA_E] = normal_dir_var(bxs,bys,bzs,exs,eys,ezs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%compare directions
% dirdir1 = zeros(3,length(binterped_x));
%     dirdir1(1,:) = dir_cp_B(1);
%     dirdir1(2,:) = dir_cp_B(2);
%     dirdir1(3,:) = dir_cp_B(3);
% dirdir2 = zeros(3,length(binterped_x)); 
%     dirdir2(1,:) = dir_cp_E(1);
%     dirdir2(2,:) = dir_cp_E(2);
%     dirdir2(3,:) = dir_cp_E(3);
% dirdir3 = zeros(3,length(binterped_x));
%     dirdir3(1,:) = dir_MVA_B(1);
%     dirdir3(2,:) = dir_MVA_B(2);
%     dirdir3(3,:) = dir_MVA_B(3);
% dirdir4 = zeros(3,length(binterped_x));
%     dirdir4(1,:) = dir_MVA_E(1);
%     dirdir4(2,:) = dir_MVA_E(2);
%     dirdir4(3,:) = dir_MVA_E(3);
dirdir5 = zeros(3,length(binterped_x));
    dirdir5(1,:) = dirdir_g(1);
    dirdir5(2,:) = dirdir_g(2);
    dirdir5(3,:) = dirdir_g(3);
dirdir = dirdir5;

% figure
% plot3(2*[0,dir_cp_B(1)],2*[0,dir_cp_B(2)],2*[0,dir_cp_B(3)],'r--')
% hold on
% plot3(2*[0,dir_cp_E(1)],2*[0,dir_cp_E(2)],2*[0,dir_cp_E(3)],'b--')
% plot3(2*[0,dir_MVA_B(1)],2*[0,dir_MVA_B(2)],2*[0,dir_MVA_B(3)],'r')
% plot3(2*[0,dir_MVA_E(1)],2*[0,dir_MVA_E(2)],2*[0,dir_MVA_E(3)],'b')
% plot3(2*[0,dirdir_g(1)],2*[0,dirdir_g(2)],2*[0,dirdir_g(3)],'k:')
% sphere(20)
% daspect([1 1 1])
% legend('cp B','cp E','MVA B','MVA E','geometry')

%%%%%%perform analysis
% 
%ppp = fliplr(ppp);
figure
subplot(6,1,1)
plot(ppp,binterped_x,'r')
hold on
plot(ppp,binterped_y,'g')
plot(ppp,binterped_z,'b')
legend('B_x','B_y','B_z')
plot(ppp,0*ppp,'k')
plot([ppp(imjdbs+bl) ppp(imjdbs+bl)],[-5 5],'k--')
plot([ppp(imjdbs) ppp(imjdbs)],[-5 5],'k--')
plot([ppp(imjdbs-bl) ppp(imjdbs-bl)],[-5 5],'k--')
set(gca,'xticklabel',[])
title('magnetic field')
axis([ppp(end) ppp(1) -5 5])

subplot(6,1,5); hold on
%plot(mmm,jpar_interped,'r')
plot(ppp,jinterped_x,'r')
plot(ppp,jinterped_y,'g')
plot(ppp,jinterped_z,'b')
hold on
%plot(mmm,jperp_interped,'g')
legend('j_x','j_y','j_z')
plot(ppp ,0*ppp,'k')
set(gca,'xticklabel',[])
title('current density')
axis([ppp(end) ppp(1) -1 1])
% 
% subplot(5,1,2)
% hold on
% plot(mmm(1:floor(length(binterped_x)/2)),jinterped_x(1:floor(length(binterped_x)/2)),'r')
% plot(mmm(1:floor(length(binterped_x)/2)),jinterped_y(1:floor(length(binterped_x)/2)),'g')
% plot(mmm(1:floor(length(binterped_x)/2)),jinterped_z(1:floor(length(binterped_x)/2)),'b')
% plot(mmm(1:floor(length(binterped_x)/2)),0*mmm(1:floor(length(binterped_x)/2)),'k')
% legend('j_x','j_y','j_z')
% set(gca,'xticklabel',[])
% title('current density')
% 
subplot(6,1,6)
hold on
plot(ppp,epar_interped,'r')
%plot(ppp,jperp_interped,'g')
plot(ppp,0*ppp,'k')
%legend('j_{||}','j_{\perp}')
%axis([31.5 35.5 -0.1 1.5])
%axis([35.5 39.5 0 1])
%set(gca,'xticklabel',[])
title('E_{||}')
axis([ppp(end) ppp(1) -0.1 0.1])

% subplot(4,1,3)
% hold on
% plot(mmm,sqrt(sum(cross(dirdir,[vinterped_x;vinterped_y;vinterped_z]).^2)),'r')
% plot(mmm,abs(dirdir(1)*vinterped_x+dirdir(2)*vinterped_y+dirdir(3)*vinterped_z),'g')
% plot(mmm,vis,'b')
% plot(mmm,0*mmm,'k')
% legend('|v_{||}|','v_{\perp}','|v|')
% %axis([31.5 35.5 0 1.2])
% %axis([35.5 39.5 -1 1.5])
% set(gca,'xticklabel',[])

% subplot(4,1,4)
% plot(mmm,(vinterped_x.*jinterped_x+vinterped_y.*jinterped_y+vinterped_z.*jinterped_z)./(vis.*jis),'r')
% hold on
% plot(mmm,jdb,'b')
% plot(mmm,0*mmm,'k')
% legend('v \cdot B','j \cdot B')

% subplot(5,1,4)
% hold on
% plot(mmm(1:floor(length(binterped_x)/2)),dirdir(1,1)*jinterped_x(1:floor(length(binterped_x)/2))+dirdir(2,1)*jinterped_y(1:floor(length(binterped_x)/2))+dirdir(3,1)*jinterped_z(1:floor(length(binterped_x)/2)),'r')
% plot(mmm(1:floor(length(binterped_x)/2)),sqrt(sum(cross(dirdir(:,1:floor(length(binterped_x)/2)),[jinterped_x(1:floor(length(binterped_x)/2));jinterped_y(1:floor(length(binterped_x)/2));jinterped_z(1:floor(length(binterped_x)/2))]).^2)),'g')
% plot(mmm(1:floor(length(binterped_x)/2)),0*mmm(1:floor(length(binterped_x)/2)),'k')
% legend('j_{||}','j_{\perp}')
% %axis([31.5 35.5 -0.1 1.5])
% %axis([35.5 39.5 0 1])
% set(gca,'xticklabel',[])
% title('j boundary normal coordinates')

%subplot(3,1,3)
% %quiver(mmm,0,binterped_x,binterped_y,'k')
% b_perp = sqrt(sum(cross(dirdir,[binterped_x;binterped_y;binterped_z]).^2));
% semilogy(mmm,abs(dirdir(1,1)*binterped_x+dirdir(2,1)*binterped_y+dirdir(3,1)*binterped_z),'r')
% hold on
% semilogy(mmm,b_perp,'g')
% legend('B_{||}','B_{\perp}')
% plot(mmm,0*mmm,'k')
% %axis([31.5 35.5 -1.1 1.1])
% %axis([35.5 39.5 -1 1.5])
% set(gca,'xticklabel',[])
% title('magnetic field boundary normal coordinates')

%figure
subplot(6,1,3)
hold on
plot(ppp,vinterped_x,'r')
plot(ppp,vinterped_y,'g')
plot(ppp,vinterped_z,'b')
plot(ppp,0*ppp,'k')
legend('v_x','v_y','v_z')
set(gca,'xticklabel',[])
title('velocity')
axis([ppp(end) ppp(1) -0.3 0.3])

subplot(6,1,2)
hold on
plot(ppp,rho_interped)
set(gca,'xticklabel',[])
title('density')
axis([ppp(end) ppp(1) 0 0.5])

subplot(6,1,4)
plot(ppp,p_interped./rho_interped.^(5/3))
title('entropy')
axis([ppp(end) ppp(1) 0 250])
set(gca,'xticklabel',[])
% 
% subplot(3,1,3)
% hold on
% plot(mmm,p_interped)
% title('pressure')

% subplot(3,1,3)
% hold on
% plot(mmm,binterped_x./sqrt(rho_interped))
% plot(mmm,binterped_y./sqrt(rho_interped))
% plot(mmm,binterped_z./sqrt(rho_interped))
% set(gca,'xticklabel',[])
% title('B/sqrt(\rho)')
[one,two] = simulation_walen(ppp,rho_interped,binterped_x,binterped_y,binterped_z,...
    vinterped_x,vinterped_y,vinterped_z,...
    binterped_x./sqrt(rho_interped),binterped_y./sqrt(rho_interped),binterped_z./sqrt(rho_interped),0.1,1.5)


