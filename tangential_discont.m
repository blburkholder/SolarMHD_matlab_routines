%define 2 vertical lines
x1 = 62;
y1 = 82;
x2 = 62;
y2 = 71;

x11 = 66;
y11 = 82;
x22 = 66;
y22 = 71;

% x1 = 21;
% y1 = 29;
% x2 = 21;
% y2 = 18;
% 
% x11 = 25;
% y11 = 29;
% x22 = 25;
% y22 = 18;
% 
%[jx,jy,jz] = get_j(nx,ny,nz,bx,by,bz,b0x,b0y,b0z,difx,dify,difz);
%[e_par,j_par,j_perp,ex,ey,ez] = get_j2(nx,ny,nz,res,jx,jy,jz,bx,by,bz,sx./rho,sy./rho,sz./rho);

nnz = 130;
xx1 = x1*ones(length(y2:0.01:y1),1);
yy1 = y2:0.01:y1;
xx2 = x11*ones(length(y22:0.01:y11),1);
yy2 = y22:0.01:y11;
jj1 = interp2(x,y,j_par(:,:,nnz),xx1,yy1');
jj2 = interp2(x,y,j_par(:,:,nnz),xx2,yy2');
max_max1 = abs(jj1) == max(abs(jj1));
x_max1 = xx1(max_max1);
y_max1 = yy1(max_max1);
j_max1 = jj1(max_max1);
max_max2 = abs(jj2) == max(abs(jj2));
x_max2 = xx2(max_max2);
y_max2 = yy2(max_max2);
j_max2 = jj2(max_max2);
mmm = x_max1:0.01:x_max2;
mslope = (y_max1 - y_max2)/(x_max1 - x_max2); 
lll = mslope*((x_max1:0.1:x_max2) - x_max1) + y_max1;
ppp = (-1/mslope)*((mmm) - x_max1) + y2;

binterped_x = interp3(x,y,z,bx,mmm,ppp,z(nnz)*ones(size(mmm)));
binterped_y = interp3(x,y,z,by,mmm,ppp,z(nnz)*ones(size(mmm)));
binterped_z = interp3(x,y,z,bz,mmm,ppp,z(nnz)*ones(size(mmm)));
jinterped_x = interp3(x,y,z,jx,mmm,ppp,z(nnz)*ones(size(mmm)));
jinterped_y = interp3(x,y,z,jy,mmm,ppp,z(nnz)*ones(size(mmm)));
jinterped_z = interp3(x,y,z,jz,mmm,ppp,z(nnz)*ones(size(mmm)));
vinterped_x = interp3(x,y,z,sx./rho,mmm,ppp,z(nnz)*ones(size(mmm)));
vinterped_y = interp3(x,y,z,sy./rho,mmm,ppp,z(nnz)*ones(size(mmm)));
vinterped_z = interp3(x,y,z,sz./rho,mmm,ppp,z(nnz)*ones(size(mmm)));
rho_interped = interp3(x,y,z,rho,mmm,ppp,z(nnz)*ones(size(mmm)));

bis = sqrt(binterped_x.^2+binterped_y.^2+binterped_z.^2);
jis = sqrt(jinterped_x.^2+jinterped_y.^2+jinterped_z.^2);
vis = sqrt(vinterped_x.^2+vinterped_y.^2+vinterped_z.^2);
jdb = (jinterped_x.*binterped_x+jinterped_y.*binterped_y+jinterped_z.*binterped_z)./(bis.*jis);
jdbs = abs(jdb(1:end-1) - jdb(2:end));
%index max j dot b slope
imjdbs = find(jdbs == max(jdbs))+1;
bl = 20;

u_jdb = imjdbs+bl;
d_jdb = imjdbs-bl;

bxs = binterped_x(d_jdb:u_jdb);
bys = binterped_y(d_jdb:u_jdb);
bzs = binterped_z(d_jdb:u_jdb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%first find DHT frame for E
% avg_K = zeros(3,3);
% avg_KV = zeros(3,1);
% for i = 1:length(bxs)
%     Ki =  [binterped_y(i)^2+binterped_z(i)^2,-binterped_x(i)*binterped_y(i),-binterped_x(i)*binterped_z(i);...
%         -binterped_x(i)*binterped_y(i),binterped_x(i)^2+binterped_z(i)^2,-binterped_y(i)*binterped_z(i);...
%         -binterped_x(i)*binterped_z(i),-binterped_y(i)*binterped_z(i),binterped_x(i)^2+binterped_y(i)^2];
%     avg_K = avg_K + Ki;
%     avg_KV = avg_KV+Ki*[vinterped_x(i);vinterped_y(i);vinterped_z(i)];
% end
% 
% avg_K = avg_K/length(bxs);
% avg_KV = avg_KV/length(bxs);
% V_DHT = avg_KV\avg_K;
% vx = -sx./rho + V_DHT(1);
% vy = -sy./rho + V_DHT(2);
% vz = -sz./rho + V_DHT(3);
% 
% [e_par,j_par,j_perp,ex,ey,ez] = get_j2(nx,ny,nz,res,jx,jy,jz,bx,by,bz,vx,vy,vz);
% edht_interped_x = interp3(x,y,z,ex,mmm,ppp,z(nnz)*ones(size(mmm)));
% edht_interped_y = interp3(x,y,z,ey,mmm,ppp,z(nnz)*ones(size(mmm)));
% edht_interped_z = interp3(x,y,z,ez,mmm,ppp,z(nnz)*ones(size(mmm)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

 figure
%B cross product method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bu_x = mean(binterped_x(imjdbs:u_jdb));
Bu_y = mean(binterped_y(imjdbs:u_jdb));
Bu_z = mean(binterped_z(imjdbs:u_jdb));
Bd_x = mean(binterped_x(d_jdb:imjdbs));
Bd_y = mean(binterped_y(d_jdb:imjdbs));
Bd_z = mean(binterped_z(d_jdb:imjdbs));

BuBd =  cross([Bu_x,Bu_y,Bu_z],[Bd_x,Bd_y,Bd_z]); 
n_hat = BuBd/norm(BuBd);
plot3(2*[0,n_hat(1)],2*[0,n_hat(2)],2*[0,n_hat(3)],'r--')
hold on
sphere(20)
daspect([1 1 1])
dirdir = zeros(3,length(binterped_x));
% dirdir(1,:) = n_hat(1);
% dirdir(2,:) = n_hat(2);
% dirdir(3,:) = n_hat(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%E cross product method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Eu_x = mean(edht_interped_x(imjdbs:u_jdb));
Eu_y = mean(edht_interped_y(imjdbs:u_jdb));
Eu_z = mean(edht_interped_z(imjdbs:u_jdb));
Ed_x = mean(edht_interped_x(d_jdb:imjdbs));
Ed_y = mean(edht_interped_y(d_jdb:imjdbs));
Ed_z = mean(edht_interped_z(d_jdb:imjdbs));

Eu_m_Ed =  [Eu_x;Eu_y;Eu_z] - [Ed_x;Ed_y;Ed_z]; 
n_hat = Eu_m_Ed'/norm(Eu_m_Ed);
dirdir = zeros(3,length(binterped_x));
plot3(2*[0,-n_hat(1)],2*[0,-n_hat(2)],2*[0,n_hat(3)],'b--')
% dirdir(1,:) = n_hat(1);
% dirdir(2,:) = n_hat(2);
% dirdir(3,:) = n_hat(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MVA (minimum variance B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bx_m = mean(bxs);
by_m = mean(bys);
bz_m = mean(bzs);
bx_m2 = bx_m^2;
bx_2m = mean(bxs.^2);
by_m2 = by_m^2;
by_2m = mean(bys.^2);
bz_m2 = bz_m^2;
bz_2m = mean(bzs.^2);
bxz_m = mean(bxs.*bzs);
bxy_m = mean(bxs.*bys);
byz_m = mean(bys.*bzs);

MB = [ bx_2m - bx_m2, bxy_m-bx_m*by_m, bxz_m-bx_m*bz_m;...
    bxy_m-bx_m*by_m, by_2m-by_m2, byz_m-by_m*bz_m;...
    bxz_m-bx_m*bz_m, byz_m-by_m*bz_m, bz_2m-bz_m2];
[v,L] = eig(MB);
n_hat = v(:,1)';
plot3(2*[0,n_hat(1)],2*[0,n_hat(2)],2*[0,n_hat(3)],'r')
dirdir(1,:) = n_hat(1);
dirdir(2,:) = n_hat(2);
dirdir(3,:) = n_hat(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MVA (maximum variance E)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Exs = edht_interped_x(d_jdb:u_jdb);
Eys = edht_interped_y(d_jdb:u_jdb);
Ezs = edht_interped_z(d_jdb:u_jdb);
Ex_m = mean(Exs);
Ey_m = mean(Eys);
Ez_m = mean(Ezs);

Ex_m2 = Ex_m^2;
Ex_2m = mean(Exs.^2);
Ey_m2 = Ey_m^2;
Ey_2m = mean(Eys.^2);
Ez_m2 = Ez_m^2;
Ez_2m = mean(Ezs.^2);
Exz_m = mean(Exs.*Ezs);
Exy_m = mean(Exs.*Eys);
Eyz_m = mean(Eys.*Ezs);

MB = [ Ex_2m - Ex_m2, Exy_m-Ex_m*Ey_m, Exz_m-Ex_m*Ez_m;...
    Exy_m-Ex_m*Ey_m, Ey_2m-Ey_m2, Eyz_m-Ey_m*Ez_m;...
    Exz_m-Ex_m*Ez_m, Eyz_m-Ey_m*Ez_m, Ez_2m-Ez_m2];
[v,L] = eig(MB);
n_hat = v(:,3)';
plot3(2*[0,n_hat(1)],2*[0,n_hat(2)],2*[0,n_hat(3)],'b')
% dirdir(1,:) = n_hat(1);
% dirdir(2,:) = n_hat(2);
% dirdir(3,:) = n_hat(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
pcolor(x,y,j_par(:,:,nnz))
shading interp
daspect([1 1 1])
hold on 
plot(xx1,yy1,'k')
plot(x_max1,y_max1,'rp')
plot(xx2,yy2,'k')
plot(x_max2,y_max2,'rp')
%line connecting the maxima
plot(x_max1:0.1:x_max2,lll,'g')
plot(mmm,ppp,'g')
%plot([0,-10*dirdir(1,1)],[0,-10*dirdir(2,1)])
title(['time=',num2str(time),' z=',num2str(z(nnz))])

figure
subplot(3,1,1)
%quiver(mmm,0,binterped_x,binterped_y,'k')
plot(mmm,sqrt(sum(cross(dirdir,[binterped_x;binterped_y;binterped_z]).^2))./bis,'k')
hold on
plot(mmm,abs(dirdir(1)*binterped_x+dirdir(2)*binterped_y+dirdir(3)*binterped_z)./bis,'k--')
plot(mmm,bis,'k.')
plot(mmm,jdb,'k-.')
legend('|B_{||}|/|B|','|B_{\perp}|/|B|','|B|','J\cdot B')
plot([mmm(find(jdbs == max(jdbs))),mmm(find(jdbs == max(jdbs)))],[-1,1])
plot([mmm(find(jdbs == max(jdbs))+bl),mmm(find(jdbs == max(jdbs))+bl)],[-1,1])
plot([mmm(find(jdbs == max(jdbs))-bl),mmm(find(jdbs == max(jdbs))-bl)],[-1,1])
%axis([31.5 35.5 -1.1 1.1])
%axis([35.5 39.5 -1 1.5])
set(gca,'xticklabel',[])

subplot(3,1,2)
%quiver(mmm,0,jinterped_x,jinterped_y,'k')
hold on
plot(mmm,sqrt(sum(cross(dirdir,[jinterped_x;jinterped_y;jinterped_z]).^2))./jis,'k')
plot(mmm,abs(dirdir(1)*jinterped_x+dirdir(2)*jinterped_y+dirdir(3)*jinterped_z)./jis,'k--')
plot(mmm,jis,'k.')
legend('|j_{||}|/|j|','|j_{\perp}|/|j|','|j|')
%axis([31.5 35.5 -0.1 1.5])
%axis([35.5 39.5 0 1])
set(gca,'xticklabel',[])

subplot(3,1,3)
%quiver(mmm,-0.3,vinterped_x,vinterped_y,'k')
hold on
plot(mmm,sqrt(sum(cross(dirdir,[vinterped_x;vinterped_y;vinterped_z]).^2))./vis,'k')
plot(mmm,abs(dirdir(1)*vinterped_x+dirdir(2)*vinterped_y+dirdir(3)*vinterped_z)./vis,'k--')
plot(mmm,vis,'k.')
plot(mmm,(vinterped_x.*jinterped_x+vinterped_y.*jinterped_y+vinterped_z.*jinterped_z)./(vis.*jis),'k-.')
legend('|v_{||}|/|v|','|v_{\perp}|/|v|','|v|','v\cdot j')
%axis([31.5 35.5 0 1.2])
%axis([35.5 39.5 -1 1.5])
set(gca,'xticklabel',[])

vxb = vinterped_x(1);
vyb = vinterped_y(1);
vzb = vinterped_z(1);
VAxb = abs(binterped_x(1))/sqrt(rho_interped(1));
VAyb = abs(binterped_y(1))/sqrt(rho_interped(1));
VAzb = abs(binterped_z(1))/sqrt(rho_interped(1));
VAb = bis(1)/sqrt(rho_interped(1));

dvx = vinterped_x - vxb;
dvy = vinterped_y - vyb;
dvz = vinterped_z - vzb;
dv = vis-sqrt(vxb^2+vyb^2+vzb^2);
dvax = abs(binterped_x)./sqrt(rho_interped) - VAxb;
dvay = abs(binterped_y)./sqrt(rho_interped) - VAyb;
dvaz = abs(binterped_z)./sqrt(rho_interped) - VAzb;
dva = bis./sqrt(rho_interped) - VAb;
%dv = sqrt(dvx.^2+dvy.^2+dvz.^2);
%dva = sqrt(dvax.^2+dvay.^2+dvaz.^2);


figure
plot(mmm,dvx,'r')
hold on
plot(mmm,dvy,'g')
plot(mmm,dvz,'b')
plot(mmm,dvax,'r--')
plot(mmm,dvay,'g--')
plot(mmm,dvaz,'b--')
plot(mmm,dva,'k.-')
plot(mmm,dv,'k.-')
legend('\Delta v_x','\Delta v_y','\Delta v_z','\Delta v_{Ax}','\Delta v_{Ay}','\Delta v_{Az}','\Delta |v_A|','\Delta |v|')

n_hat = [bx_m,by_m,bz_m]/sqrt(bx_m^2+by_m^2+bz_m^2);
dirdir(1,:) = n_hat(1);
dirdir(2,:) = n_hat(2);
dirdir(3,:) = n_hat(3);

figure
plot(mmm,(dirdir(1,:).*dvx+dirdir(2,:).*dvy+dirdir(3,:).*dvz))
hold on
plot(mmm,(dirdir(1,:).*dvax+dirdir(2,:).*dvay+dirdir(3,:).*dvaz))
legend('Proj_B \Delta v','Proj_B \Delta v_a')





