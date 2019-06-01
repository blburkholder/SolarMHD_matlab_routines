tx1 = 21;
ty1 = 72;
tx2 = 21;
ty2 = 62;

tx11 = 23;
ty11 = 72;
tx22 = 23;
ty22 = 62;
txy = [tx1 ty1 tx2 ty2 tx11 ty11 tx22 ty22];

y1 = 60;
y2 = 70;
y11 = 60;
y22 = 70;
tyz = [y1 y2 y11 y22];

nnz = 45;
bl = 350;
rl = 5;

[~,~,~,~,~,~,~,mmm,ppp,nnn,dirdir] = tangential_discont(x,y,z,nx,ny,nz,bx,by,bz,b0x,b0y,b0z,difx,dify,difz,sx,sy,sz,res,rho,p,time,txy,tyz,nnz,bl,rl);
orang = [246/255 103/255 51/255];
purp = [82/255 45/255 128/255];

binterped_x = interp3(x,y,z,permute(bx,[2 1 3]),mmm,ppp,nnn);
binterped_y = interp3(x,y,z,permute(by,[2 1 3]),mmm,ppp,nnn);
binterped_z = interp3(x,y,z,permute(bz,[2 1 3]),mmm,ppp,nnn);
% jinterped_x = interp3(x,y,z,permute(jx,[2 1 3]),mmm,ppp,nnn);
% jinterped_y = interp3(x,y,z,permute(jy,[2 1 3]),mmm,ppp,nnn);
% jinterped_z = interp3(x,y,z,permute(jz,[2 1 3]),mmm,ppp,nnn);
vinterped_x = interp3(x,y,z,permute(sx./rho,[2 1 3]),mmm,ppp,nnn);
vinterped_y = interp3(x,y,z,permute(sy./rho,[2 1 3]),mmm,ppp,nnn);
vinterped_z = interp3(x,y,z,permute(sz./rho,[2 1 3]),mmm,ppp,nnn);
rho_interped = interp3(x,y,z,permute(rho,[2 1 3]),mmm,ppp,nnn);
p_interped = interp3(x,y,z,permute(p,[2 1 3]),mmm,ppp,nnn);

b_perp = sqrt(sum(cross(dirdir,[binterped_x;binterped_y;binterped_z]).^2));
bpar = abs(dirdir(1,1)*binterped_x+dirdir(2,1)*binterped_y+dirdir(3,1)*binterped_z);

figure
subplot(3,2,1); hold on
fill([ppp(1) ppp(bpar == min(bpar)) ppp(bpar == min(bpar)) ppp(1)],[-0.5 -0.5  2 2],orang)
fill([ppp(end) ppp(bpar == min(bpar)) ppp(bpar == min(bpar)) ppp(end)],[-0.5 -0.5  2 2],purp)
plot(ppp,0*ppp,'k')
p1 = plot(ppp,binterped_x,'c','LineWidth',3);
p2 = plot(ppp,binterped_y,'g','LineWidth',3);
p3 = plot(ppp,binterped_z,'b','LineWidth',3);
p4 = plot(ppp,p_interped./(binterped_x.^2+binterped_y.^2+binterped_z.^2),'k:','LineWidth',3);
legend([p1,p2,p3,p4],'B_x','B_y','B_z','\beta','Location','southwest')
% plot([ppp(imjdbs+bl) ppp(imjdbs+bl)],[-5 5],'k--')
% plot([ppp(imjdbs) ppp(imjdbs)],[-5 5],'k--')
% plot([ppp(imjdbs-bl) ppp(imjdbs-bl)],[-5 5],'k--')
set(gca,'xticklabel',[])
title('magnetic field')
axis([ppp(end) ppp(1) -0.5 2])
set(gca,'FontSize',16)

subplot(3,2,2); hold on
fill([ppp(1) ppp(bpar == min(bpar)) ppp(bpar == min(bpar)) ppp(1)],[.0001 .0001  2.5 2.5],orang)
fill([ppp(end) ppp(bpar == min(bpar)) ppp(bpar == min(bpar)) ppp(end)],[.0001 .0001  2.5 2.5],purp)
p1 = semilogy(ppp,bpar,'c','LineWidth',3);
p2 = semilogy(ppp,b_perp,'g','LineWidth',3);
%set(gca,'YScale','log')
legend([p1,p2],'B_{n}','B_{t}')
plot(ppp,0*ppp,'k')
axis([ppp(end) ppp(1) 0.0001 2.5])
title('magnetic field boundary normal coordinates')
set(gca,'xticklabel',[])
set(gca,'FontSize',16)

% subplot(3,2,2)
% plot(ppp,abs(binterped_x./sqrt(rho_interped)),'r')
% hold on
% plot(ppp,abs(binterped_y./sqrt(rho_interped)),'g')
% plot(ppp,abs(binterped_z./sqrt(rho_interped)),'b')
% legend('v_{Ax}','v_{Ay}','v_{Az}','Location','northeast')
% plot(ppp,0*ppp,'k')
% set(gca,'xticklabel',[])
% title('Alfven velocity')
%axis([ppp(end) ppp(1) 0 10])

subplot(3,2,3)
hold on
fill([ppp(1) ppp(bpar == min(bpar)) ppp(bpar == min(bpar)) ppp(1)],[-.1 -.1  0.3 0.3],orang)
fill([ppp(end) ppp(bpar == min(bpar)) ppp(bpar == min(bpar)) ppp(end)],[-.1 -.1  0.3 0.3],purp)
p1 = plot(ppp,vinterped_x,'c','LineWidth',3);
p2 = plot(ppp,vinterped_y,'g','LineWidth',3);
p3 = plot(ppp,vinterped_z,'b','LineWidth',3);
plot(ppp,0*ppp,'k')
legend([p1,p2,p3],'v_x','v_y','v_z','Location','northeast')
set(gca,'xticklabel',[])
title('velocity')
axis([ppp(end) ppp(1) -0.1 0.3])
set(gca,'FontSize',16)

subplot(3,2,4)
hold on
fill([ppp(1) ppp(bpar == min(bpar)) ppp(bpar == min(bpar)) ppp(1)],[0 0  0.1 0.1],orang)
fill([ppp(end) ppp(bpar == min(bpar)) ppp(bpar == min(bpar)) ppp(end)],[0 0  0.1 0.1],purp)
plot(ppp,rho_interped,'k','LineWidth',3)
set(gca,'xticklabel',[])
title('density')
axis([ppp(end) ppp(1) 0 0.1])
set(gca,'FontSize',16)

subplot(3,2,5); hold on
fill([ppp(1) ppp(bpar == min(bpar)) ppp(bpar == min(bpar)) ppp(1)],[0 0  1.5 1.5],orang)
fill([ppp(end) ppp(bpar == min(bpar)) ppp(bpar == min(bpar)) ppp(end)],[0 0  1.5 1.5],purp)
plot(ppp,p_interped,'k','LineWidth',3)
%plot(ppp,binterped_x.^2+binterped_y.^2+binterped_z.^2,'k-.')
%plot(ppp,p_interped+binterped_x.^2+binterped_y.^2+binterped_z.^2,'k-')
title('thermal pressure')
%legend('thermal','B^2','total','Location','northeast')
axis([ppp(end) ppp(1) 0 1.5])
xlabel('x')
set(gca,'FontSize',16)

subplot(3,2,6); hold on
fill([ppp(1) ppp(bpar == min(bpar)) ppp(bpar == min(bpar)) ppp(1)],[60 60  110 110],orang)
fill([ppp(end) ppp(bpar == min(bpar)) ppp(bpar == min(bpar)) ppp(end)],[60 60 110 110],purp)
plot(ppp,p_interped./rho_interped.^(5/3),'k','LineWidth',3)
title('entropy')
axis([ppp(end) ppp(1) 60 110])
xlabel('x')
set(gca,'FontSize',16)
 
% figure
% subplot(3,1,2)
% hold on
% plot(ppp,epar_interped,'r')
% %plot(ppp,jperp_interped,'g')
% plot(ppp,0*ppp,'k')
% %legend('j_{||}','j_{\perp}')
% %axis([31.5 35.5 -0.1 1.5])
% %axis([35.5 39.5 0 1])
% %set(gca,'xticklabel',[])
% title('E_{||}')
% %axis([ppp(end) ppp(1) -0.01 0.01])
% set(gca,'xticklabel',[])

% subplot(3,1,1); hold on
% %plot(mmm,jpar_interped,'r')
% plot(ppp,jinterped_x,'r')
% plot(ppp,jinterped_y,'g')
% plot(ppp,jinterped_z,'b')
% hold on
% %plot(mmm,jperp_interped,'g')
% legend('j_x','j_y','j_z')
% plot(ppp ,0*ppp,'k')
% title('current density')
% %axis([ppp(end) ppp(1) -0.5 0.5])
% set(gca,'xticklabel',[])

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
%quiver(mmm,0,binterped_x,binterped_y,'k')
% b_perp = sqrt(sum(cross(dirdir,[binterped_x;binterped_y;binterped_z]).^2));
% semilogy(ppp,abs(dirdir(1,1)*binterped_x+dirdir(2,1)*binterped_y+dirdir(3,1)*binterped_z),'r')
% hold on
% semilogy(ppp,b_perp,'g')
% legend('B_{n}','B_{t}')
% plot(ppp,0*ppp,'k')
% %axis([ppp(end) ppp(1) 0.01 2.5])
% title('magnetic field boundary normal coordinates')
% xlabel('x')

% subplot(3,1,3)
% hold on
% plot(mmm,binterped_x./sqrt(rho_interped))
% plot(mmm,binterped_y./sqrt(rho_interped))
% plot(mmm,binterped_z./sqrt(rho_interped))
% set(gca,'xticklabel',[])
% title('B/sqrt(\rho)')
[one,two] = simulation_walen(ppp,binterped_x,binterped_y,binterped_z,...
    vinterped_x,vinterped_y,vinterped_z,...
    abs(binterped_x),abs(binterped_y),abs(binterped_z),0.5,1.5);
%    abs(binterped_x)./sqrt(rho_interped),abs(binterped_y)./sqrt(rho_interped),abs(binterped_z)./sqrt(rho_interped),0.5,3);
% subplot(3,2,6)
% plot(ppp,p_interped./rho_interped.^(5/3))
% title('entropy')
%axis([-2 2 0.9 1.4])

yu1 = 1;
yu2 = 800;
b1 = sqrt(binterped_x(yu1)^2+binterped_y(yu1)^2+binterped_z(yu1)^2);
b2 = sqrt(binterped_x(yu2)^2+binterped_y(yu2)^2+binterped_z(yu2)^2);
rho1 = rho_interped(yu1);
rho2 = rho_interped(yu2);
vout = sqrt(b1*b2*(b1+b2)/(rho1*b2+rho2*b1)/(4*pi));

E_rate = 0.1*b1*b2/(b1+b2)*vout
