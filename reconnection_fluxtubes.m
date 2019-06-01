%mhd10_3.9s3.mat
tx1 = 35;
ty1 = 65;
tx2 = 35;
ty2 = 55;

tx11 = 33;
ty11 = 65;
tx22 =33;
ty22 = 55;
txy = [tx1 ty1 tx2 ty2 tx11 ty11 tx22 ty22];

y1 = 55;
y2 = 65;
y11 = 55;
y22 = 65;
tyz = [y1 y2 y11 y22];

nnz = 45;
bl = 150;
rl = 5;
fart = 100:600;

%extendo-o-run mhd10
% tx1 = 39;
% ty1 = 20;
% tx2 = 39;
% ty2 = 10;
% 
% tx11 = 37;
% ty11 = 20;
% tx22 =37;
% ty22 = 10;
% txy = [tx1 ty1 tx2 ty2 tx11 ty11 tx22 ty22];
% 
% y1 = 10;
% y2 = 20;
% y11 = 10;
% y22 = 20;
% tyz = [y1 y2 y11 y22];
% 
% nnz = 35;
% bl = 150;
% rl = 5;
% fart = 1:length(binterped_x);

[brotx,broty,brotz,jx,jy,jz,epar,mmm,ppp,nnn,dirdir] =...
 tangential_discont(x,y,z,nx,ny,nz,bx,by,bz,b0x,b0y,b0z,difx,dify,difz,sx,sy,sz,res,rho,p,time,txy,tyz,nnz,bl,rl);
  %[jx,jy,jz] = get_j(nx,ny,nz,bx,by,bz,b0x,b0y,b0z,difx,dify,difz);
  %[e_par,j_par,j_perp,~,~,~] = get_j2(nx,ny,nz,res,jx,jy,jz,bx,by,bz,...
  %   sx./rho,sy./rho,sz./rho);
%epar = interp3(x,y,z,permute(e_par,[2 1 3]),mmm,ppp,nnn);

binterped_x = interp3(x,y,z,permute(bx,[2 1 3]),mmm,ppp,nnn);
binterped_y = interp3(x,y,z,permute(by,[2 1 3]),mmm,ppp,nnn);
binterped_z = interp3(x,y,z,permute(bz,[2 1 3]),mmm,ppp,nnn);
vinterped_x = interp3(x,y,z,permute(sx./rho,[2 1 3]),mmm,ppp,nnn);
vinterped_y = interp3(x,y,z,permute(sy./rho,[2 1 3]),mmm,ppp,nnn);
vinterped_z = interp3(x,y,z,permute(sz./rho,[2 1 3]),mmm,ppp,nnn);
rho_interped = interp3(x,y,z,permute(rho,[2 1 3]),mmm,ppp,nnn);
p_interped = interp3(x,y,z,permute(p,[2 1 3]),mmm,ppp,nnn);
%jinterped_x = interp3(x,y,z,permute(jx,[2 1 3]),mmm,ppp,nnn);
%jinterped_y = interp3(x,y,z,permute(jy,[2 1 3]),mmm,ppp,nnn);
%jinterped_z = interp3(x,y,z,permute(jz,[2 1 3]),mmm,ppp,nnn);

figure;
subplot(4,2,1); hold on
plot(ppp(fart),binterped_x(fart),'c','LineWidth',3)
plot(ppp(fart),binterped_y(fart),'g','LineWidth',3)
plot(ppp(fart),binterped_z(fart),'b','LineWidth',3)
plot(ppp(fart),p_interped(fart)./(binterped_x(fart).^2+binterped_y(fart).^2+binterped_z(fart).^2),'k:','LineWidth',3);
legend('B_x','B_y','B_z','\beta')
plot(ppp(fart),0*ppp(fart),'k')
title('magnetic field')
axis([ppp(fart(1)) ppp(fart(end)) -1 3])
set(gca,'xticklabel',[])

subplot(4,2,2); hold on
plot(ppp(fart),brotx(fart),'c','LineWidth',3)
plot(ppp(fart),broty(fart),'g','LineWidth',3)
plot(ppp(fart),brotz(fart),'b','LineWidth',3)
legend('B_{min}','B_{int}','B_{max}')
plot(ppp(fart),0*ppp(fart),'k')
title('B variance coordinates')
axis([ppp(fart(1)) ppp(fart(end)) -1 2])
set(gca,'xticklabel',[])

subplot(4,2,3); hold on
plot(ppp(fart),jx(fart),'c','LineWidth',3)
plot(ppp(fart),jy(fart),'g','LineWidth',3)
plot(ppp(fart),jz(fart),'b','LineWidth',3)
j = sqrt(jx(fart).^2+jy(fart).^2+jz(fart).^2);
plot(ppp(fart),j,'k:','LineWidth',3)
legend('j_x','j_y','j_z','|j|')
plot(ppp(fart),0*ppp(fart),'k')
title('current density')
axis([ppp(fart(1)) ppp(fart(end)) -0.4 0.4])
set(gca,'xticklabel',[])

subplot(4,2,4); hold on
plot(ppp(fart),epar(fart),'c','LineWidth',3)
plot(ppp(fart),0*ppp(fart),'k')
title('E_{||}')
axis([ppp(fart(1)) ppp(fart(end)) -.001 .004])
set(gca,'xticklabel',[])

subplot(4,2,5); hold on
plot(ppp(fart),vinterped_x(fart),'c','LineWidth',3)
plot(ppp(fart),vinterped_y(fart),'g','LineWidth',3)
plot(ppp(fart),vinterped_z(fart),'b','LineWidth',3)
legend('v_x','v_y','v_z')
plot(ppp(fart),0*ppp(fart),'k')
title('velocity')
axis([ppp(fart(1)) ppp(fart(end)) -0.1 0.1])
set(gca,'xticklabel',[])

subplot(4,2,6); hold on
plot(ppp(fart),rho_interped(fart),'b','LineWidth',3)
title('density')
axis([ppp(fart(1)) ppp(fart(end)) 0.05 0.15])
set(gca,'xticklabel',[])

subplot(4,2,7); hold on
plot(ppp(fart),p_interped(fart),'k','LineWidth',3)
title('thermal pressure')
axis([ppp(fart(1)) ppp(fart(end)) 0 3])
xlabel('x')

subplot(4,2,8); hold on
S = p_interped./rho_interped.^(5/3);
plot(ppp(fart),S(fart),'k','LineWidth',3)
title('entropy')
axis([ppp(fart(1)) ppp(fart(end)) 80 90])
%m1 = plot([ppp(j == max(j)) ppp(j == max(j))],[78 86],'r','LineWidth',2);
%m2 = plot([ppp(S == max(S)) ppp(S == max(S))],[78 86],'k','LineWidth',2);
%m3 = plot([ppp(epar == max(epar)) ppp(epar == max(epar))],[78 84],'c','LineWidth',2);
%m4 = plot([ppp(rho_interped == max(rho_interped)) ppp(rho_interped == max(rho_interped))],[78 84],'b','LineWidth',2);
%hleg = legend([m1,m2,m3,m4],'|j|','S','E_{||}','\rho');
%htitle = get(hleg,'Title');
%set(htitle,'String','maxima')
xlabel('x')

yu1 = 150;
yu2 = 550;
b1 = abs(brotz(yu1));
b2 = abs(brotz(yu2));
rho1 = rho_interped(yu1);
rho2 = rho_interped(yu2);
vout = sqrt(b1*b2*(b1+b2)/(rho1*b2+rho2*b1));

% plot([ppp(yu1) ppp(yu1)],[-2 4])
% plot([ppp(yu2) ppp(yu2)],[-2 4])

E_rate = 0.1*b1*b2/(b1+b2)*vout
