function [xf,yf] = velocity_perturbation(x,y,xi,yi,tf)
nx = length(x);
ny = length(y);

vx0 = 1.0;
vy0 = 0;
vz0 = 0;
vw0 = 0;

rho0 = 0.072;

argx1 = zeros(nx,ny);
argx2 = zeros(nx,ny); 
argy1 = zeros(nx,ny); 
argy2 = zeros(nx,ny);
argz1 = zeros(nx,ny); 
argz2 = zeros(nx,ny); 
argw1 = zeros(nx,ny); 
argw2 = zeros(nx,ny); 

vprof1 = zeros(nx,ny);  
vprof2 = zeros(nx,ny);  
vprof3 = zeros(nx,ny); 
vprof4 = zeros(nx,ny);  

vxprof = zeros(nx,ny);  
vyprof = zeros(nx,ny);      
vzprof = zeros(nx,ny);  

  lx1=-20;
  lx2=20;
  cx1=0;
  cx2=-92.8;
  ly1=-12;
  ly2=12;
  cy1=-45;
  cy2=-85;
  lz1=-12;
  lz2=12;
  cz1=40;
  cz2=-80;
  lw1=-12;
  lw2=12;
  cw1=0;
  cw2=-40; 

for iy = 1:ny
    for ix = 1:nx
        argx1(ix,iy)  = (y(iy)-x(ix)+cx1)/lx1;
        argx2(ix,iy)  = (y(iy)+x(ix)+cx2)/lx2;
        argy1(ix,iy)  = (y(iy)-x(ix)+cy1)/ly1;
        argy2(ix,iy)  = (y(iy)+x(ix)+cy2)/ly2;
        argz1(ix,iy)  = (y(iy)-x(ix)+cz1)/lz1;
        argz2(ix,iy)  = (y(iy)+x(ix)+cz2)/lz2;
        argw1(ix,iy)  = (y(iy)-x(ix)+cw1)/lw1;
        argw2(ix,iy)  = (y(iy)+x(ix)+cw2)/lw2;		
    end
end
for iy = 1:ny
    for ix = 1:nx
        vprof1(ix,iy)  = 1./cosh(argx1(ix,iy))/cosh(argx2(ix,iy));
        vprof2(ix,iy)  = 1./cosh(argy1(ix,iy))/cosh(argy2(ix,iy));
        vprof3(ix,iy)  = 1./cosh(argz1(ix,iy))/cosh(argz2(ix,iy));
        vprof4(ix,iy)  = 1./cosh(argw1(ix,iy))/cosh(argw2(ix,iy));
    end
end
for iy = 1:ny
    for ix = 1:nx
        vxprof(ix,iy) = -vx0*(1/lx1*tanh(argx1(ix,iy)) +...
            1/lx2*tanh(argx2(ix,iy)) )*vprof1(ix,iy) +...
            vy0*(1/ly1*tanh(argy1(ix,iy)) + 1/ly2*tanh(argy2(ix,iy)))*vprof2(ix,iy)...
            -vz0*(1/lz1*tanh(argz1(ix,iy)) + 1/lz2*tanh(argz2(ix,iy)))*vprof3(ix,iy)...
            +vw0*(1/lw1*tanh(argw1(ix,iy)) + 1/lw2*tanh(argw2(ix,iy)))*vprof4(ix,iy); 
        vyprof(ix,iy) = -vx0*(1./lx1*tanh(argx1(ix,iy)) -... 
            1/lx2*tanh(argx2(ix,iy)))*vprof1(ix,iy) +...
            vy0*(1/ly1*tanh(argy1(ix,iy)) - 1/ly2*tanh(argy2(ix,iy)))*vprof2(ix,iy)...
            - vz0*(1/lz1*tanh(argz1(ix,iy)) - 1/lz2*tanh(argz2(ix,iy)))*vprof3(ix,iy)...
            + vw0*(1/lw1*tanh(argw1(ix,iy)) - 1/lw2*tanh(argw2(ix,iy)))*vprof4(ix,iy); 
        vzprof(ix,iy) = 0.0;
    end
end

rhoprof = 5.*(1-tanh(-2))+rho0;

rho = zeros(nx,ny);
for iy = 1:ny
    for ix = 1:nx
        rho(ix,iy) = rhoprof;
    end
end

sx = zeros(nx,ny);
sy = zeros(nx,ny);
for iy = 1:ny
    for ix = 1:nx
        sx(ix,iy) = (rhoprof-rho0)*vxprof(ix,iy);
        sy(ix,iy) = (rhoprof-rho0)*vyprof(ix,iy);
    end
end

vx = zeros(nx,ny);
vy = zeros(nx,ny);
for iy = 1:ny
    for ix =1:nx
        rhoinv = 1./rho(ix,iy);
        vx(ix,iy) = sx(ix,iy)*rhoinv;
        vy(ix,iy) = sy(ix,iy)*rhoinv;
    end
end

%this shows that the velocity is divergence free
% div = zeros(nx-2,ny-2);
% for i = 2:nx-1
%     for j = 2:ny-1
%         div(i-1,j-1) = vxprof(i+1,j) - vxprof(i-1,j) + vyprof(i,j+1) - vyprof(i,j-1);
%     end
% end
% figure
% pcolor(div)
% shading interp
% colorbar

% figure
% colormap(jet)
% pcolor(x,y,(vy'.^2 + vx'.^2));
% shading interp
% hold on
% plot([x(1) x(end-1)],[y(127) y(127)],'w')
% plot([x(127) x(127)],[y(2) y(end-1)],'w')
% %h = streamslice(x,y,vx',vy');
% [qx,qy] = meshgrid(x,y);
% h = quiver(qx(1:10:end,1:10:end),qy(1:10:end,1:10:end),vx(1:10:end,1:10:end)',vy(1:10:end,1:10:end)');
% set(h,'Color','k')
% daspect([1 1 1])

%hold on
start_point = [xi(:),yi(:)];

t = 0;
%makes the maximum step length a quarter of the
%fastest travel time across a grid cell
if max(max(vx)) > max(max(vy))
    dt = 0.25*(x(2) - x(1))/max(max(vx));
else
    dt = 0.25*(y(2) - y(1))/max(max(vy));
end
%dt = 0.1;
nstep = floor(tf/dt);
%scatter(start_point(:,1),start_point(:,2),10,'k.')
for k = 1:2
    for i = 1:nstep

        vxx = interp2(x,y,vx',start_point(:,1),start_point(:,2));
        vyy = interp2(x,y,vy',start_point(:,1),start_point(:,2));

        start_point(:,1) = start_point(:,1) + vxx*dt;
        start_point(:,2) = start_point(:,2) + vyy*dt;
        
%        scatter(start_point(:,1),start_point(:,2),10,'k.')
        t = t + dt;
    end
    nstep = 1;
    dt = tf - t;
end
%scatter(start_point(:,1),start_point(:,2),10,'w.')
%title(['tf = ', num2str(t)])
%daspect([1,1,1])
xf = reshape(start_point(:,1),[size(xi)])';
yf = reshape(start_point(:,2),[size(yi)])';

% figure
% plot(x,vx(:,2))
% hold on
% plot(x,vx(:,end-1))
% plot(x,vx(2,:))
% plot(x,vx(end-1,:))
% plot([x(127),x(127)],[-0.005 0.005])
% title('vx line symmetry')
% 
% figure
% plot(y,vy(:,2))
% hold on
% plot(y,vy(:,end-1))
% plot(y,vy(2,:))
% plot(y,vy(end-1,:))
% plot([y(127),y(127)],[-0.005 0.005])
% title('vy line symmetry')
% plot(x,vx(end-1,:))
