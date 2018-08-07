function [xf,yf] = velocity_integrator(x,y,xi,yi,tf,vx,vy)
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
t
xf = reshape(start_point(:,1),[size(xi)]);
yf = reshape(start_point(:,2),[size(yi)]);



