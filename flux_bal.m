%prepare3(0)
%prepare3(10)

load('mhd0.mat')
net = zeros(nz-2,1);
tot = zeros(nz-2,1);
for i = 2:nz-1
    [net(i-1),tot(i-1)] = flux(bz(:,:,i),x(2)-x(1),y(2)-y(1));
end
    
net1 = zeros(nz-2,1);
tot1 = zeros(nz-2,1);
load('mhd10.mat')
for i = 2:nz-1
    [net1(i-1),tot1(i-1)] = flux(bz(:,:,i),x(2)-x(1),y(2)-y(1));
end

f = figure;
ax = axes('Parent',f);

subplot(3,2,2)
semilogy(z(2:end-1),tot1)
title('total flux t_f')
subplot(3,2,1)
plot(z(2:end-1),net1)
title('net flux t_f')
subplot(3,2,5)
plot(z(2:end-1),net - net1)
title('net flux diff t_i - t_f')
subplot(3,2,6)
plot(z(2:end-1),tot - tot1)
title('total flux diff t_i - t_f')
subplot(3,2,3)
plot(z(2:end-1),net)
title('net flux t_i')
subplot(3,2,4)
semilogy(z(2:end-1),tot)
title('total flux t_i')





