load('mhd0.mat')

net1 = zeros(nz-2,1);
tot1 = zeros(nz-2,1);
for i = 2:nz-1
    [net1(i-1),tot1(i-1)] = flux(bz(:,:,i),x(2)-x(1),y(2)-y(1));
end
    
net2 = zeros(nz-2,1);
tot2 = zeros(nz-2,1);

load('mhd3.mat')

for i = 2:nz-1
    [net2(i-1),tot2(i-1)] = flux(bz(:,:,i),x(2)-x(1),y(2)-y(1));
end

net3 = zeros(nz-2,1);
tot3 = zeros(nz-2,1);

load('mhd6.mat')

for i = 2:nz-1
    [net3(i-1),tot3(i-1)] = flux(bz(:,:,i),x(2)-x(1),y(2)-y(1));
end

net4 = zeros(nz-2,1);
tot4 = zeros(nz-2,1);

load('mhd9.mat')

for i = 2:nz-1
    [net4(i-1),tot4(i-1)] = flux(bz(:,:,i),x(2)-x(1),y(2)-y(1));
end

net5 = zeros(nz-2,1);
tot5 = zeros(nz-2,1);

load('mhd12.mat')

for i = 2:nz-1
    [net5(i-1),tot5(i-1)] = flux(bz(:,:,i),x(2)-x(1),y(2)-y(1));
end

f = figure;
ax = axes('Parent',f);

subplot(3,2,2)
semilogy(z(2:end-1),tot1)
title('total flux t_0')
subplot(3,2,1)
plot(z(2:end-1),net1)
title('net flux t_0')

net = net1;
tot = tot1;

subplot(3,2,5)
plot(z(2:end-1),net - net1)
title('net flux diff t_i - t_0')
subplot(3,2,6)
plot(z(2:end-1),tot - tot1)
title('total flux diff t_i - t_0')
subplot(3,2,3)
plot(z(2:end-1),net)
title('net flux t_i')
subplot(3,2,4)
semilogy(z(2:end-1),tot)
title('total flux t_i')

i = 0;
b = uicontrol('Parent',f,'Style','slider','Position',[81,0,419,23],...
              'value',i, 'min',0, 'max',5);
b.Callback  = @(objHandle,~)  update_plot2(z,net1,tot1,net2,tot2,net3,tot3,net4,tot4,net5,tot5,get(objHandle,'Value'))








