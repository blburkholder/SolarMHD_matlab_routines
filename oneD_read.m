%function [a] = oneD_read() 

fid = fopen('m1bin','rb');

hr1=fread(fid, 1, 'int32');            % Read the first record start tag. Returns hr1 = 72
nx = fread(fid, 1, 'int32');             % Read the second record start tag. Returns hr2 = 12
hr2=fread(fid, 1, 'int32');             % Read the second record start tag. Returns hr2 = 12
hr3=fread(fid, 1, 'int32');             % Read the second record start tag. Returns hr2 = 12
iout=fread(fid, 1, 'int32');             % Read the second record start tag. Returns hr2 = 12
iend=fread(fid, 1, 'int32');             % Read the second record start tag. Returns hr2 = 12
intart=fread(fid, 1, 'int32');             % Read the second record start tag. Returns hr2 = 12
intu=fread(fid, 1, 'int32');             % Read the second record start tag. Returns hr2 = 12
ivisrho=fread(fid, 1, 'int32');             % Read the second record start tag. Returns hr2 = 12
ivisv=fread(fid, 1, 'int32');             % Read the second record start tag. Returns hr2 = 12
ivisu=fread(fid, 1, 'int32');             % Read the second record start tag. Returns hr2 = 12
igrid=fread(fid, 1, 'int32');             % Read the second record start tag. Returns hr2 = 12
idiag=fread(fid, 1, 'int32');             % Read the second record start tag. Returns hr2 = 12
init=fread(fid, 1, 'int32');             % Read the second record start tag. Returns hr2 = 12

hr4=fread(fid, 1, 'int32');             % Read the second record start tag. Returns hr2 = 12
hr5=fread(fid, 1, 'int32');             % Read the second record start tag. Returns hr2 = 12
dt=fread(fid, 1, 'float32');             % Read the second record start tag. Returns hr2 = 12
gamma=fread(fid, 1, 'float32');             % Read the second record start tag. Returns hr2 = 12
resist=fread(fid, 1, 'float32');             % Read the second record start tag. Returns hr2 = 12
visrho=fread(fid, 1, 'float32');             % Read the second record start tag. Returns hr2 = 12
visv=fread(fid, 1, 'float32');             % Read the second record start tag. Returns hr2 = 12
visu=fread(fid, 1, 'float32');             % Read the second record start tag. Returns 
xmax =fread(fid, 1, 'float32');             % Read the second record start tag. Returns hr2 = 12
hr = fread(fid,1,'int32');
hr = fread(fid,1,'int32');
timestep = fread(fid,1,'int32');
time = fread(fid,1,'int32');

hr = fread(fid,1,'int32');
hr = fread(fid,1,'float32');
bx = zeros(2001,iend/iout);
by = zeros(2001,iend/iout);
bz = zeros(2001,iend/iout);
vx = zeros(2001,iend/iout);
vy = zeros(2001,iend/iout);
vz = zeros(2001,iend/iout);
rho = zeros(2001,iend/iout);
p = zeros(2001,iend/iout);


for i = 1:iend/iout
    a = fread(fid,nx*9,'float32');
    rho(:,i) = a(2002:4002);
    p(:,i) = a(4003:6003);
    vx(:,i) = a(6004:8004);
    vy(:,i) = a(8005:10005);
    vz(:,i) = a(10006:12006);
    bx(:,i) = a(12007:14007);
    by(:,i) = a(14008:16008);
    bz(:,i) = a(16009:18009);

    timestep = fread(fid,1,'int32');
    time = fread(fid,1,'int32');

    hr = fread(fid,1,'int32');
    hr = fread(fid,1,'int32');
    hr = fread(fid,1,'int32');
    hr = fread(fid,1,'int32');
end
x = a(1:2001);
fclose(fid);

    figure
    subplot(5,1,1); hold on
    plot(x,bx(:,1),'r')
    plot(x,by(:,1),'g')
    plot(x,bz(:,1),'b')
    title('magnetic field')
    legend('B_x','B_y','B_z')
    axis([-2 2 -1.1 1.1])
    set(gca,'xticklabel',[])

    subplot(5,1,2);
    plot(x,rho(:,1))
    title('density')
    axis([-2 2 0.8 1.5])
    set(gca,'xticklabel',[])
    
    subplot(5,1,5); hold on
    plot(x,abs(bx(:,1))./sqrt(rho(:,1)),'r')
    plot(x,abs(by(:,1))./sqrt(rho(:,1)),'g')
    plot(x,abs(bz(:,1))./sqrt(rho(:,1)),'b')
    title('Alfven velocity')
    xlabel('x')
    legend('v_{Ax}','v_{Ay}','v_{Az}')
    axis([-2 2 -0.1 1.1])

    subplot(5,1,3); hold on
    plot(x,vx(:,1),'r')
    plot(x,vy(:,1),'g')
    plot(x,vz(:,1),'b')
    title('velocity')
    set(gca,'xticklabel',[])
    legend('v_x','v_y','v_z')
    axis([-2 2 -1 1])

    subplot(5,1,4)
    plot(x,p(:,1)./rho(:,1).^(5/3))
    title('entropy')
    axis([-2 2 0.9 1.2])
    set(gca,'xticklabel',[])

[ind1,ind2] = simulation_walen(x,rho(:,end),bx(:,end),by(:,end),bz(:,end),vx(:,end),vy(:,end),vz(:,end),...
        abs(bx(:,end))./rho(:,end),abs(by(:,end))./rho(:,end),abs(bz(:,end))./rho(:,end),0.9,1.1);
subplot(5,1,4)
plot(x,p(:,end)./rho(:,end).^(5/3))
title('entropy')
axis([-2 2 0.9 1.2])


