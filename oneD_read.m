function [x,rho,p,temp,vx,vy,vz,bx,by,bz] = oneD_read(fname) 

fid = fopen(fname,'rb');

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
temp = zeros(2001,iend/iout);


for i = 1:iend/iout
    a = fread(fid,nx*9,'float32');
    rho(:,i) = a(2002:4002);
    p(:,i) = a(4003:6003);
    temp(:,i) = p(:,i)./rho(:,i);
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




