function [x,y,rho,p,vx,vy,vz,bx,by,bz] = twoD_read(fname) 
fid = fopen(fname,'rb');

hr1 = fread(fid, 1, 'int32');           % Read the first record start tag. Returns hr1 = 72
nx = fread(fid, 1, 'int32');             % Read the second record start tag. Returns hr2 = 12
ny = fread(fid, 1, 'int32');             % Read the second record start tag. Returns hr2 = 12
time = fread(fid, 1, 'float32');             % Read the second record start tag. Returns hr2 = 12
hr1 = fread(fid, 1, 'int32');           % Read the first record start tag. Returns hr1 = 72

hr1 = fread(fid, 1, 'int32');         % Read the first record start tag. Returns hr1 = 72
hr1 = fread(fid, 1, 'int32');          % Read the first record start tag. Returns hr1 = 72
hr1 = fread(fid, 1, 'int32');           % Read the first record start tag. Returns hr1 = 72
x = fread(fid, nx, 'float32');           % Read the first record start tag. Returns hr1 = 72
hr1 = fread(fid, nx, 'float32');           % Read the first record start tag. Returns hr1 = 72
hr1 = fread(fid, nx, 'float32');           % Read the first record start tag. Returns hr1 = 72
hr1 = fread(fid, nx, 'float32');           % Read the first record start tag. Returns hr1 = 72
hr1 = fread(fid, nx, 'float32');           % Read the first record start tag. Returns hr1 = 72
hr1 = fread(fid, nx, 'float32');          % Read the first record start tag. Returns hr1 = 72
hr1 = fread(fid, nx, 'float32');           % Read the first record start tag. Returns hr1 = 72

y = fread(fid, ny, 'float32');           % Read the first record start tag. Returns hr1 = 72
hr1 = fread(fid, ny, 'float32');           % Read the first record start tag. Returns hr1 = 72
hr1 = fread(fid, ny, 'float32');           % Read the first record start tag. Returns hr1 = 72
hr1 = fread(fid, ny, 'float32');           % Read the first record start tag. Returns hr1 = 72
hr1 = fread(fid, ny, 'float32');           % Read the first record start tag. Returns hr1 = 72
hr1 = fread(fid, ny, 'float32');          % Read the first record start tag. Returns hr1 = 72
hr1 = fread(fid, ny, 'float32');           % Read the first record start tag. Returns hr1 = 72

hr1 = fread(fid, 1, 'float32');           % Read the first record start tag. Returns hr1 = 72
hr1 = fread(fid, 1, 'float32');           % Read the first record start tag. Returns hr1 = 72
bx = fread(fid, nx*ny, 'float32');   bx = reshape(bx,[nx,ny])';        % Read the first record start tag. Returns hr1 = 72
by = fread(fid, nx*ny, 'float32');   by = reshape(by,[nx,ny])';        % Read the first record start tag. Returns hr1 = 72
bz = fread(fid, nx*ny, 'float32');   bz = reshape(bz,[nx,ny])';        % Read the first record start tag. Returns hr1 = 72

hr1 = fread(fid, 1, 'int32');           % Read the first record start tag. Returns hr1 = 72
hr1 = fread(fid, 1, 'int32');           % Read the first record start tag. Returns hr1 = 72
sx = fread(fid, nx*ny, 'float32');   sx = reshape(sx,[nx,ny])';        % Read the first record start tag. Returns hr1 = 72
sy = fread(fid, nx*ny, 'float32');   sy = reshape(sy,[nx,ny])';        % Read the first record start tag. Returns hr1 = 72
sz = fread(fid, nx*ny, 'float32');   sz = reshape(sz,[nx,ny])';        % Read the first record start tag. Returns hr1 = 72

hr1 = fread(fid, 1, 'int32');           % Read the first record start tag. Returns hr1 = 72
hr1 = fread(fid, 1, 'int32');           % Read the first record start tag. Returns hr1 = 72
rho = fread(fid, nx*ny, 'float32');   rho = reshape(rho,[nx,ny])';        % Read the first record start tag. Returns hr1 = 72
u = fread(fid, nx*ny, 'float32');   u = reshape(u,[nx,ny])';        % Read the first record start tag. Returns hr1 = 72

vx = sx./rho;
vy = sy./rho;
vz = sz./rho;
p = 2*u.^(5/3);

fclose(fid);
end