function [time]=read3(file_number);
global nx ny nz  x y z difx dify difz
global bx by bz b0x b0y b0z sx sy sz rho u prof res lambda

fnh='magtap';  fnt=num2str(file_number);
for i=1:5
fn=[fnh,fnt];
[fid,message] = fopen (fn,'r','n');  
if fid~=-1, break,  end
fnh=[fnh,'0'];
end
fprintf(['input file is ',fn,'\n'])

hi='int64'; rl='float32'; it='int32';
   
%line 1 read
           hr1 = fread(fid,1,hi);
%check 32 or 64 system
           if hr1 ~=24,  hi='int32';  frewind(fid);  hr1=fread(fid,1,hi);   end;  hi
            nx = fread(fid,1,it);
            ny = fread(fid,1,it);
            nz = fread(fid,1,it);
          time = fread(fid,1,rl);
             f = fread(fid,2,it);
           hr2 = fread(fid,1,hi);
           fprintf(['time=',num2str(time,'%7.2f'),'\n'])
           fprintf(['nx=',num2str(nx,'%4d'),' ny=',num2str(ny,'%4d'),' nz=',num2str(nz,'%4d'),'\n'])
           
           
%line 2 read
      %line 2 read           
           hr1 = fread(fid,1,hi);
             x = fread(fid,nx,rl);
          difx = fread(fid,nx,rl);
         ddifx = fread(fid,nx,rl);
        ddifpx = fread(fid,nx,rl);
        ddifmx = fread(fid,nx,rl);
        meanpx = fread(fid,nx,rl);
        meanmx = fread(fid,nx,rl);
             y = fread(fid,ny,rl);
          dify = fread(fid,ny,rl);
         ddify = fread(fid,ny,rl);
        ddifpy = fread(fid,ny,rl);
        ddifmy = fread(fid,ny,rl);
        meanpy = fread(fid,ny,rl);
        meanmy = fread(fid,ny,rl);
             z = fread(fid,nz,rl);
          difz = fread(fid,nz,rl);
         ddifz = fread(fid,nz,rl);
        ddifpz = fread(fid,nz,rl);
        ddifmz = fread(fid,nz,rl);
        meanpz = fread(fid,nz,rl);
        meanmz = fread(fid,nz,rl);

           hr2 = fread(fid,1,hi);
           xmax=max(x);xmin=min(x);
           ymax=max(y);ymin=min(y);
           zmax=max(z);zmin=min(z);
       fprintf(['xmin=',num2str(xmin,'%7.2f'),' xmax=',num2str(xmax,'%7.2f'),'\n']);
       fprintf(['ymin=',num2str(ymin,'%7.2f'),' ymax=',num2str(ymax,'%7.2f'),'\n']);
       fprintf(['zmin=',num2str(zmin,'%7.2f'),' zmax=',num2str(zmax,'%7.2f'),'\n']);
       fprintf(['read: \n']);
       nv=nx*ny*nz;        
%line 3 read
           hr1 = fread(fid,1,hi);
            bx = fread(fid,nv,rl); bx = reshape(bx,nx,ny,nz); fprintf(['bx, ']); bx = permute(bx, [2 1 3]);
            by = fread(fid,nv,rl); by = reshape(by,nx,ny,nz); fprintf(['by, ']); by = permute(by, [2 1 3]);
            bz = fread(fid,nv,rl); bz = reshape(bz,nx,ny,nz); fprintf(['bz, ']); bz = permute(bz, [2 1 3]);
           hr2 = fread(fid,1,hi);
%line 4 read
           hr1 = fread(fid,1,hi);
            b0x = fread(fid,nv,rl);  b0x = reshape(b0x,nx,ny,nz); fprintf(['b0x, ']); b0x = permute(b0x, [2 1 3]);
            b0y = fread(fid,nv,rl);  b0y = reshape(b0y,nx,ny,nz); fprintf(['b0y, ']); b0y = permute(b0y, [2 1 3]);
            b0z = fread(fid,nv,rl);  b0z = reshape(b0z,nx,ny,nz); fprintf(['b0z, ']); b0z = permute(b0z, [2 1 3]);
           hr2 = fread(fid,1,hi);
%line 5 read
           hr1 = fread(fid,1,hi);
            sx = fread(fid,nv,rl);  sx = reshape(sx,nx,ny,nz); fprintf(['sx, ']); sx = permute(sx, [2 1 3]);
            sy = fread(fid,nv,rl);  sy = reshape(sy,nx,ny,nz); fprintf(['sy, ']); sy = permute(sy, [2 1 3]);
            sz = fread(fid,nv,rl);  sz = reshape(sz,nx,ny,nz); fprintf(['sz, ']); sz = permute(sz, [2 1 3]);
           hr2 = fread(fid,1,hi);

%line 6 read
           hr1 = fread(fid,1,hi);
           rho  = fread(fid,nv,rl);  rho = reshape( rho,nx,ny,nz); fprintf(['rho, ']); rho = permute(rho, [2 1 3]);
             u  = fread(fid,nv,rl);    u = reshape(   u,nx,ny,nz); fprintf(['u, ']); u = permute(u, [2 1 3]);
           res  = fread(fid,nv,rl);  res = reshape( res,nx,ny,nz); fprintf(['res, ']); res = permute(res, [2 1 3]);
          prof  = fread(fid,nv,rl); prof = reshape(prof,nx,ny,nz); fprintf(['prof, ']); prof = permute(prof, [2 1 3]);
            hr2 = fread(fid,1,hi);

%line 7 read
           hr1 = fread(fid,1,hi);
           lambda  = fread(fid,nv,rl);  lambda = reshape( lambda,nx,ny,nz); fprintf(['lambda, ']); lambda = permute(lambda, [2 1 3]);
           hr2 = fread(fid,1,hi);
   fprintf(['\n']);

   
   
   fclose(fid);           

           
           
           
 
  
  
  
   