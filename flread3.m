function [time]=flread3(file_number);
global nflx nfly nflz
global fepo

fnh='fltkm';  fnt=num2str(file_number);
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
           time = fread(fid,1,rl);
           hr2 = fread(fid,1,hi);

           hr1 = fread(fid,1,hi);
           nflx = fread(fid,1,it);
           nfly = fread(fid,1,it);
           nflz = fread(fid,1,it);
           hr2 = fread(fid,1,hi);

           fprintf(['time=',num2str(time,'%7.2f'),'\n'])                    
           fprintf(['nflx=',num2str(nflx,'%7.2f'),'\n']);
           fprintf(['nfly=',num2str(nfly,'%7.2f'),'\n']);
           fprintf(['nflz=',num2str(nflz,'%7.2f'),'\n']);
           fprintf(['read: \n']);
           nv=nflx*nfly*nflz*3;        
%line 3 read
           hr1 = fread(fid,1,hi);
           fepo = fread(fid,nv,rl);
           fepo = reshape(fepo,nflx,nfly,nflz,3);
           hr2 = fread(fid,1,hi); 
   fprintf(['\n']);
   fclose(fid);  