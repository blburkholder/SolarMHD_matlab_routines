function [time]=mentrop_read3(file_number);
global pox poy poz xpo ypo zp0 zp1 zps potp potm pgam
global ftlength ftvol ftmass ftjpar ftjperp ftemag ftekin
global ftethe pmaxft pminft zmaxft zminft bminft bbaseft
global xfin1 yfin1 zfin1 bzfin1 dist1

fnh='mentrop';  fnt=num2str(file_number);
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
pox = fread(fid,1,it);
poy = fread(fid,1,it);
poz = fread(fid,1,it);
time = fread(fid,1,rl);
hr2 = fread(fid,1,hi);
fprintf(['time=',num2str(time,'%7.2f'),'\n'])
fprintf(['pox=',num2str(pox,'%4d'),' poy=',num2str(poy,'%4d'),' poz=',num2str(poz,'%4d'),'\n'])
  
   hr1 = fread(fid,1,hi);
     xpo = fread(fid,pox,rl);
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     ypo = fread(fid,poy,rl);
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     zp0 = fread(fid,1,rl);
     zp1 = fread(fid,1,rl);
     zps = fread(fid,1,rl);
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     potp = fread(fid,pox*poy,rl); potp = reshape(potp,pox,poy)';
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     potm = fread(fid,pox*poy,rl); potm = reshape(potm,pox,poy)';
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     pgam = fread(fid,pox*poy,rl); pgam = reshape(pgam,pox,poy)';
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     ftlength = fread(fid,pox*poy,rl); ftlength = reshape(ftlength,pox,poy)';
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     ftvol = fread(fid,pox*poy,rl); ftvol = reshape(ftvol,pox,poy)';
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     ftmass = fread(fid,pox*poy,rl); ftmass = reshape(ftmass,pox,poy)';
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     ftjpar = fread(fid,pox*poy,rl); ftjpar = reshape(ftjpar,pox,poy)';
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     ftjperp = fread(fid,pox*poy,rl); ftjperp = reshape(ftjperp,pox,poy)';
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     ftemag = fread(fid,pox*poy,rl); ftemag = reshape(ftemag,pox,poy)';
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     ftekin = fread(fid,pox*poy,rl); ftekin = reshape(ftekin,pox,poy)';
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     ftethe = fread(fid,pox*poy,rl); ftethe = reshape(ftethe,pox,poy)';
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     pmaxft = fread(fid,pox*poy,rl); pmaxft = reshape(pmaxft,pox,poy)';
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     pminft = fread(fid,pox*poy,rl); pminft = reshape(pminft,pox,poy)';
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     zmaxft = fread(fid,pox*poy,rl); zmaxft = reshape(zmaxft,pox,poy)';
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     zminft = fread(fid,pox*poy,rl); zminft = reshape(zminft,pox,poy)';
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     bminft = fread(fid,pox*poy,rl); bminft = reshape(bminft,pox,poy)';
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     bbaseft = fread(fid,pox*poy,rl); bbaseft = reshape(bbaseft,pox,poy)';
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     xfin1 = fread(fid,pox*poy,rl); xfin1 = reshape(xfin1,pox,poy)';
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     yfin1 = fread(fid,pox*poy,rl); yfin1 = reshape(yfin1,pox,poy)';
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     zfin1 = fread(fid,pox*poy,rl); zfin1 = reshape(zfin1,pox,poy)';
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     bzfin1 = fread(fid,pox*poy,rl); bzfin1 = reshape(bzfin1,pox,poy)';
   hr2 = fread(fid,1,hi);

   hr1 = fread(fid,1,hi);
     dist1 = fread(fid,pox*poy,rl); dist1 = reshape(dist1,pox,poy)';
   hr2 = fread(fid,1,hi);

   fprintf(['\n']);
   fclose(fid); 
end