function prepare2(fn)
global nx ny nz x y z difx dify difz
global bx by bz b0x b0y b0z sx sy sz rho u prof res lambda
%fn=input('please gieve the file number. ');
time=read3(fn);
gamma=5/3;
fprintf(['smooth begin.\n'])
%smooth the original data
bx=smooth3(bx); by=smooth3(by); bz=smooth3(bz);
b0x=smooth3(b0x); b0y=smooth3(b0y); b0z=smooth3(b0z);
sx=smooth3(sx); sy=smooth3(sy); sz=smooth3(sz);
rho=smooth3(rho); u=smooth3(u); res=smooth3(res); 
%got new physics quantity
fprintf(['New Physics Quantity:'])
S=(u./rho); S=S.^gamma;                    fprintf(['S,']);
p=2.*u.^gamma;                             fprintf(['p,']);
%j=sqrt(jx.*jx+jy.*jy+jz.*jz);              fprintf(['j,']);
%jc0=sqrt(p.*rho);                          fprintf(['j0,']);
B=sqrt(bx.*bx+by.*by+bz.*bz);              fprintf(['b,']);
%EB=res.*(jx.*bx+jy.*by+jz.*bz);            fprintf(['EB,']);
%Epara=EB./(B+eps);                         fprintf(['Epara,']);
%Edif=res.*j;                               fprintf(['Edif.\n']);



pname=['mhd',int2str(fn),'.mat'];
save(pname,'time','nx','ny','nz','x','y','z', 'difx','dify', 'difz',...
            'bx', 'by', 'bz', ...
            'b0x', 'b0y', 'b0z', ...
            'sx', 'sy', 'sz', ...
            'rho', 'u', 'prof', 'res','lambda',...
             'S','p')
