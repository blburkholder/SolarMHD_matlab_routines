function flprepare2(fn)
global nflx nfly nflz
global fepo
%fn=input('please gieve the file number. ');
time=flread3(fn);

pname=['flt',int2str(fn),'.mat'];
save(pname,'time','nflx','nfly','nflz','fepo')
