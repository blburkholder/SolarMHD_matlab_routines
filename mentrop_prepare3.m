function mentrop_prepare3(fn)
global pox poy poz xpo ypo zpo zp0 zp1 potp potm pgam
global ftlength ftvol ftmass ftjpar ftjperp ftemag ftekin ftepar
global ftethe pmaxft pminft zmaxft zminft bminft bbaseft
global xfin1 yfin1 zfin1 xfin2 yfin2 zfin2 bzfin1 bzfin2 dist1 dist2

time=mentrop_read4(fn);
pname=['flt_mentrop',int2str(fn),'.mat'];
save(pname,'time','pox', 'poy', 'poz', 'xpo', 'ypo', 'zpo',...
    'zp0', 'zp1', 'potp', 'potm', 'pgam', 'ftlength', 'ftvol',...
    'ftmass', 'ftjpar', 'ftjperp', 'ftemag', 'ftekin','ftepar',...
    'ftethe', 'pmaxft', 'pminft', 'zmaxft', 'zminft', 'bminft',...
    'bbaseft', 'xfin1', 'yfin1', 'zfin1','xfin2', 'yfin2',...
    'zfin2', 'bzfin1','bzfin2','dist1','dist2')