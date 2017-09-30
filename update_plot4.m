function update_plot4(y,z,rho,ix,ny,nz)
r = rho(:,:,floor(ix));
%r = reshape(r(:), [ny nz]);
%pcolor(y,z,log10(r'))
pcolor(y,z,r)
shading interp
title(strcat('slice',num2str(floor(ix))));
h = colorbar;
%set(h, 'limits', [min(min(r)) max(max(r))])
caxis([min(min(r)) max(max(r))])
end