  load('mentrop10.mat')
  ftepar(ftepar > 0.3) = 0.3;
  ftepar(ftepar < -0.02) = -0.02;
  %ftepar(ftepar == 0) = 0.001;
  figure
  %scatter(xpo(:),ypo(:),[],ftepar(:),'.')
  pcolor(xpo,ypo,ftepar)
  colormap(jet)
  shading interp
  colorbar
  daspect([1 1 1])
  xlabel('X')
  ylabel('Y')
  title('$\int E_{||}$','interpreter','latex')



