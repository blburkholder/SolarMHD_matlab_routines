  load('mentrop6.mat')
  ftepar(ftepar > 0.4) = 0.4;
  ftepar(ftepar < -0.05) = -0.05;
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



