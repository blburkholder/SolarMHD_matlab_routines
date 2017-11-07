function update_plot(x,y,ll,iz)
pcolor(x,y,ll(:,:,floor(iz)))
shading interp
title(strcat('slice',num2str(floor(iz))));
colorbar
end