function update_plot3(x,y,null,iz)
pcolor(x,y,null(:,:,floor(iz))')
shading interp
title(strcat('slice',num2str(floor(iz))));
end