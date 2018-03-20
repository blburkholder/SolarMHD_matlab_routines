[jx,jy,jz] = get_j(nx,ny,nz,bx,by,bz,b0x,b0y,b0z,difx,dify,difz);
figure
platypus = patch(isosurface(x,y,z,rho,0.1));
platypus.FaceVertexCData = platypus.Vertices(:,end);
platypus.FaceColor = 'flat';
platypus.EdgeColor = 'none';
daspect([1 1 1])
axis([0 100 0 100 0 70])

%calculate JxB forces
jxb_x = jy.*bz - jz.*by;
jxb_y = jz.*bx - jx.*bz;
jxb_z = jx.*by - jy.*bx;
%calculate pressure gradient forces -grad(P)
gp_z = (p(2:end-1,2:end-1,1:end-2) - p(2:end-1,2:end-1,3:end))/(2*(x(2)-x(1)));
gp_y = (p(2:end-1,1:end-2,2:end-1) - p(2:end-1,3:end,2:end-1))/(2*(x(2)-x(1)));
gp_x = (p(1:end-2,2:end-1,2:end-1) - p(3:end,2:end-1,2:end-1))/(2*(x(2)-x(1)));
 
%project forces onto field lines
%jxb_p = (jxb_x.*bx + jxb_y.*by + jxb_z.*bz)./sqrt(bx.^2 + by.^2 + bz.^2);
gp_p = (gp_x.*bx(2:end-1,2:end-1,2:end-1) + gp_y.*by(2:end-1,2:end-1,2:end-1) +...
    gp_z.*bz(2:end-1,2:end-1,2:end-1))./sqrt(bx(2:end-1,2:end-1,2:end-1).^2 +...
    by(2:end-1,2:end-1,2:end-1).^2 + bz(2:end-1,2:end-1,2:end-1).^2);

beta = p./(bx.^2+by.^2+bz.^2);

s1 = 1;
s2 = 60;

for i = 40:10:80
    figure
    pcolor(x,y,jxb_z(:,:,i))
    %pcolor(x(2:end-1),y(2:end-1),sz(2:end-1,2:end-1,i)./rho(2:end-1,2:end-1,i))
    shading interp
    colorbar
    title(strcat('jxb_z',num2str(z(i))))

%     figure
%     pcolor(x(2:end-1),y(2:end-1),gp_z(:,:,i))
%     shading interp
%     colorbar
%     title(strcat('\nabla p_{z}',num2str(z(i))))
%   saveas(gca,strcat('assorted_figs/curtain',num2str(floor(10*z(i))),'.png'))

%     figure
%     pcolor(beta(:,:,i))
%     shading interp
%     colorbar
%     title(strcat('z=',num2str(z(i))))
%     %saveas(gca,strcat('assorted_figs/beeta',num2str(floor(10*z(i))),'.png'))
%     pause
end