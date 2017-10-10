%[jx,jy,jz,P] = get_j(nx,ny,nz,bx,by,bz,b0x,b0y,b0z,res,difx,dify,difz,0);
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
gp_x = p(2:end-1,2:end-1,1:end-2) - p(2:end-1,2:end-1,3:end)/(2*(x(2)-x(1)));
gp_y = p(2:end-1,1:end-2,2:end-1) - p(2:end-1,3:end,2:end-1)/(2*(x(2)-x(1)));
gp_z = p(1:end-2,2:end-1,2:end-1) - p(3:end,2:end-1,2:end-1)/(2*(x(2)-x(1)));
 
%project forces onto field lines
% jxb_p = (jxb_x.*bx + jxb_y.*by + jxb_z.*bz)./sqrt(bx.^2 + by.^2 + bz.^2);
gp_p = (gp_x.*bx(2:end-1,2:end-1,2:end-1) + gp_y.*by(2:end-1,2:end-1,2:end-1) +...
    gp_z.*bz(2:end-1,2:end-1,2:end-1))./sqrt(bx(2:end-1,2:end-1,2:end-1).^2 +...
    by(2:end-1,2:end-1,2:end-1).^2 + bz(2:end-1,2:end-1,2:end-1).^2);

beta = p./(bx.^2+by.^2+bz.^2);

s1 = 1;
s2 = 60;

for i = s1:2:s2
    figure
    pcolor(jxb_x(:,:,i).^2+jxb_y(:,:,i).^2+jxb_z(:,:,i).^2)
    shading interp
    colorbar
    title(strcat('jxb z=',num2str(z(i))))

    figure
    pcolor(gp_p(:,:,i))
    shading interp
    colorbar
    title(strcat('\nabla p z=',num2str(z(i))))
%    saveas(gca,strcat('assorted_figs/curtain',num2str(floor(10*z(i))),'.png'))

%     figure
%     pcolor(beta(:,:,i))
%     shading interp
%     colorbar
%     title(strcat('z=',num2str(z(i))))
%     %saveas(gca,strcat('assorted_figs/beeta',num2str(floor(10*z(i))),'.png'))
%     pause
end