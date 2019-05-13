function [brotated_x,brotated_y,brotated_z,jinterped_x,jinterped_y,jinterped_z,epar_interped,mmm,ppp,nnn,dirdir] =...
    tangential_discont(x,y,z,nx,ny,nz,bx,by,bz,b0x,b0y,b0z,difx,dify,difz,sx,sy,sz,res,rho,p,time,txy,tyz,nnz,bl)
    tx1 = txy(1);    ty1 = txy(2);    tx2 = txy(3);    ty2 = txy(4);    tx11 = txy(5);    ty11 = txy(6);    tx22 = txy(7);    ty22 = txy(8);

    %%%%%%%calculate current density and fac density for simulation volume
     [jx,jy,jz] = get_j(nx,ny,nz,bx,by,bz,b0x,b0y,b0z,difx,dify,difz);
     [e_par,j_par,j_perp,~,~,~] = get_j2(nx,ny,nz,res,jx,jy,jz,bx,by,bz,...
        sx./rho,sy./rho,sz./rho);
     j_par = j_par.^2+j_perp.^2;

    %find maximum current density along vertical lines
    txx1 = tx1*ones(length(ty2:0.01:ty1),1);
    tyy1 = ty2:0.01:ty1;
    txx2 = tx11*ones(length(ty22:0.01:ty11),1);
    tyy2 = ty22:0.01:ty11;
    tjj1 = interp2(x,y,e_par(:,:,nnz)',txx1,tyy1');
    tjj2 = interp2(x,y,e_par(:,:,nnz)',txx2,tyy2');
    %tjj1 = interp2(x,y,j_par(:,:,nnz),txx1,tyy1');
    %tjj2 = interp2(x,y,j_par(:,:,nnz),txx2,tyy2');
    tmax_max1 = abs(tjj1) == max(abs(tjj1));
    tx_max1 = txx1(tmax_max1);
    ty_max1 = tyy1(tmax_max1);
    tj_max1 = tjj1(tmax_max1);
    tmax_max2 = abs(tjj2) == max(abs(tjj2));
    tx_max2 = txx2(tmax_max2);
    ty_max2 = tyy2(tmax_max2);
    tj_max2 = tjj2(tmax_max2);

    %find line connecting maximum current along 2 vertical lines to define a
    %normal to the boundary in the 2D x-y plane 
    mmm = (tx_max1-3):0.01:(tx_max2+3);
    %mmm = tx_max1:0.01:tx_max2;
    mslope1 = (ty_max1 - ty_max2)/(tx_max1 - tx_max2); 
    lll = mslope1*((mmm) - tx_max1) + ty_max1;
    %ppp = (-1/mslope1)*((mmm) - tx_max1) + ty1;
    ppp = (-1/mslope1)*((mmm)- (tx_max2+tx_max1)/2) + (ty_max1+ty_max2)/2; 

    diffx = mmm(end)-mmm(1);
    diffy = ppp(end)-ppp(1);
    magg = sqrt(diffx^2+diffy^2);
    dirx = diffx/magg;
    diry = diffy/magg;

    %make plot showing direction in 2D plane at given height
    [xxx,yyy] = meshgrid(x,y);

    figure
    pcolor(x,y,e_par(:,:,nnz)')
    shading interp
    daspect([1 1 1])
    shading interp
    daspect([1 1 1])
    hold on 
    %h = streamslice(x,y,bx(:,:,nnz),by(:,:,nnz),3);
    %set(h,'Color','k')
    plot(txx1,tyy1,'k')
    plot(tx_max1,ty_max1,'rp')
    plot(txx2,tyy2,'k')
    plot(tx_max2,ty_max2,'rp')
    %line connecting the maxima
    plot(mmm,lll,'g')
    plot(mmm,ppp,'g')
    title(['time=',num2str(time),' z=',num2str(z(nnz))])
    colorbar
    set(gca,'xdir','reverse')
    view([90 90])
    xlabel('y')
    ylabel('x')

    %find normal to boundary in 3d geometrically?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_y0 = ty2*mslope1 + tx_max1;
    x_yf = -(y(end-1)-ty2)*mslope1+tx_max1;
    if x_y0 > x_yf
        mm = x_yf:0.01:x_y0;
    else
        mm = x_y0:0.01:x_yf;
    end
    pp = (-1/mslope1)*((mm) - tx_max1) + ty2;

    xx = zeros(length(mm),153);
    yy = zeros(length(mm),153);
    zz = zeros(length(mm),153);
    for i = 1:153
        xx(:,i) = mm;
        yy(:,i) = pp;
        zz(:,i) = z(i)*ones(size(mm));
    end
    %plot(mm,pp)
    ang_slice = interp3(x,y,z,permute(e_par,[2 1 3]),xx,yy,zz);
    %ang_slice = interp3(x,y,z,j_par,xx,yy,zz);
    bzang_slice = interp3(x,y,z,bz,xx,yy,zz);

    %define 2 horizontal lines
    y1 = tyz(1);
    z1 = z(nnz)+1.5;
    y2 = tyz(2);
    z2 = z(nnz)+1.5;

    y11 = tyz(3);
    z11 = z(nnz)-1.5;
    y22 = tyz(4);
    z22 = z(nnz)-1.5;

    zz1 = z1*ones(length(y1:0.01:y2),1);
    yy1 = y1:0.01:y2;
    zz2 = z11*ones(length(y11:0.01:y22),1);
    yy2 = y11:0.01:y22;
    jj1 = interp2(pp,z,abs(ang_slice'),yy1',zz1);
    jj2 = interp2(pp,z,abs(ang_slice'),yy2',zz2);
    max_max1 = jj1 == max(jj1);
    z_max1 = zz1(max_max1);
    y_max1 = yy1(max_max1);
    j_max1 = jj1(max_max1);
    max_max2 = jj2 == max(jj2);
    z_max2 = zz2(max_max2);
    y_max2 = yy2(max_max2);
    j_max2 = jj2(max_max2);

    nnn = z_max2:(abs(z_max2 - z_max1)/(length(mmm)-1)):z_max1;
    mslope = (y_max1 - y_max2)/(z_max1 - z_max2); 
    iii = mslope*((z_max2:0.1:z_max1) - z_max1) + y_max1;
    qqq = (-1/mslope)*((nnn) - z_max1) + y1;

    figure
    %pcolor(pp,z(20:end),ang_slice(:,20:end)')
    pcolor(pp,z(nnz-20:end),ang_slice(:,nnz-20:end)')
    shading interp
    daspect([1 1 1])
    hold on
    reso = 10;
    qz = z(1:65);
    qu = 0*bzang_slice;
    qv = bzang_slice;
    %h = quiver(pp(1:2*reso:end),qz(1:reso:end),qu(1:2*reso:end,1:reso:end)',qv(1:2*reso:end,1:reso:end)');
    %set(h,'Color','k')
    plot(yy1,zz1,'k')
    plot(y_max1,z_max1,'rp')
    plot(yy2,zz2,'k')
    plot(y_max2,z_max2,'rp')
    %line connecting the maxima
    %plot(iii,10:0.1:15,'g','LineWidth',0.1)
    plot(iii,z_max2:0.1:z_max1,'g')
    colorbar
    xlabel('x-y cut')
    ylabel('z')
    title('E_{||}')

    risetot = sqrt((mmm(1)-mmm(end))^2+(ppp(1)-ppp(end))^2)*mslope;
    z_one = z(nnz) - risetot/2;
    z_two = z(nnz) + risetot/2;

    if z_one == z_two
        nnn = z(nnz)*ones(size(mmm));
    else
        nnn = z_one:-(z_one-z_two)/(length(mmm)-1):z_two;
    end

    dirdir_g = [mmm(1)-mmm(end),ppp(1)-ppp(end),nnn(1)-nnn(end)];
    dirdir_g = dirdir_g/norm(dirdir_g);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %interpolated quantities along 3d line

    % binterped_x = interp3(x,y,z,bx,mmm,ppp,nnn);
    % binterped_y = interp3(x,y,z,by,mmm,ppp,nnn);
    % binterped_z = interp3(x,y,z,bz,mmm,ppp,nnn);
    % jinterped_x = interp3(x,y,z,jx,mmm,ppp,nnn);
    % jinterped_y = interp3(x,y,z,jy,mmm,ppp,nnn);
    % jinterped_z = interp3(x,y,z,jz,mmm,ppp,nnn);
    % jpar_interped = interp3(x,y,z,j_par,mmm,ppp,nnn);
    % jperp_interped = interp3(x,y,z,j_perp,mmm,ppp,nnn);
    % vinterped_x = interp3(x,y,z,sx./rho,mmm,ppp,nnn);
    % vinterped_y = interp3(x,y,z,sy./rho,mmm,ppp,nnn);
    % vinterped_z = interp3(x,y,z,sz./rho,mmm,ppp,nnn);
    % rho_interped = interp3(x,y,z,rho,mmm,ppp,nnn);
    % p_interped = interp3(x,y,z,p,mmm,ppp,nnn);
    % epar_interped = interp3(x,y,z,e_par,mmm,ppp,nnn);

    binterped_x = interp3(x,y,z,permute(bx,[2 1 3]),mmm,ppp,nnn);
    binterped_y = interp3(x,y,z,permute(by,[2 1 3]),mmm,ppp,nnn);
    binterped_z = interp3(x,y,z,permute(bz,[2 1 3]),mmm,ppp,nnn);
    jinterped_x = interp3(x,y,z,permute(jx,[2 1 3]),mmm,ppp,nnn);
    jinterped_y = interp3(x,y,z,permute(jy,[2 1 3]),mmm,ppp,nnn);
    jinterped_z = interp3(x,y,z,permute(jz,[2 1 3]),mmm,ppp,nnn);
    jpar_interped = interp3(x,y,z,permute(j_par,[2 1 3]),mmm,ppp,nnn);
    jperp_interped = interp3(x,y,z,permute(j_perp,[2 1 3]),mmm,ppp,nnn);
    vinterped_x = interp3(x,y,z,permute(sx./rho,[2 1 3]),mmm,ppp,nnn);
    vinterped_y = interp3(x,y,z,permute(sy./rho,[2 1 3]),mmm,ppp,nnn);
    vinterped_z = interp3(x,y,z,permute(sz./rho,[2 1 3]),mmm,ppp,nnn);
    rho_interped = interp3(x,y,z,permute(rho,[2 1 3]),mmm,ppp,nnn);
    p_interped = interp3(x,y,z,permute(p,[2 1 3]),mmm,ppp,nnn);
    epar_interped = interp3(x,y,z,permute(e_par,[2 1 3]),mmm,ppp,nnn);

    bis = sqrt(binterped_x.^2+binterped_y.^2+binterped_z.^2);
    jis = sqrt(jinterped_x.^2+jinterped_y.^2+jinterped_z.^2);
    vis = sqrt(vinterped_x.^2+vinterped_y.^2+vinterped_z.^2);

    jdb = (jinterped_x.*binterped_x+jinterped_y.*binterped_y+...
        jinterped_z.*binterped_z)./(bis.*jis);

    %index max j
    imjdbs = find(abs(epar_interped) == max(abs(epar_interped)));
    %%grid cells away from "middle" to perform averages
    u_jdb = imjdbs+bl;
    d_jdb = imjdbs-bl;

    bxs = binterped_x(d_jdb:u_jdb);
    bys = binterped_y(d_jdb:u_jdb);
    bzs = binterped_z(d_jdb:u_jdb);
    vxs = vinterped_x(d_jdb:u_jdb);
    vys = vinterped_y(d_jdb:u_jdb);
    vzs = vinterped_z(d_jdb:u_jdb);

%     bxs = binterped_x;
%     bys = binterped_y;
%     bzs = binterped_z;
%     vxs = vinterped_x;
%     vys = vinterped_y;
%     vzs = vinterped_z;

    [rb1,rb2,rv1,rv2] = riemann_problemator(bxs(1),bxs(end),bys(1),bys(end),bzs(1),bzs(end),vxs(1),vys(1),vzs(1),vxs(end),vys(end),vzs(end))

    % bxs(1)
    % bxs(end)
    % bys(1)
    % bys(end)
    % bzs(1)
    % bzs(end)
    % vxs(1)
    % vys(1)
    % vzs(1)
    % vxs(end)
    % vys(end)
    % vzs(end)

    %do the hard work
    %V_DHT = get_DHT_frame(bxs,bys,bzs,vxs,vys,vzs);
    %i changed the 3 lines below
    %vinterped_x = vinterped_x - V_DHT(1);
    %vinterped_y = vinterped_y - V_DHT(2);
    %vinterped_z = vinterped_z - V_DHT(3);
    % [e_par,j_par,j_perp,ex,ey,ez] = get_j2(nx,ny,nz,res,jx,jy,jz,bx,by,bz,vx,vy,vz);
    % edht_interped_x = interp3(x,y,z,ex,mmm,ppp,nnn);
    % edht_interped_y = interp3(x,y,z,ey,mmm,ppp,nnn);
    % edht_interped_z = interp3(x,y,z,ez,mmm,ppp,nnn);
    % exs = edht_interped_x(d_jdb:u_jdb);
    % eys = edht_interped_y(d_jdb:u_jdb);
    % ezs = edht_interped_z(d_jdb:u_jdb);
    % [dir_cp_B,dir_cp_E] = normal_dir_cp(binterped_x,binterped_y,binterped_z,...
    %     edht_interped_x,edht_interped_y,edht_interped_z,imjdbs,u_jdb,d_jdb);
    [MVA1,MVA2,MVA3] = normal_dir_var(bxs,bys,bzs);

    brotated_x = zeros(size(binterped_x));
    brotated_y = zeros(size(binterped_x));
    brotated_z = zeros(size(binterped_x));
    for kle = 1:length(brotated_x)
        rotnew = [MVA1(:),MVA2(:),MVA3(:)]*[binterped_x(kle);binterped_y(kle);binterped_z(kle)];
        brotated_x(kle) = rotnew(1);
        brotated_y(kle) = rotnew(2);
        brotated_z(kle) = rotnew(3);
    end
    figure; hold on
    plot(ppp,brotated_x,'r')
    plot(ppp,brotated_y,'g')
    plot(ppp,brotated_z,'b')
    plot([ppp(1) ppp(end)],[0 0],'k')
    %axis([ppp(end) ppp(1) -2 2])
    %axis([ppp(1) ppp(end) -2 2])
    title('variance coordinates')
    legend('B_{min}','B_{int}','B_{max}')
    xlabel('x')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dirdir5 = zeros(3,length(binterped_x));
        dirdir5(1,:) = dirdir_g(1);
        dirdir5(2,:) = dirdir_g(2);
        dirdir5(3,:) = dirdir_g(3);
    dirdir = dirdir5;
end