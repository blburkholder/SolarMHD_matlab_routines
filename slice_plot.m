%clear all

%magnetic field lines solid color

re = 5;
% figure
%hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%from top of box%%%%%%%%%%%%%%%%%%%%%%%%
%{
[cx,cy,cz] = meshgrid(min(x):res:max(x),min(y):res:max(y),max(z));
streamline(stream3(y,x,z,by,bx,bz,cx(:),cy(:),cz(:)));
streamline(stream3(y,x,z,-by,-bx,-bz,cx(:),cy(:),cz(:)));
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%from bottom of box%%%%%%%%%%%%%%%%%%%%%%%%
%{
[cx,cy,cz] = meshgrid(min(x):re:max(x),min(y):re:max(y),min(z));
h = streamline(stream3(y,x,z,by,bx,bz,cx(:),cy(:),cz(:)));
set(h,'Color','g')
h = streamline(stream3(y,x,z,-by,-bx,-bz,cx(:),cy(:),cz(:)));
set(h,'Color','m')
%}

%magnetic field lines colored wrt height
%{
res = 5;
%resolution
figure
hold on
colormap(jet)
%%%%%%%%%%%%%%%%%%%%%%%%%%%from top of box%%%%%%%%%%%%%%%%%%%%%%%%
[cx,cy,cz] = meshgrid(min(x):res:max(x),min(y):res:max(y),max(z));
color_stream_pos = stream3(y,x,z,by,bx,bz,cx(:),cy(:),cz(:));
color_stream_neg = stream3(y,x,z,-by,-bx,-bz,cx(:),cy(:),cz(:));
maximum = max(z);
for i = 1:length(color_stream_pos)
        try
        patch([(color_stream_pos{i}(:,1))' nan],[(color_stream_pos{i}(:,2))' nan],...
            [(color_stream_pos{i}(:,3))' nan],[abs((color_stream_pos{i}(:,3)/maximum)') nan],...
            'EdgeColor','interp','FaceColor','none');
        catch

        end
end

for i = 1:length(color_stream_neg)
        try
        patch([(color_stream_neg{i}(:,1))' nan],[(color_stream_neg{i}(:,2))' nan],...
            [(color_stream_neg{i}(:,3))' nan],[1-abs((color_stream_neg{i}(:,3)/maximum)') nan],...
            'EdgeColor','interp','FaceColor','none');
        catch

        end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%bottom%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
colormap(jet)
[cx,cy,cz] = meshgrid(min(x):res:max(x),min(y):res:max(y),min(z));
color_stream_pos = stream3(y,x,z,by,bx,bz,cx(:),cy(:),cz(:));
color_stream_neg = stream3(y,x,z,-by,-bx,-bz,cx(:),cy(:),cz(:));
for i = 1:length(color_stream_pos)
        try
        patch([(color_stream_pos{i}(:,1))' nan],[(color_stream_pos{i}(:,2))' nan],...
            [(color_stream_pos{i}(:,3))' nan],[abs((color_stream_pos{i}(:,3)/maximum)') nan],...
            'EdgeColor','interp','FaceColor','none');
        catch

        end
end

for i = 1:length(color_stream_neg)
        try
        patch([(color_stream_neg{i}(:,1))' nan],[(color_stream_neg{i}(:,2))' nan],...
            [(color_stream_neg{i}(:,3))' nan],[1-abs((color_stream_neg{i}(:,3)/maximum)') nan],...
            'EdgeColor','interp','FaceColor','none');
        catch

        end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%isosurfaces of pressure
%or whatever else by replacing p
%h = patch(isosurface(y,x,z,p,0.1))
%set(h,'FaceAlpha',0.1)
%set(h,'EdgeAlpha',0.1)
%set(h,'EdgeColor','b')

%slices
%{
h = slice(x,y,z,bx.^2+by.^2+bz.^2,[],[],[min(z)]);
colormap(jet)
for i = 1:length(h)
    h(i).FaceColor = 'interp';
    h(i).FaceAlpha = 1;
    h(i).EdgeColor = 'none';
    h(i).DiffuseStrength = 0.8;
end
quiver(x(1:10:end),y(1:10:end),sy(1:10:end,1:10:end,1),sx(1:10:end,1:10:end,1),'Color','w')
%}

%in plane 3D vector field
%{
whereslice = 1;
figure
%pcolor(v(:,:,whereslice))
pcolor(p(:,:,whereslice))
shading interp
colorbar
hold on
%[cx,cy] = meshgrid(1:10:93,1:10:93);
quiver(1:10:259,1:10:259,sx(1:10:end,1:10:end,whereslice),sy(1:10:end,1:10:end,whereslice),'Color','w')
%}

%flux tube volumes
% ex = 0:0.3625:93.5250;
% wy = 0:0.3625:93.5250;
%{
ex = 0:93;
wy = 0:93;
[sx,sy,sz] = meshgrid(ex,wy,1);
ft = stream3(y,x,z,by,bx,bz,sx(:),sy(:),sz(:));
ftv = zeros(length(ex)^2,1);
for i = 1:length(ft)
    if ~isempty(ft{i})
        xval = interp3(x,y,z,bx,ft{i}(:,1),ft{i}(:,2),ft{i}(:,3));
        yval = interp3(x,y,z,by,ft{i}(:,1),ft{i}(:,2),ft{i}(:,3));
        zval = interp3(x,y,z,bz,ft{i}(:,1),ft{i}(:,2),ft{i}(:,3));
        for j = 2:length(ft{i}(:,1))
            ftv(i) = ftv(i) + sqrt((ft{i}(j,1)-ft{i}(j-1,1)).^2 + (ft{i}(j,2)-ft{i}(j-1,2)).^2 + (ft{i}(j,3)-ft{i}(j-1,3)).^2)/...
                sqrt(xval(j)^2+yval(j)^2+zval(j)^2);
        end
    end
end
ftv = reshape(ftv,[length(ex),length(ex)]);
figure
h = pcolor(ftv);
set(h,'edgecolor','none');
colorbar
%}

%current_density
%%%%%%%%%%%%%%%%%%%%%this is how antonius does it%%%%%%%%%%%%%%%%%%%%

% jx = zeros(nx,ny,nz);
% jy = zeros(nx,ny,nz);
% jz = zeros(nx,ny,nz);
% bjx = zeros(nx,ny,nz);
% bjy = zeros(nx,ny,nz);
% bjz = zeros(nx,ny,nz);
% for ix = 1:nx
%     for iy = 1:ny
%         for iz = 1:nz
%             bjx(ix,iy,iz) = bx(ix,iy,iz) - b0x(ix,iy,iz);
%             bjy(ix,iy,iz) = by(ix,iy,iz) - b0y(ix,iy,iz);
%             bjz(ix,iy,iz) = bz(ix,iy,iz) - b0z(ix,iy,iz);
%         end
%     end
% end
% for ix = 2:nx-1
%     for iy = 2:ny-1
%         for iz = 2:nz-1
%             jx(ix,iy,iz) = dify(iy)*(bjz(ix,iy+1,iz)-bjz(ix,iy-1,iz)) -...
%                 difz(iz)*(bjy(ix,iy,iz+1)-bjy(ix,iy,iz-1));
%             jy(ix,iy,iz) = difz(iz)*(bjx(ix,iy,iz+1)-bjx(ix,iy,iz-1)) -...
%                 difx(ix)*(bjz(ix+1,iy,iz)-bjz(ix-1,iy,iz));
%             jz(ix,iy,iz) = difx(ix)*(bjy(ix+1,iy,iz)-bjy(ix-1,iy,iz)) -...
%                 dify(iy)*(bjx(ix,iy+1,iz)-bjx(ix,iy-1,iz));
%         end
%     end
% end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
jx1 = zeros(nx,ny,nz);
jy1 = zeros(nx,ny,nz);
jz1 = zeros(nx,ny,nz);
for ix = 2:nx-1
    for iy = 2:ny-1
        for iz = 2:nz-1
            jx1(ix,iy,iz) = dify(iy)*(bz(ix,iy+1,iz)-bz(ix,iy-1,iz)) -...
                difz(iz)*(by(ix,iy,iz+1)-by(ix,iy,iz-1));
            jy1(ix,iy,iz) = difz(iz)*(bx(ix,iy,iz+1)-bx(ix,iy,iz-1)) -...
                difx(ix)*(bz(ix+1,iy,iz)-bz(ix-1,iy,iz));
            jz1(ix,iy,iz) = difx(ix)*(by(ix+1,iy,iz)-by(ix-1,iy,iz)) -...
                dify(iy)*(bx(ix,iy+1,iz)-bx(ix,iy-1,iz));
        end
    end
end

jx2 = zeros(nx,ny,nz);
jy2 = zeros(nx,ny,nz);
jz2 = zeros(nx,ny,nz);
for ix = 2:nx-1
    for iy = 2:ny-1
        for iz = 2:nz-1
            jx2(ix,iy,iz) = dify(iy)*(b0z(ix,iy+1,iz)-b0z(ix,iy-1,iz)) -...
                difz(iz)*(b0y(ix,iy,iz+1)-b0y(ix,iy,iz-1));
            jy2(ix,iy,iz) = difz(iz)*(b0x(ix,iy,iz+1)-b0x(ix,iy,iz-1)) -...
                difx(ix)*(b0z(ix+1,iy,iz)-b0z(ix-1,iy,iz));
            jz2(ix,iy,iz) = difx(ix)*(b0y(ix+1,iy,iz)-b0y(ix-1,iy,iz)) -...
                dify(iy)*(b0x(ix,iy+1,iz)-b0x(ix,iy-1,iz));
        end
    end
end
%}

%parallel electric fields
%{
e_ll = zeros(nx,ny,nz);
for ix = 2:nx-1
    for iy = 2:ny-1
        for iz = 2:nz-1
            e_ll(ix,iy,iz) = res(ix,iy,iz)*(jx(ix,iy,iz)*bx(ix,iy,iz) +...
                jy(ix,iy,iz)*by(ix,iy,iz) + jz(ix,iy,iz)*bz(ix,iy,iz))/...
                sqrt(bx(ix,iy,iz)^2 + by(ix,iy,iz)^2 + bz(ix,iy,iz)^2); 
        end
    end
end
%}

%specific entropy
%sp_ent = p./(rho).^(5/3);
 
%electric field

% ex = zeros(nx,ny,nz);
% ey = zeros(nx,ny,nz);
% ez = zeros(nx,ny,nz);
% for ix = 2:nx-1
%     for iy = 2:ny-1
%         for iz = 2:nz-1
%             ex(ix,iy,iz) = res(ix,iy,iz)*jx(ix,iy,iz) - (sy(ix,iy,iz)*bz(ix,iy,iz) - sz(ix,iy,iz)*by(ix,iy,iz))/rho(ix,iy,iz);
%             ey(ix,iy,iz) = res(ix,iy,iz)*jy(ix,iy,iz) - (sz(ix,iy,iz)*bx(ix,iy,iz) - sx(ix,iy,iz)*bz(ix,iy,iz))/rho(ix,iy,iz);
%             ez(ix,iy,iz) = res(ix,iy,iz)*jz(ix,iy,iz) - (sx(ix,iy,iz)*by(ix,iy,iz) - sy(ix,iy,iz)*bx(ix,iy,iz))/rho(ix,iy,iz);
%         end
%     end
% end
% %}
% 
% %E cross B perturbation
% ExB_px = zeros(nx,ny,nz);
% ExB_py = zeros(nx,ny,nz);
% ExB_pz = zeros(nx,ny,nz);
% for ix = 2:nx-1
%     for iy = 2:ny-1
%         for iz = 2:nz-1
%             ExB_px(ix,iy,iz) = ey(ix,iy,iz)*(bz(ix,iy,iz)-b0z(ix,iy,iz)) - ez(ix,iy,iz)*(by(ix,iy,iz)-b0y(ix,iy,iz));
%             ExB_py(ix,iy,iz) = ez(ix,iy,iz)*(bx(ix,iy,iz)-b0x(ix,iy,iz)) - ex(ix,iy,iz)*(bz(ix,iy,iz)-b0z(ix,iy,iz)); 
%             ExB_pz(ix,iy,iz) = ex(ix,iy,iz)*(by(ix,iy,iz)-b0y(ix,iy,iz)) - ey(ix,iy,iz)*(bx(ix,iy,iz)-b0x(ix,iy,iz));
%         end
%     end
% end


%generator?
%edj = ex.*jx + ey.*jy + ez.*jz;

%Field aligned current
%{
fac = zeros(nx,ny,nz);
%fac1 = zeros(nx,ny,nz);
%fac2 = zeros(nx,ny,nz);
for ix = 2:nx-1
    for iy = 2:ny-1
        for iz = 2:nz-1
            BBB = sqrt(bx(ix,iy,iz)^2 + by(ix,iy,iz)^2 + bz(ix,iy,iz)^2);
            dir_bx = bx(ix,iy,iz)/BBB;
            dir_by = by(ix,iy,iz)/BBB;
            dir_bz = bz(ix,iy,iz)/BBB;
            fac(ix,iy,iz) = jx(ix,iy,iz)*dir_bx + jy(ix,iy,iz)*dir_by + jz(ix,iy,iz)*dir_bz;
            %fac1(ix,iy,iz) = jx1(ix,iy,iz)*dir_bx + jy1(ix,iy,iz)*dir_by + jz1(ix,iy,iz)*dir_bz;
            %fac2(ix,iy,iz) = jx2(ix,iy,iz)*dir_bx + jy2(ix,iy,iz)*dir_by + jz2(ix,iy,iz)*dir_bz;
        end
    end
end
%}

%%%%%%%%%%%%%%%can only do one of these at a time%%%%%%%%%%%%%%%%%%%%%%%%%
%this one for E parallel
% e_ll = ExB_pz.^2 + ExB_py.^2 + ExB_px.^2;
% iz = 2;
% f = figure;
% ax = axes('Parent',f);
% h = pcolor(ax,x,y,e_ll(:,:,iz)');
% shading interp
% colorbar
% b = uicontrol('Parent',f,'Style','slider','Position',[81,0,419,23],...
%               'value',iz, 'min',2, 'max',nz-1);
% set(b, 'SliderStep', [1/(nz-1) , 10/(nz-1) ]);
% b.Callback  = @(objHandle,~)  update_plot(x,y,e_ll,get(objHandle,'Value'))
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iz = 2;
f = figure;
ax = axes('Parent',f);
r = sx(1:re:end,1:re:end,iz);
%r = reshape(r(:), [nx ny]);
%h = pcolor(ax,y,x,log10(r'));
h = pcolor(ax,y(1:re:end),x(1:re:end),r);
hold on
shading interp
colorbar
b = uicontrol('Parent',f,'Style','slider','Position',[81,0,419,23],...
              'value',iz, 'min',2, 'max',nz-1);
set(b, 'SliderStep', [1/(nz-1) , 10/(nz-1) ]);
b.Callback  = @(objHandle,~)  update_plot4(y(1:re:end),x(1:re:end),sx(1:re:end,1:re:end,:),get(objHandle,'Value'),ny,nx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%{
bzz = bz;
bzz(find(bzz<0)) = 0;
bzz(find(bzz>0)) = 1;
iz = 2;
f = figure;
ax = axes('Parent',f);
h = pcolor(ax,x,y,bzz(:,:,iz)');
shading interp
colorbar
b = uicontrol('Parent',f,'Style','slider','Position',[81,0,419,23],...
              'value',iz, 'min',2, 'max',nz-1);
set(b, 'SliderStep', [1/(nz-1) , 10/(nz-1) ]);
b.Callback  = @(objHandle,~)  update_plot(x,y,bzz,get(objHandle,'Value'))
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this one for field aligned currents
%{
iz = 2;
l = figure;
ax1 = axes('Parent',l);
h1 = pcolor(ax1,x,y,fac(:,:,iz)');
shading interp
colorbar
b1 = uicontrol('Parent',l,'Style','slider','Position',[81,0,419,23],...
              'value',iz, 'min',2, 'max',nz-1);
b1.Callback  = @(objHandle,~)  update_plot(x,y,fac,floor(get(objHandle,'Value'))) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}


%{
figure
colormap(jet)
for iz = 2:nz-1
    subplot(3,1,1)
    pcolor(x,y,fac(:,:,iz));
    shading interp
    colorbar
    title(strcat('slice ',num2str(iz)));
    subplot(3,1,2)
    pcolor(x,y,fac1(:,:,iz));
    shading interp
    colorbar
    subplot(3,1,3)
    pcolor(x,y,fac2(:,:,iz));
    shading interp
    colorbar
    pause(0.01)
end
%}

