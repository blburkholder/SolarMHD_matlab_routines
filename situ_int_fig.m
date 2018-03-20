% % space slices out without the data
% figure
% % create the slices 
% jpar1 = j_par;
% jmax = max(max(j_par(:,:,35)));
% jpar1(jpar1 > jmax) = jmax;
% jmin = min(min(j_par(:,:,35)));
% jpar1(jpar1 < jmin) = jmin;
% h = slice(x,y,z,j_par/max([abs(jmin),jmax]),[],[],[35]);
% h.FaceColor = 'interp';
% h.EdgeColor = 'none';
% 
% hold on
% jpar2 = j_par;
% jmax = max(max(j_par(:,:,20)));
% jpar2(jpar2 > jmax) = jmax;
% jmin = min(min(j_par(:,:,20)));
% jpar2(jpar2 < jmin) = jmin;
% h = slice(x,y,z,jpar2/max([abs(jmin),jmax]),[],[],[20]);
% 
% jpar2 = j_par;
% jmax = max(max(j_par(:,:,50)));
% jpar2(jpar2 > jmax) = jmax;
% jmin = min(min(j_par(:,:,50)));
% jpar2(jpar2 < jmin) = jmin;
% h = slice(x,y,z,jpar2/max([abs(jmin),jmax]),[],[],[50]);
% 
% ree = 30;
% [ssx,ssy,ssz] = meshgrid(x(2:ree:end-1),y(2:ree:end-1),z(end-1));
% streamline(stream3(x,y,z,bx,by,bz,ssx,ssy,ssz))
% streamline(stream3(x,y,z,-bx,-by,-bz,ssx,ssy,ssz))
% daspect([1 1 0.3])
% 
% load('mentrop10_t.mat')
% ftjpar(ftjpar > 5) = 5;
% ftjpar(ftjpar < -5) = -5;
% jpar_cube = zeros(1001,1001,153);
% jpar_cube(:,:,end-2) = ftjpar;
% jpar_cube(:,:,end-1) = ftjpar;
% jpar_cube(:,:,end) = ftjpar;
% jmax = max(max(jpar_cube(:,:,end-1)));
% jmin = min(min(jpar_cube(:,:,end-1)));
% h = slice(xpo,ypo,z,jpar_cube/max([abs(jmin),jmax]),[],[],[z(end-1)]);
% 
% load('mentrop10.mat')
% ftjpar(ftjpar > 5) = 5;
% ftjpar(ftjpar < -5) = -5;
% jmax = max(max(ftjpar));
% jmin = min(min(ftjpar));
% h = pcolor(xpo,ypo,ftjpar/max([abs(jmin),jmax]));
% 
% shading interp
% caxis([-0.1 0.1])

ree = 15;
zzz = 120;
%[ssx,ssy,ssz] = meshgrid(x(60:ree:219),y(195:ree:233),z(zzz));
[ssx,ssy,ssz] = meshgrid(x(2:ree:end-1),y(2:ree:end-1),z(zzz));
flines1 = perio_stream3(x,y,z,bx,by,bz,ssx,ssy,ssz);
flines2 = perio_stream33(x,y,z,bx,by,bz,ssx,ssy,ssz);
%flines1 = stream3(x,y,z,bx,by,bz,ssx,ssy,ssz);
%flines2 = stream3(x,y,z,-bx,-by,-bz,ssx,ssy,ssz);

conj_x = zeros(2,length(flines1));
conj_y = zeros(2,length(flines1));
conj_z = zeros(2,length(flines1));

conj_vx = zeros(2,length(flines1));
conj_vy = zeros(2,length(flines1)); 

for i = 1:length(flines1)
    conj_x(:,i) = [flines1{i}(end,1),flines2{i}(end,1)];
    conj_y(:,i) = [flines1{i}(end,2),flines2{i}(end,2)];
    conj_z(:,i) = [flines1{i}(end,3),flines2{i}(end,3)];
end

conj_vx(1,:) = interp3(x,y,z,sx./rho,conj_x(1,:),conj_y(1,:),conj_z(1,:));
conj_vx(2,:) = interp3(x,y,z,sx./rho,conj_x(2,:),conj_y(2,:),conj_z(2,:));
conj_vy(1,:) = interp3(x,y,z,sy./rho,conj_x(1,:),conj_y(1,:),conj_z(1,:));
conj_vy(2,:) = interp3(x,y,z,sy./rho,conj_x(2,:),conj_y(2,:),conj_z(2,:));

figure
pcolor(x,y,j_par(:,:,zzz));
hold on
quiver(ssx,ssy,reshape(conj_vx(1,:),size(ssx)),reshape(conj_vy(1,:),size(ssx)),'k')
quiver(ssx,ssy,reshape(conj_vx(2,:),size(ssx)),reshape(conj_vy(2,:),size(ssx)),'k')
%h = scatter(ssx(:),ssy(:));
%cm = colormap;
%h.CData = zeros(size(h.XData));

shading interp
daspect([1 1 1])
 title(['z=',num2str(z(zzz))])
