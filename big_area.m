wx = [-flipud(x(4:end));x];
wvx = [-fliplr(flipud(sx(:,3:end,:)./rho(:,3:end,:))),sx(:,2:end,:)./rho(:,2:end,:)];
wvy = [-fliplr(flipud(sy(:,3:end,:)./rho(:,3:end,:))),sy(:,2:end,:)./rho(:,2:end,:)];
wvz = [fliplr(flipud(sz(:,3:end,:)./rho(:,3:end,:))),sz(:,2:end,:)./rho(:,2:end,:)];
wbx = [-fliplr(flipud(bx(:,3:end,:))),bx(:,2:end,:)];
wby = [-fliplr(flipud(by(:,3:end,:))),by(:,2:end,:)];
wbz = [fliplr(flipud(bz(:,3:end,:))),bz(:,2:end,:)];
wrho = [fliplr(flipud(rho(:,3:end,:))),rho(:,2:end,:)];

nn = nulls(end-4,:);
nflines = 5;
ss = 0.1;
vvvec = (-ss:2*ss/nflines:ss);
[xx1,yy1,zz1] = meshgrid(nn(2)+vvvec,nn(1)+vvvec,nn(3)+vvvec);

plines1 = stream3(wx,y,z,wbx,wby,wbz,xx1,yy1,zz1,[0.01,100]); %spine
plines2 = stream3(wx,y,z,-wbx,-wby,-wbz,xx1,yy1,zz1,[0.01,100]); %fan
%plines1 = stream3(wx,y,z,wbx,wby,wbz,xx1,yy1,zz1); %spine
%plines2 = stream3(wx,y,z,-wbx,-wby,-wbz,xx1,yy1,zz1); %fan

%figure; hold on;
xyzfin = zeros(length(plines2),3);
 for i = 1:length(plines2)
     xyzfin(i,:) = plines2{i}(end,:); 
 end
%scatter3(xyzfin(:,1),xyzfin(:,2),xyzfin(:,3))
%h = streamline(plines1);
%set(h,'Color','r')
[normal,~,point] = affine_fit(xyzfin);
scatter3(point(1),point(2),point(3))
numflines = 50;
%quiver3(point(1),point(2),point(3),normal(1)/3,normal(2)/3,normal(3)/3,'r','linewidth',2)

[X,Y] = meshgrid(min(xyzfin(:,1)):(max(xyzfin(:,1))-min(xyzfin(:,1)))/numflines:max(xyzfin(:,1)),...
    min(xyzfin(:,2)):(max(xyzfin(:,2))-min(xyzfin(:,2)))/numflines:max(xyzfin(:,2)));
Z = -(normal(1)/normal(3)*X+normal(2)/normal(3)*Y-dot(normal,point)/normal(3));
%surf(X,Y, - (normal(1)/normal(3)*X+normal(2)/normal(3)*Y-dot(normal,point)/normal(3)),'facecolor','red','facealpha',0.5);
%streamline(plines2)

j_plines1 = stream3(wx,y,z,wbx,wby,wbz,X(:),Y(:),Z(:)); %spine
j_plines2 = stream3(wx,y,z,-wbx,-wby,-wbz,X(:),Y(:),Z(:)); %fan
e_plines1 = j_plines1; e_plines2 = j_plines2;

if ~exist('j_par')
 [jx,jy,jz] = get_j(nx,ny,nz,bx,by,bz,b0x,b0y,b0z,difx,dify,difz);
 [e_par,j_par,j_perp,ex,ey,ez] = get_j2(nx,ny,nz,res,jx,jy,jz,bx,by,bz,sx./rho,sy./rho,sz./rho);
end
%j_tot = sqrt(j_par.^2 + j_perp.^2);
%wj_tot = [fliplr(flipud(j_tot(:,3:end,:))),j_tot(:,2:end,:)];
wj_par = [fliplr(flipud(j_par(:,3:end,:))),j_par(:,2:end,:)];
we_par = [fliplr(flipud(e_par(:,3:end,:))),e_par(:,2:end,:)];

b1 = 0; b2 = 0;
q1 = 0; q2 = 0;
r1 = 0; r2 = 0;
p2_jpar = zeros(length(j_plines2),1);
p1_jpar = zeros(length(j_plines2),1);
p2_epar = zeros(length(j_plines2),1);
p1_epar = zeros(length(j_plines2),1);
for i = 1:length(j_plines2)
    if ~isempty(j_plines2{i})
        q1 = interp3(wx(2:end-1),y(2:end-1),z(2:end-1),abs(wj_par(2:end-1,2:end-1,2:end-1)),j_plines2{i}(:,1),j_plines2{i}(:,2),j_plines2{i}(:,3));
        r1 = interp3(wx(2:end-1),y(2:end-1),z(2:end-1),abs(wj_par(2:end-1,2:end-1,2:end-1)),j_plines1{i}(:,1),j_plines1{i}(:,2),j_plines1{i}(:,3));
        if max(max([q1;r1])) > b1
            b1 = max(max([q1;r1]));
        end

        q2 = interp3(wx(2:end-1),y(2:end-1),z(2:end-1),abs(we_par(2:end-1,2:end-1,2:end-1)),j_plines2{i}(:,1),j_plines2{i}(:,2),j_plines2{i}(:,3));
        r2 = interp3(wx(2:end-1),y(2:end-1),z(2:end-1),abs(we_par(2:end-1,2:end-1,2:end-1)),j_plines1{i}(:,1),j_plines1{i}(:,2),j_plines1{i}(:,3));
        if max(max([q2;r2])) > b2
            b2 = max(max([q2;r2]));
        end

        j_plines2{i}(:,4) = q1;
        j_plines1{i}(:,4) = r1;
        e_plines2{i}(:,4) = q2;
        e_plines1{i}(:,4) = r2;
    end
i
end
maximum1 = b1;
maximum2 = b2;

figure; hold on;
for i = 1:length(j_plines2)
    try
    patch([(j_plines2{i}(:,1))' nan],[(j_plines2{i}(:,2))' nan],...
        [(j_plines2{i}(:,3))' nan],[(j_plines2{i}(:,4)/maximum1)' nan],...
        'EdgeColor','interp','FaceColor','none');
    catch

    end
    try
    patch([(j_plines1{i}(:,1))' nan],[(j_plines1{i}(:,2))' nan],...
        [(j_plines1{i}(:,3))' nan],[(j_plines1{i}(:,4)/maximum1)' nan],...
        'EdgeColor','interp','FaceColor','none');
    catch

    end
end
daspect([1 1 1])
view(3)
title('j_{||}')
colormap(hsv)
colorbar

figure; hold on;
for i = 1:length(e_plines2)
    try
    patch([(e_plines2{i}(:,1))' nan],[(e_plines2{i}(:,2))' nan],...
        [(e_plines2{i}(:,3))' nan],[(e_plines2{i}(:,4)/maximum2)' nan],...
        'EdgeColor','interp','FaceColor','none');
    catch

    end
    try
    patch([(e_plines1{i}(:,1))' nan],[(e_plines1{i}(:,2))' nan],...
        [(e_plines1{i}(:,3))' nan],[(e_plines2{i}(:,4)/maximum2)' nan],...
        'EdgeColor','interp','FaceColor','none');
    catch

    end
end
daspect([1 1 1])
view(3)
title('E_{||}')
colormap(hsv)
colorbar


plot3(mmm,ppp,nnn,'LineWidth',3,'Color','k')
%isosurface(wx,y,z,wj_tot,0.01)
