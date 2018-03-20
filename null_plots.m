%RUN NULL TEST FIRST! OR ELSE...

%final argument determines what P is
% 1 parallel electric field
% 2 field aligned current
% 3 perpendicular current
%[jx,jy,jz,P] = get_j(nx,ny,nz,bx,by,bz,b0x,b0y,b0z,res,difx,dify,difz,0);
nus = nulls;
for i = 1:length(nulls)-1
    for j = 1:length(nulls)
        if norm(nus(i,:) - nus(j,:)) < 1 && i ~= j
            nus(j,:) = [0,0,0];
        end
    end
end

figure
eps = 0.00;
for i = 1:length(nus)
    if nus(i,1) < x(end-1)-eps && nus(i,1) > x(2)+eps && nus(i,2) < y(end-1)-eps && nus(i,2) > y(2)+eps && nus(i,3) < z(end-1)-eps && nus(i,3) > z(2)+eps && sum(nus(i,:)) > 0
        n = 300;
        spfl = ones(n,3);
        spfl(:,1) = spfl(:,1)*nus(i,1)+0.5*randn(n,1);
        spfl(:,2) = spfl(:,2)*nus(i,2)+0.5*randn(n,1);
        spfl(:,3) = spfl(:,3)*nus(i,3)+0.5*randn(n,1);
    %     zz = z(2:6);
    %     xx = x(194:200); 
    %     yy = y(102:107);
%          figure
%         plines = streamline(stream3(xx,yy,zz,bx(102:107,194:200,2:6),...
%             by(102:107,194:200,2:6),bz(102:107,194:200,2:6),spfl(:,1),spfl(:,2),spfl(:,3)));
%         hold on
%         flines = streamline(stream3(xx,yy,zz,-bx(102:107,194:200,2:6),...
%             -by(102:107,194:200,2:6),-bz(102:107,194:200,2:6),spfl(:,1),spfl(:,2),spfl(:,3)));
%         set(flines,'Color','r');
    %     isosurface(xx,yy,zz,e_ll(102:107,194:200,2:6),0.03)
    %     isosurface(xx,yy,zz,e_ll(102:107,194:200,2:6),-0.03)
    %     isosurface(xx,yy,zz,e_ll(102:107,194:200,2:6),0.05)
    %     isosurface(xx,yy,zz,e_ll(102:107,194:200,2:6),-0.05)

%     plines = streamline(stream3(x(2:end-1),y(2:end-1),z(2:end-1),bx(2:end-1,2:end-1,2:end-1),...
%         by(2:end-1,2:end-1,2:end-1),bz(2:end-1,2:end-1,2:end-1),spfl(:,2),spfl(:,1),spfl(:,3),[0.01,1000]));
%     hold on
%     flines = streamline(stream3(x(2:end-1),y(2:end-1),z(2:end-1),-bx(2:end-1,2:end-1,2:end-1),...
%         -by(2:end-1,2:end-1,2:end-1),-bz(2:end-1,2:end-1,2:end-1),spfl(:,2),spfl(:,1),spfl(:,3),[0.01,1000]));
%     set(flines,'Color','r');

        plines = streamline(stream3(x(2:end-1),y(2:end-1),z(2:end-1),bx(2:end-1,2:end-1,2:end-1),...
            by(2:end-1,2:end-1,2:end-1),bz(2:end-1,2:end-1,2:end-1),spfl(:,2),spfl(:,1),spfl(:,3)));
        hold on
        flines = streamline(stream3(x(2:end-1),y(2:end-1),z(2:end-1),-bx(2:end-1,2:end-1,2:end-1),...
            -by(2:end-1,2:end-1,2:end-1),-bz(2:end-1,2:end-1,2:end-1),spfl(:,2),spfl(:,1),spfl(:,3)));
        set(flines,'Color','r');
    end
end
%     isosurface(x(2:end-1),y(2:end-1),z(2:end-1),rho(2:end-1,2:end-1,2:end-1),0.1)
%     a=camlight('right');
%    pcolor(x(2:end-1),y(2:end-1),bx(2:end-1,2:end-1,1).^2+by(2:end-1,2:end-1,1).^2)
%    shading interp
%ftvol(ftvol > 300) = 300;
%ftjpar(ftjpar < -100) = median(median(ftjpar));
%ftmass(ftmass > 500) = median(median(ftmass));

pcolor(x,y,bz(:,:,2))
shading interp
h = streamslice(x,y,bx(:,:,2),by(:,:,2));
set(h,'Color','k')
daspect([1 1 1])
colorbar
%view(3)
%title(strcat('pos net 42 time = ',num2str(time)))

