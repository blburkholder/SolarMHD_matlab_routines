%RUN NULL TEST FIRST! OR ELSE...

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
% 
% e_ll = zeros(nx,ny,nz);
% for ix = 2:nx-1
%     for iy = 2:ny-1
%         for iz = 2:nz-1
%             e_ll(ix,iy,iz) = res(ix,iy,iz)*(jx(ix,iy,iz)*bx(ix,iy,iz) +...
%                 jy(ix,iy,iz)*by(ix,iy,iz) + jz(ix,iy,iz)*bz(ix,iy,iz))/...
%                 sqrt(bx(ix,iy,iz)^2 + by(ix,iy,iz)^2 + bz(ix,iy,iz)^2); 
%         end
%     end
% end

figure
eps = 0.00;
for i = 1:length(nulls)
    if nulls(i,1) < x(end-1)-eps && nulls(i,1) > x(2)+eps && nulls(i,2) < y(end-1)-eps && nulls(i,2) > y(2)+eps && nulls(i,3) < z(end-1)-eps && nulls(i,3) > z(2)+eps
        n = 100;
        spfl = ones(n,3);
        spfl(:,1) = spfl(:,1)*nulls(i,1)+0.1*randn(n,1)-0.05;
        spfl(:,2) = spfl(:,2)*nulls(i,2)+0.1*randn(n,1)-0.05;
        spfl(:,3) = spfl(:,3)*nulls(i,3)+0.1*randn(n,1)-0.05;
    %     zz = z(2:6);
    %     xx = x(194:200); 
    %     yy = y(102:107);
%          figure
    %     plines = streamline(stream3(xx,yy,zz,by(102:107,194:200,2:6),...
    %         bx(102:107,194:200,2:6),bz(102:107,194:200,2:6),spfl(:,1),spfl(:,2),spfl(:,3)));
    %     hold on
    %     flines = streamline(stream3(xx,yy,zz,-by(102:107,194:200,2:6),...
    %         -bx(102:107,194:200,2:6),-bz(102:107,194:200,2:6),spfl(:,1),spfl(:,2),spfl(:,3)));
    %     set(flines,'Color','r');
    %     isosurface(xx,yy,zz,e_ll(102:107,194:200,2:6),0.03)
    %     isosurface(xx,yy,zz,e_ll(102:107,194:200,2:6),-0.03)
    %     isosurface(xx,yy,zz,e_ll(102:107,194:200,2:6),0.05)
    %     isosurface(xx,yy,zz,e_ll(102:107,194:200,2:6),-0.05)

%     plines = streamline(stream3(y(2:end-1),x(2:end-1),z(2:end-1),by(2:end-1,2:end-1,2:end-1),...
%         bx(2:end-1,2:end-1,2:end-1),bz(2:end-1,2:end-1,2:end-1),spfl(:,1),spfl(:,2),spfl(:,3),[0.01,100]));
%     hold on
%     flines = streamline(stream3(y(2:end-1),x(2:end-1),z(2:end-1),-by(2:end-1,2:end-1,2:end-1),...
%         -bx(2:end-1,2:end-1,2:end-1),-bz(2:end-1,2:end-1,2:end-1),spfl(:,1),spfl(:,2),spfl(:,3),[0.01,100]));
%   set(flines,'Color','r');

        plines = streamline(stream3(y(2:end-1),x(2:end-1),z(2:end-1),by(2:end-1,2:end-1,2:end-1),...
            bx(2:end-1,2:end-1,2:end-1),bz(2:end-1,2:end-1,2:end-1),spfl(:,1),spfl(:,2),spfl(:,3)));
        hold on
        flines = streamline(stream3(y(2:end-1),x(2:end-1),z(2:end-1),-by(2:end-1,2:end-1,2:end-1),...
            -bx(2:end-1,2:end-1,2:end-1),-bz(2:end-1,2:end-1,2:end-1),spfl(:,1),spfl(:,2),spfl(:,3)));
        set(flines,'Color','r');
    end
end
%     isosurface(y(2:end-1),x(2:end-1),z(2:end-1),rho(2:end-1,2:end-1,2:end-1),0.1)
%     a=camlight('right');
%    pcolor(y(2:end-1),x(2:end-1),bx(2:end-1,2:end-1,1).^2+by(2:end-1,2:end-1,1).^2)
%    shading interp
ftvol(ftvol > 300) = 300;
%ftjpar(ftjpar < -100) = median(median(ftjpar));
%ftmass(ftmass > 500) = median(median(ftmass));

pcolor(xpo,ypo,ftvol)
shading interp
daspect([1 1 1])
colorbar
%title(strcat('pos net 42 time = ',num2str(time)))

