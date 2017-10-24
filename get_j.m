function [jx,jy,jz] = get_j(nx,ny,nz,bx,by,bz,b0x,b0y,b0z,difx,dify,difz)
    jx = zeros(nx,ny,nz);
    jy = zeros(nx,ny,nz);
    jz = zeros(nx,ny,nz);
    bjx = zeros(nx,ny,nz);
    bjy = zeros(nx,ny,nz);
    bjz = zeros(nx,ny,nz);
    for ix = 1:nx
        for iy = 1:ny
            for iz = 1:nz
                bjx(ix,iy,iz) = bx(ix,iy,iz) - b0x(ix,iy,iz);
                bjy(ix,iy,iz) = by(ix,iy,iz) - b0y(ix,iy,iz);
                bjz(ix,iy,iz) = bz(ix,iy,iz) - b0z(ix,iy,iz);
            end
        end
    end
    for ix = 2:nx-1
        for iy = 2:ny-1
            for iz = 2:nz-1
                jx(ix,iy,iz) = dify(iy)*(bjz(ix,iy+1,iz)-bjz(ix,iy-1,iz)) -...
                    difz(iz)*(bjy(ix,iy,iz+1)-bjy(ix,iy,iz-1));
                jy(ix,iy,iz) = difz(iz)*(bjx(ix,iy,iz+1)-bjx(ix,iy,iz-1)) -...
                    difx(ix)*(bjz(ix+1,iy,iz)-bjz(ix-1,iy,iz));
                jz(ix,iy,iz) = difx(ix)*(bjy(ix+1,iy,iz)-bjy(ix-1,iy,iz)) -...
                    dify(iy)*(bjx(ix,iy+1,iz)-bjx(ix,iy-1,iz));
            end
        end
    end
end