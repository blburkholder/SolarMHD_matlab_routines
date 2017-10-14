function [jx,jy,jz,P] = get_j(nx,ny,nz,bx,by,bz,b0x,b0y,b0z,res,difx,dify,difz,tag)
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

    %E parallel
    if tag == 1
        P = zeros(nx,ny,nz);
        for ix = 2:nx-1
            for iy = 2:ny-1
                for iz = 2:nz-1
                    P(ix,iy,iz) = res(ix,iy,iz)*(jx(ix,iy,iz)*bx(ix,iy,iz) +...
                        jy(ix,iy,iz)*by(ix,iy,iz) + jz(ix,iy,iz)*bz(ix,iy,iz))/...
                        sqrt(bx(ix,iy,iz)^2 + by(ix,iy,iz)^2 + bz(ix,iy,iz)^2); 
                end
            end
        end
    %field aligned current
    elseif tag == 2
        P = zeros(nx,ny,nz);
        for ix = 2:nx-1
            for iy = 2:ny-1
                for iz = 2:nz-1
                    BBB = sqrt(bx(ix,iy,iz)^2 + by(ix,iy,iz)^2 + bz(ix,iy,iz)^2);
                    dir_bx = bx(ix,iy,iz)/BBB;
                    dir_by = by(ix,iy,iz)/BBB;
                    dir_bz = bz(ix,iy,iz)/BBB;
                    P(ix,iy,iz) = jx(ix,iy,iz)*dir_bx + jy(ix,iy,iz)*dir_by + jz(ix,iy,iz)*dir_bz;
                end
            end
        end
    %perpendicular current
    elseif tag == 3
        P = zeros(nx,ny,nz);
        for ix = 2:nx-1
            for iy=2:ny-1
                for iz = 2:nz-1
                    P(ix,iy,iz) = sqrt((jy(ix,iy,iz)*bz(ix,iy,iz) -...
                        jz(ix,iy,iz)*by(ix,iy,iz))^2 + (jz(ix,iy,iz)*bx(ix,iy,iz) -...
                        jx(ix,iy,iz)*bz(ix,iy,iz))^2 + (jx(ix,iy,iz)*by(ix,iy,iz) -...
                        jy(ix,iy,iz)*bx(ix,iy,iz))^2)/sqrt((bx(ix,iy,iz))^2 +...
                        (by(ix,iy,iz))^2 + (bz(ix,iy,iz))^2);
                end
            end
        end
    else
        P = 0;
    end
end