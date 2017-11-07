function [e_par,j_par,j_perp,ex,ey,ez] = get_j2(nx,ny,nz,res,jx,jy,jz,bx,by,bz,sx,sy,sz,rho)
  
    %parallel electric field
    e_par = zeros(nx,ny,nz);
    for ix = 2:nx-1
        for iy = 2:ny-1
            for iz = 2:nz-1
                e_par(ix,iy,iz) = res(ix,iy,iz)*(jx(ix,iy,iz)*bx(ix,iy,iz) +...
                    jy(ix,iy,iz)*by(ix,iy,iz) + jz(ix,iy,iz)*bz(ix,iy,iz))/...
                    sqrt(bx(ix,iy,iz)^2 + by(ix,iy,iz)^2 + bz(ix,iy,iz)^2); 
            end
        end
    end

    %field aligned current
    j_par = zeros(nx,ny,nz);
    for ix = 2:nx-1
        for iy = 2:ny-1
            for iz = 2:nz-1
                BBB = sqrt(bx(ix,iy,iz)^2 + by(ix,iy,iz)^2 + bz(ix,iy,iz)^2);
                dir_bx = bx(ix,iy,iz)/BBB;
                dir_by = by(ix,iy,iz)/BBB;
                dir_bz = bz(ix,iy,iz)/BBB;
                j_par(ix,iy,iz) = jx(ix,iy,iz)*dir_bx + jy(ix,iy,iz)*dir_by + jz(ix,iy,iz)*dir_bz;
            end
        end
    end

    %perpendicular current
    j_perp = zeros(nx,ny,nz);
    for ix = 2:nx-1
        for iy=2:ny-1
            for iz = 2:nz-1
                j_perp(ix,iy,iz) = sqrt((jy(ix,iy,iz)*bz(ix,iy,iz) -...
                    jz(ix,iy,iz)*by(ix,iy,iz))^2 + (jz(ix,iy,iz)*bx(ix,iy,iz) -...
                    jx(ix,iy,iz)*bz(ix,iy,iz))^2 + (jx(ix,iy,iz)*by(ix,iy,iz) -...
                    jy(ix,iy,iz)*bx(ix,iy,iz))^2)/sqrt((bx(ix,iy,iz))^2 +...
                    (by(ix,iy,iz))^2 + (bz(ix,iy,iz))^2);
            end
        end
    end

    %electric field
    ex = zeros(nx,ny,nz);
    ey = zeros(nx,ny,nz);
    ez = zeros(nx,ny,nz);
    for ix = 2:nx-1
        for iy = 2:ny-1
            for iz = 2:nz-1
%                 ex(ix,iy,iz) = res(ix,iy,iz)*jx(ix,iy,iz) - (sy(ix,iy,iz)*bz(ix,iy,iz) - sz(ix,iy,iz)*by(ix,iy,iz))/rho(ix,iy,iz);
%                 ey(ix,iy,iz) = res(ix,iy,iz)*jy(ix,iy,iz) - (sz(ix,iy,iz)*bx(ix,iy,iz) - sx(ix,iy,iz)*bz(ix,iy,iz))/rho(ix,iy,iz);
%                 ez(ix,iy,iz) = res(ix,iy,iz)*jz(ix,iy,iz) - (sx(ix,iy,iz)*by(ix,iy,iz) - sy(ix,iy,iz)*bx(ix,iy,iz))/rho(ix,iy,iz);
                ex(ix,iy,iz) =  -(sy(ix,iy,iz)*bz(ix,iy,iz) - sz(ix,iy,iz)*by(ix,iy,iz))/rho(ix,iy,iz);
                ey(ix,iy,iz) =  -(sz(ix,iy,iz)*bx(ix,iy,iz) - sx(ix,iy,iz)*bz(ix,iy,iz))/rho(ix,iy,iz);
                ez(ix,iy,iz) =  -(sx(ix,iy,iz)*by(ix,iy,iz) - sy(ix,iy,iz)*bx(ix,iy,iz))/rho(ix,iy,iz);
            end
        end
    end
end


    
