%this one spends less time computing field lines
%but more time doing other stuff

flippers = gle-elg;
inds = find(abs(flippers)>1);
coor_x = coor_x1(off:reso:end-off,off:reso:end-off);
coor_y = coor_y1(off:reso:end-off,off:reso:end-off);
fcoor_x = coor_x(inds);
fcoor_y = coor_y(inds);

%for d = 1:files+1
for d = 5
    load(strcat('mhd',num2str(d-1),'.mat'));
    figure
    hold on
    
    [cor_x,cor_y] = velocity_perturbation(x,y,fcoor_x,fcoor_y,time);
    plines = perio_stream3(x(2:end-1),y(2:end-1),z(2:end-1),...
        bx(2:end-1,2:end-1,2:end-1),by(2:end-1,2:end-1,2:end-1),...
        bz(2:end-1,2:end-1,2:end-1),cor_x,cor_y,0*ones(size(cor_y)));

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

    fac = zeros(nx,ny,nz);
    for ix = 2:nx-1
        for iy = 2:ny-1
            for iz = 2:nz-1
                BBB = sqrt(bx(ix,iy,iz)^2 + by(ix,iy,iz)^2 + bz(ix,iy,iz)^2);
                dir_bx = bx(ix,iy,iz)/BBB;
                dir_by = by(ix,iy,iz)/BBB;
                dir_bz = bz(ix,iy,iz)/BBB;
                fac(ix,iy,iz) = jx(ix,iy,iz)*dir_bx + jy(ix,iy,iz)*dir_by + jz(ix,iy,iz)*dir_bz;
            end
        end
    end

    for i = 1:length(plines)
        plines{i}(:,4) = interp3(x(2:end-1),y(2:end-1),z(2:end-1),fac(2:end-1,2:end-1,2:end-1),plines{i}(:,1),plines{i}(:,2),plines{i}(:,3));
    end

    maximum = max(max(max(fac)));
    for i = 1:length(plines)
        if plines{i}(end,3) > 5
            scatter3(plines{i}(end,1),plines{i}(end,2),plines{i}(end,3),'bp')
        else
            scatter3(plines{i}(end,1),plines{i}(end,2),plines{i}(end,3),'rp')
        end
        try

        patch([(plines{i}(:,1))' nan],[(plines{i}(:,2))' nan],...
            [(plines{i}(:,3))' nan],[abs((plines{i}(:,4)/maximum)') nan],...
            'EdgeColor','interp','FaceColor','none');
        catch
        
        end
    end

    axis([0 93 0 93 0 65])
    daspect([1 1 1])
    view(3)
    colorbar
end

    