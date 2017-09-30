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


% fac = zeros(nx,ny,nz);
% for ix = 2:nx-1
%     for iy = 2:ny-1
%         for iz = 2:nz-1
%             BBB = sqrt(bx(ix,iy,iz)^2 + by(ix,iy,iz)^2 + bz(ix,iy,iz)^2);
%             dir_bx = bx(ix,iy,iz)/BBB;
%             dir_by = by(ix,iy,iz)/BBB;
%             dir_bz = bz(ix,iy,iz)/BBB;
%             fac(ix,iy,iz) = jx(ix,iy,iz)*dir_bx + jy(ix,iy,iz)*dir_by + jz(ix,iy,iz)*dir_bz;
%         end
%     end
% end

% perp_current = zeros(nx,ny,nz);
% for ix = 2:nx-1
%     for iy=2:ny-1
%         for iz = 2:nz-1
%             perp_current(ix,iy,iz) = sqrt((jy(ix,iy,iz)*bz(ix,iy,iz) -...
%                 jz(ix,iy,iz)*by(ix,iy,iz))^2 + (jz(ix,iy,iz)*bx(ix,iy,iz) -...
%                 jx(ix,iy,iz)*bz(ix,iy,iz))^2 + (jx(ix,iy,iz)*by(ix,iy,iz) -...
%                 jy(ix,iy,iz)*bx(ix,iy,iz))^2)/sqrt((bx(ix,iy,iz))^2 +...
%                 (by(ix,iy,iz))^2 + (bz(ix,iy,iz))^2);
%         end
%     end
% end

%load('mentrop5.mat')
%prop = log(abs(ftekin)+eps);
n = nz-1;
%bz_n = bz./(sqrt(bx.^2+by.^2+bz.^2));
%theta = acos(bz_n(2:end-1,2:end-1,n));

%prop = bx(2:end-1,2:end-1,n);
prop = ftvol;
fline_project = abs(sz)./rho;
%fline_project = perp_current;

%boundary gradient min & max
grad_min = 20;
grad_max = 50;
%field line drawing resolution
reso = 50;
%flux tube integrator quantity maximum
maxi = 150;
%how many times to split the boundaries
slines = 5;

[row,col] = size(prop);
bs = zeros(count,2);
count = 0;
for i = 3:row-2
    for j = 3:col-2
        dx_p = sqrt((prop(i,j) - prop(i+1,j))^2);
        dx_m = sqrt((prop(i,j) - prop(i-1,j))^2);
        dy_p = sqrt((prop(i,j) - prop(i,j+1))^2);
        dy_m = sqrt((prop(i,j) - prop(i,j-1))^2);
        dd_pp = sqrt((prop(i,j) - prop(i+1,j+1))^2);
        dd_pn = sqrt((prop(i,j) - prop(i+1,j-1))^2);
        dd_np = sqrt((prop(i,j) - prop(i-1,j+1))^2);
        dd_nn = sqrt((prop(i,j) - prop(i-1,j-1))^2);
        if dx_p > grad_min || dx_m > grad_min || dy_p > grad_min || dy_m > grad_min ||...
            dd_pp > grad_min || dd_pn > grad_min || dd_np > grad_min || dd_nn > grad_min
%         if (dx_p > grad_min && dx_p < grad_max) && (dx_m > grad_min && dx_m < grad_max)...
%             && (dy_p > grad_min && dy_p < grad_max) && (dy_m > grad_min && dy_m < grad_max)...
%             && (dd_pp > grad_min && dd_pp < grad_max) && (dd_pn > grad_min && dd_pn < grad_max)...
%             && (dd_np > grad_min && dd_np < grad_max) && (dd_nn > grad_min && dd_nn < grad_max)
            count = count + 1;
            bs(count,1) = xpo(i);
            bs(count,2) = ypo(j);
            %bs(count,1) = x(i);
            %bs(count,2) = y(j);
        end
    end
end
count

% load('mentrop5.mat')
% prop = bbaseft;

prop(prop > maxi) = maxi;
prop(prop < -maxi) = -maxi;
prop = prop/max(max(prop));

for i = 1:slines
    figure
    hold on
    pcolor(ypo,xpo,prop)
    %pcolor(y(2:end-1),x(2:end-1),prop)
    %streamslice(y,x,by(:,:,n),bx(:,:,n))
    shading interp
     scatter(bs((i-1)*floor(count/slines)+1:ceil(i*count/slines),2),...
        bs((i-1)*floor(count/slines)+1:ceil(i*count/slines),1),[2],'w.')
%      scatter3(bs((i-1)*floor(count/slines)+1:ceil(i*count/slines),2),...
%          bs((i-1)*floor(count/slines)+1:ceil(i*count/slines),1),...
%         z(end-1)*ones(length(bs((i-1)*floor(count/slines)+1:ceil(i*count/slines),1)),1),[2],'w.')
    colorbar

    24
    %tubey tubey from bottom
    % flines = streamtube(perio_stream3(y(2:end-1),x(2:end-1),z(2:end-1),...
    %     -by(2:end-1,2:end-1,2:end-1),-bx(2:end-1,2:end-1,2:end-1),...
    %     -bz(2:end-1,2:end-1,2:end-1),bs(1:reso:end,2),bs(1:reso:end,1),zeros(length(bs(1:reso:end,1)),1)),0.4);
    % shading interp
    % view(3)

    %fl from bottom
    % flines = streamline(perio_stream3(y(2:end-1),x(2:end-1),z(2:end-1),...
    %     -by(2:end-1,2:end-1,2:end-1),-bx(2:end-1,2:end-1,2:end-1),...
    %     -bz(2:end-1,2:end-1,2:end-1),bs(1:reso:end,2),bs(1:reso:end,1),zeros(length(bs(1:reso:end,1)),1)));
    % set(flines,'Color','r');
    % view(3)
    %top
    flines = streamline(perio_stream3(y(2:end-1),x(2:end-1),z(2:end-1),...
        -by(2:end-1,2:end-1,2:end-1),-bx(2:end-1,2:end-1,2:end-1),...
        -bz(2:end-1,2:end-1,2:end-1),bs((i-1)*floor(count/slines)+1:reso:ceil(i*count/slines),2),...
        bs((i-1)*floor(count/slines)+1:reso:ceil(i*count/slines),1),...
        z(2)*ones(length(bs((i-1)*floor(count/slines)+1:reso:ceil(i*count/slines),1)),1)));
    set(flines,'Color','r');
    view(3)

    %fl from anywhere without perio_stream3
%     streamline(stream3(y(2:end-1),x(2:end-1),z(2:end-1),...
%         by(2:end-1,2:end-1,2:end-1),bx(2:end-1,2:end-1,2:end-1),...
%         bz(2:end-1,2:end-1,2:end-1),bs(1:reso:end,2),bs(1:reso:end,1),z(n)*ones(length(bs(1:reso:end,1)),1),));
%     flines = streamline(stream3(y(2:end-1),x(2:end-1),z(2:end-1),...
%         -by(2:end-1,2:end-1,2:end-1),-bx(2:end-1,2:end-1,2:end-1),...
%         -bz(2:end-1,2:end-1,2:end-1),bs(1:reso:end,2),bs(1:reso:end,1),z(n)*ones(length(bs(1:reso:end,1)),1)));
%     set(flines,'Color','r');
    % view(3)

    %fl colored wrt whatever the f you want
    %bottom
%     plines = perio_stream3(y(2:end-1),x(2:end-1),z(2:end-1),...
%          by(2:end-1,2:end-1,2:end-1),bx(2:end-1,2:end-1,2:end-1),...
%          bz(2:end-1,2:end-1,2:end-1),bs((i-1)*floor(count/slines)+1:...
%             reso:ceil(i*count/slines),2),bs((i-1)*floor(count/slines)+1:...
%             reso:ceil(i*count/slines),1),zeros(length(bs((i-1)*...
%             floor(count/slines)+1:reso:ceil(i*count/slines),1)),1));
    %top
%     plines = perio_stream3(y(2:end-1),x(2:end-1),z(2:end-1),...
%          by(2:end-1,2:end-1,2:end-1),bx(2:end-1,2:end-1,2:end-1),...
%          bz(2:end-1,2:end-1,2:end-1),bs((i-1)*floor(count/slines)+1:...
%             reso:ceil(i*count/slines),2),bs((i-1)*floor(count/slines)+1:...
%             reso:ceil(i*count/slines),1),z(end-1)*ones(length(bs((i-1)*...
%             floor(count/slines)+1:reso:ceil(i*count/slines),1)),1));

%     conj = zeros(length(plines),3);
%     for l = 1:length(plines)
%         conj(l,1) = plines{l}(end,1);
%         conj(l,2) = plines{l}(end,2);
%         conj(l,3) = plines{l}(end,3);
%     end
%     scatter3(conj(:,1),conj(:,2),conj(:,3),'rp')


%     b = 0;
%     q = 0;
%     for i = 1:length(plines)
%         if ~isempty(plines{i})
%             q = interp3(y(2:end-1),x(2:end-1),z(2:end-1),fline_project(2:end-1,2:end-1,2:end-1),plines{i}(:,1),plines{i}(:,2),plines{i}(:,3));
%             if max(q) > b
%                 b = max(q);
%             end
%             plines{i}(:,4) = q;
%         end
%     end
%     %dont try to use a quantity that can go negative
%     maximum = b;
%     for i = 1:length(plines)
%             try
%             patch([(plines{i}(:,1))' nan],[(plines{i}(:,2))' nan],...
%                 [(plines{i}(:,3))' nan],[(plines{i}(:,4)/maximum)' nan],...
%                 'EdgeColor','interp','FaceColor','none');
%             catch
% 
%             end
%     end

    daspect([1 1 1])
end




