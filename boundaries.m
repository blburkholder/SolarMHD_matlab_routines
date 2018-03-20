%final argument determines what P is
% 1 parallel electric field
% 2 field aligned current
% 3 perpendicular current
%[jx,jy,jz,P] = get_j(nx,ny,nz,bx,by,bz,b0x,b0y,b0z,res,difx,dify,difz,0);

%load('mentrop5.mat')

%property which will paint bottom/top boundary
% prop = log(abs(ftekin)+eps);
% prop = bx(2:end-1,2:end-1,n);
prop = zfin1;

%quantity to project on field lines (cannot be negative)
fline_project = abs(sz)./rho;
%fline_project = perp_current;

%tag to decide whether to color flines
project_P_onto_fline = 0;

%boundary gradient min & max
grad_min = 1;
grad_max = 50;
%field line drawing resolution
reso = 41;
%prop quantity maximum in case of extraneous values
maxi = 150;
%how many times to split the boundaries
slines = 1;
%where to start flines
n = 2;

count = 500;
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
            bs(count,1) = xpo(j);
            bs(count,2) = ypo(i);
            %bs(count,1) = x(i);
            %bs(count,2) = y(j);
        end
    end
end
count

%use this when incorporating bottom & top boundaries
% load('mentrop5.mat')
% prop = bbaseft;

prop(prop > maxi) = maxi;
prop(prop < -maxi) = -maxi;
prop = prop/max(max(prop));

for j = 1:slines
    figure
    hold on

%     pcolor(ypo,xpo,prop)
%     pcolor(y(2:end-1),x(2:end-1),prop)
%     streamslice(y,x,by(:,:,n),bx(:,:,n))
%     shading interp
%     colormap(jet)
%     colorbar
%     daspect([1 1 1])
    scatter(bs((j-1)*floor(count/slines)+1:ceil(j*count/slines),1),...
        bs((j-1)*floor(count/slines)+1:ceil(j*count/slines),2),[2],'k.')
%     scatter3(bs((j-1)*floor(count/slines)+1:ceil(j*count/slines),2),...
%         bs((j-1)*floor(count/slines)+1:ceil(j*count/slines),1),...
%     z(end-1)*ones(length(bs((j-1)*floor(count/slines)+1:ceil(j*count/slines),1)),1),[2],'k.')

    plines = perio_stream3(x(2:end-1),y(2:end-1),z(2:end-1),...
        bx(2:end-1,2:end-1,2:end-1),by(2:end-1,2:end-1,2:end-1),...
        bz(2:end-1,2:end-1,2:end-1),bs((j-1)*floor(count/slines)+1:reso:ceil(j*count/slines),1),...
        bs((j-1)*floor(count/slines)+1:reso:ceil(j*count/slines),2),...
        z(n)*ones(length(bs((j-1)*floor(count/slines)+1:reso:ceil(j*count/slines),1)),1));

    if project_P_onto_fline == 0
        h = streamtube(plines,0.2);
        shading interp
        view(3)
    else
        b = 0;
        q = 0;
        for i = 1:length(plines)
            if ~isempty(plines{i})
                q = interp3(x(2:end-1),y(2:end-1),z(2:end-1),fline_project(2:end-1,2:end-1,2:end-1),plines{i}(:,1),plines{i}(:,2),plines{i}(:,3));
                if max(q) > b
                    b = max(q);
                end
                plines{i}(:,4) = q;
            end
        end
        maximum = b;
        for i = 1:length(plines)
            try
            patch([(plines{i}(:,1))' nan],[(plines{i}(:,2))' nan],...
                [(plines{i}(:,3))' nan],[(plines{i}(:,4)/maximum)' nan],...
                'EdgeColor','interp','FaceColor','none');
            catch

            end
        end
    end

    %fl without perio_stream3
%     streamline(stream3(y(2:end-1),x(2:end-1),z(2:end-1),...
%         by(2:end-1,2:end-1,2:end-1),bx(2:end-1,2:end-1,2:end-1),...
%         bz(2:end-1,2:end-1,2:end-1),bs((i-1)*floor(count/slines)+1:reso:ceil(i*count/slines),2),...
%         bs((i-1)*floor(count/slines)+1:reso:ceil(i*count/slines),1),...
%         z(n)*ones(length(bs((i-1)*floor(count/slines)+1:reso:ceil(i*count/slines),1)),1)));
%     flines = streamline(stream3(y(2:end-1),x(2:end-1),z(2:end-1),...
%         -by(2:end-1,2:end-1,2:end-1),-bx(2:end-1,2:end-1,2:end-1),...
%         -bz(2:end-1,2:end-1,2:end-1),bs((i-1)*floor(count/slines)+1:reso:ceil(i*count/slines),2),...
%         bs((i-1)*floor(count/slines)+1:reso:ceil(i*count/slines),1),...
%         z(n)*ones(length(bs((i-1)*floor(count/slines)+1:reso:ceil(i*count/slines),1)),1)));
%     set(flines,'Color','r');
%     view(3)

% %this will put red stars at the conjugate footpoints
%     conj = zeros(length(plines),3);
%     for l = 1:length(plines)
%         conj(l,1) = plines{l}(end,1);
%         conj(l,2) = plines{l}(end,2);
%         conj(l,3) = plines{l}(end,3);
%     end
%     hold on
%     %scatter3(conj(:,1),conj(:,2),conj(:,3),'rp')
%     scatter(conj(:,1),conj(:,2),'rp')

    daspect([1 1 1])
end




