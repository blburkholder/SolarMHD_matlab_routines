files = 10;
%s = 1-bottom 2-middle 3-top
s = 1;
%DO this if you need TO%%
%  for i = 0:files
%     flprepare2(i)
% %      prepare3(i)
% end

%choose 1 to do all files in run
%choose 0 to just compare first and last
doall = 0;

% choose 1 to use fortran fluid tracer
% choose 0 to use 2D matlab fluid tracer (only relevant for s = 1)
xuan = 0;

load('flt0.mat')
load('mhd0.mat')

%reolution
reso = 1;
%stay away from boundaries
off = 5;

[ssx,ssy] = size(fepo(off:reso:end-off,off:reso:end-off,1,1));
track = zeros(ssx,ssy,2,files+1);
z1 = fepo;

slice = fepo(:,:,s,:);
coor_x = slice(:,:,1,1);
coor_y = slice(:,:,1,2);
coor_z = slice(:,:,1,3);

pos1 = perio_stream3(x(2:end-1),y(2:end-1),z(2:end-1),bx(2:end-1,2:end-1,2:end-1),...
    by(2:end-1,2:end-1,2:end-1),bz(2:end-1,2:end-1,2:end-1),coor_x(off:reso:end-off,off:reso:end-off),...
    coor_y(off:reso:end-off,off:reso:end-off),...
    coor_z(off:reso:end-off,off:reso:end-off));

conj = zeros(length(pos1),2);
for i = 1:length(pos1)-1
        if ((pos1{i}(end,1) < x(1)) || (pos1{i}(end,1) > x(end)) ||...
            (pos1{i}(end,2) < y(1)) || (pos1{i}(end,2) > y(end)) ||...
            (pos1{i}(end,3) < z(1)) || (pos1{i}(end,3) > z(end)))
            conj(i,1) = nan;
        else
            %conj(i,1) = interp3(x,y,z,bz,pos{i}(end,1),pos{i}(end,2),pos{i}(end,3));
            conj(i,1) = pos1{i}(end,3);
        end
end

figure
cx = coor_x(off:reso:end-off,off:reso:end-off);
cy = coor_y(off:reso:end-off,off:reso:end-off);
h = scatter(cx(:),cy(:),[10],conj(:,1),'filled');
colorbar

if xuan
    for i = 1:files+1
        slice = fepo(:,:,s,:);

        coor_x = slice(:,:,1);
        coor_y = slice(:,:,2);

        track(:,:,1,i) = coor_x(off:reso:end-off,off:reso:end-off);
        track(:,:,2,i) = coor_y(off:reso:end-off,off:reso:end-off);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if doall
            figure

            coor_z = slice(off:reso:end-off,off:reso:end-off,3);
            pos2 = perio_stream3(x(2:end-1),y(2:end-1),z(2:end-1),...
                bx(2:end-1,2:end-1,2:end-1),by(2:end-1,2:end-1,2:end-1),...
                bz(2:end-1,2:end-1,2:end-1),coor_x(off:reso:end-off,off:reso:end-off),...
                coor_y(off:reso:end-off,off:reso:end-off),coor_z);

            for k = 1:length(pos2)-1
                ps = size(pos2{k});
                if (ps(1) + ps(2) < 3)
                    conj(k,2) = nan;
                elseif (pos2{k}(end,1) < x(2)) || (pos2{k}(end,1) > x(end-1)) ||...
                    (pos2{k}(end,2) < y(2)) || (pos2{k}(end,2) > y(end-1))
                    conj(k,2) = nan;
                else
                    conj(k,2) = pos2{k}(end,3);
                end
            end

            cx = coor_x(off:reso:end-off,off:reso:end-off);
            cy = coor_y(off:reso:end-off,off:reso:end-off);
            h = scatter(cx(:),cy(:),[10],conj(:,2),'filled');
            colorbar
            figure
            %cx = coor_x(off:reso:end-off,off:reso:end-off);
            %cy = coor_y(off:reso:end-off,off:reso:end-off);
            scatter(cx(:),cy(:),[10],(conj(:,2) - conj(:,1)),'filled');
        end        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if i ~= files+1
            if doall
                load(strcat('mhd',num2str(i),'.mat'))
            end
            load(['flt',num2str(i),'.mat'])
        end
    end
end

if ~doall
    load('mhd10.mat')
end

figure
if xuan
    if ~doall
        coor_z = slice(off:reso:end-off,off:reso:end-off,3);
        pos2 = perio_stream3(x(2:end-1),y(2:end-1),z(2:end-1),...
            bx(2:end-1,2:end-1,2:end-1),by(2:end-1,2:end-1,2:end-1),...
            bz(2:end-1,2:end-1,2:end-1),coor_x(off:reso:end-off,off:reso:end-off),...
            coor_y(off:reso:end-off,off:reso:end-off),coor_z);
    end
else
    [coor_x,coor_y] = velocity_perturbation(x,y,coor_x(off:reso:end-off,off:reso:end-off),...
            coor_y(off:reso:end-off,off:reso:end-off),time);
    pos2 = perio_stream3(x(2:end-1),y(2:end-1),z(2:end-1),...
        bx(2:end-1,2:end-1,2:end-1),by(2:end-1,2:end-1,2:end-1),...
        bz(2:end-1,2:end-1,2:end-1),coor_x,coor_y, 0.1*ones(size(coor_x)));
end

if ~doall
    for i = 1:length(pos2)-1
        %checking for field lines only 1 long
        ps = size(pos2{i});
        if (ps(1) + ps(2) < 3)
            conj(i,2) = nan;
        else
            conj(i,2) = pos2{i}(end,3);
        end
    end
end

if xuan
    if ~doall
        cx = coor_x(off:reso:end-off,off:reso:end-off);
        cy = coor_y(off:reso:end-off,off:reso:end-off);
        h = scatter(cx(:),cy(:),[10],conj(:,2),'filled');
        colorbar
    end
else
    pcolor(coor_x,coor_y,reshape(conj(:,2),[size(coor_x)]));
    shading interp
    colorbar
end

if ~doall
    figure
    title(gca,['z=',num2str(z1(1,1,s,3))]);
end

if xuan
    cx = coor_x(off:reso:end-off,off:reso:end-off);
    cy = coor_y(off:reso:end-off,off:reso:end-off);
    scatter(cx(:),cy(:),[10],(conj(:,2) - conj(:,1)),'filled');
    colorbar
else
    elg = reshape(conj(:,1),[size(coor_x)]);
    gle = reshape(conj(:,2),[size(coor_x)]);
    pcolor(coor_x, coor_y, gle-elg');
    shading interp
    colorbar
end

hold on
reso = 5;
if xuan
    for i = 1:files+1
        tx = track(1:reso:end,1:reso:end,1,i);
        ty = track(1:reso:end,1:reso:end,2,i);
        k = scatter(tx(:),ty(:),4,'k.');
        k.LineWidth = 0.1;
    end
    k = scatter(tx(:),ty(:),4,'rp');
    k.LineWidth = 0.1;
else
    slice = z1(:,:,s,:);
    coor_x = slice(:,:,1,1);
    coor_y = slice(:,:,1,2);
    velocity_perturbation(x,y,coor_x(1:7*reso:end,1:7*reso:end),coor_y(1:7*reso:end,1:7*reso:end),time);
end
