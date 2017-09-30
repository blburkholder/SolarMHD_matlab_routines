%for i = 0:10
%load(strcat('mhd',num2str(i),'.mat'))

[sbx,sby,sbz] = size(bx);
nulls = [];
ex = 3;
for k = 1:8 %cut the volume into 8 subvolumes for speed
    if k == 1
        b1x = bx(2:floor(sbx/2)+ex,2:floor(sby/2)+ex,2:floor(sbz/2)+ex);
        b1y = by(2:floor(sbx/2)+ex,2:floor(sby/2)+ex,2:floor(sbz/2)+ex);
        b1z = bz(2:floor(sbx/2)+ex,2:floor(sby/2)+ex,2:floor(sbz/2)+ex);
        x1 = x(2:floor(length(x)/2)+ex);
        y1 = y(2:floor(length(y)/2)+ex);
        z1 = z(2:floor(length(z)/2)+ex);
    elseif k == 2
        b1x = bx(floor(sbx/2)-ex:end-1,floor(sby/2)-ex:end-1,floor(sbz/2)-ex:end-1);
        b1y = by(floor(sbx/2)-ex:end-1,floor(sby/2)-ex:end-1,floor(sbz/2)-ex:end-1);
        b1z = bz(floor(sbx/2)-ex:end-1,floor(sby/2)-ex:end-1,floor(sbz/2)-ex:end-1);
        x1 = x(floor(length(x)/2)-ex:end-1);
        y1 = y(floor(length(y)/2)-ex:end-1);
        z1 = z(floor(length(z)/2)-ex:end-1);
    elseif k == 3
        b1x = bx(floor(sbx/2)-ex:end-1,floor(sby/2)-ex:end-1,2:floor(sbz/2)+ex);
        b1y = by(floor(sbx/2)-ex:end-1,floor(sby/2)-ex:end-1,2:floor(sbz/2)+ex);
        b1z = bz(floor(sbx/2)-ex:end-1,floor(sby/2)-ex:end-1,2:floor(sbz/2)+ex);
        x1 = x(floor(length(x)/2)-ex:end-1);
        y1 = y(floor(length(y)/2)-ex:end-1);
        z1 = z(2:floor(length(z)/2)+ex);
    elseif k == 4
        b1x = bx(floor(sbx/2)-ex:end-1,2:floor(sby/2)+ex,floor(sbz/2)-ex:end-1);
        b1y = by(floor(sbx/2)-ex:end-1,2:floor(sby/2)+ex,floor(sbz/2)-ex:end-1);
        b1z = bz(floor(sbx/2)-ex:end-1,2:floor(sby/2)+ex,floor(sbz/2)-ex:end-1);
        x1 = x(floor(length(x)/2)-ex:end-1);
        y1 = y(2:floor(length(y)/2)+ex);
        z1 = z(floor(length(z)/2)-ex:end-1);
    elseif k == 5
        b1x = bx(2:floor(sbx/2)+ex,floor(sby/2)-ex:end-1,floor(sbz/2)-ex:end-1);
        b1y = by(2:floor(sbx/2)+ex,floor(sby/2)-ex:end-1,floor(sbz/2)-ex:end-1);
        b1z = bz(2:floor(sbx/2)+ex,floor(sby/2)-ex:end-1,floor(sbz/2)-ex:end-1);
        x1 = x(2:floor(length(x)/2)+ex);
        y1 = y(floor(length(y)/2)-ex:end-1);
        z1 = z(floor(length(z)/2)-ex:end-1);
    elseif k == 6
        b1x = bx(2:floor(sbx/2)+ex,2:floor(sby/2)+ex,floor(sbz/2)-ex:end-1);
        b1y = by(2:floor(sbx/2)+ex,2:floor(sby/2)+ex,floor(sbz/2)-ex:end-1);
        b1z = bz(2:floor(sbx/2)+ex,2:floor(sby/2)+ex,floor(sbz/2)-ex:end-1);
        x1 = x(2:floor(length(x)/2)+ex);
        y1 = y(2:floor(length(y)/2)+ex);
        z1 = z(floor(length(z)/2)-ex:end-1);
    elseif k == 7
        b1x = bx(2:floor(sbx/2)+ex,floor(sby/2)-ex:end-1,2:floor(sbz/2)+ex);
        b1y = by(2:floor(sbx/2)+ex,floor(sby/2)-ex:end-1,2:floor(sbz/2)+ex);
        b1z = bz(2:floor(sbx/2)+ex,floor(sby/2)-ex:end-1,2:floor(sbz/2)+ex);
        x1 = x(2:floor(length(x)/2)+ex);
        y1 = y(floor(length(y)/2)-ex:end-1);
        z1 = z(2:floor(length(z)/2)+ex);
    else
        b1x = bx(floor(sbx/2)-ex:end-1,2:floor(sby/2)+ex,2:floor(sbz/2)+ex);
        b1y = by(floor(sbx/2)-ex:end-1,2:floor(sby/2)+ex,2:floor(sbz/2)+ex);
        b1z = bz(floor(sbx/2)-ex:end-1,2:floor(sby/2)+ex,2:floor(sbz/2)+ex);
        x1 = x(floor(length(x)/2)-ex:end-1);
        y1 = y(2:floor(length(y)/2)+ex);
        z1 = z(2:floor(length(z)/2)+ex);   
    end
    for g = 1:8
        if g == 1
            blx = b1x(1:floor(length(x1)/2)+ex,1:floor(length(y1)/2)+ex,1:floor(length(z1)/2)+ex);
            bly = b1y(1:floor(length(x1)/2)+ex,1:floor(length(y1)/2)+ex,1:floor(length(z1)/2)+ex);
            blz = b1z(1:floor(length(x1)/2)+ex,1:floor(length(y1)/2)+ex,1:floor(length(z1)/2)+ex);
            xx = x1(1:floor(length(x1)/2)+ex);
            yy = y1(1:floor(length(y1)/2)+ex);
            zz = z1(1:floor(length(z1)/2)+ex);   
        elseif g == 2
            blx = b1x(1:floor(length(x1)/2)+ex,floor(length(y1)/2)-ex:end,1:floor(length(z1)/2)+ex);
            bly = b1y(1:floor(length(x1)/2)+ex,floor(length(y1)/2)-ex:end,1:floor(length(z1)/2)+ex);
            blz = b1z(1:floor(length(x1)/2)+ex,floor(length(y1)/2)-ex:end,1:floor(length(z1)/2)+ex);
            xx = x1(1:floor(length(x1)/2)+ex);
            yy = y1(floor(length(y1)/2)-ex:end);
            zz = z1(1:floor(length(z1)/2)+ex);  
        elseif g == 3
            blx = b1x(1:floor(length(x1)/2)+ex,1:floor(length(y1)/2)+ex,floor(length(z1)/2)-ex:end);
            bly = b1y(1:floor(length(x1)/2)+ex,1:floor(length(y1)/2)+ex,floor(length(z1)/2)-ex:end);
            blz = b1z(1:floor(length(x1)/2)+ex,1:floor(length(y1)/2)+ex,floor(length(z1)/2)-ex:end);
            xx = x1(1:floor(length(x1)/2)+ex);
            yy = y1(1:floor(length(y1)/2)+ex);
            zz = z1(floor(length(z1)/2)-ex:end);  
        elseif g == 4
            blx = b1x(1:floor(length(x1)/2)+ex,floor(length(y1)/2)-ex:end,floor(length(z1)/2)-ex:end);
            bly = b1y(1:floor(length(x1)/2)+ex,floor(length(y1)/2)-ex:end,floor(length(z1)/2)-ex:end);
            blz = b1z(1:floor(length(x1)/2)+ex,floor(length(y1)/2)-ex:end,floor(length(z1)/2)-ex:end);
            xx = x1(1:floor(length(x1)/2)+ex);
            yy = y1(floor(length(y1)/2)-ex:end);
            zz = z1(floor(length(z1)/2)-ex:end);  
        elseif g == 5
            blx = b1x(floor(length(x1)/2)-ex:end,1:floor(length(y1)/2)+ex,1:floor(length(z1)/2)+ex);
            bly = b1y(floor(length(x1)/2)-ex:end,1:floor(length(y1)/2)+ex,1:floor(length(z1)/2)+ex);
            blz = b1z(floor(length(x1)/2)-ex:end,1:floor(length(y1)/2)+ex,1:floor(length(z1)/2)+ex);
            xx = x1(floor(length(x1)/2)-ex:end);
            yy = y1(1:floor(length(y1)/2)+ex);
            zz = z1(1:floor(length(z1)/2)+ex);   
        elseif g == 6
            blx = b1x(floor(length(x1)/2)-ex:end,floor(length(y1)/2)-ex:end,1:floor(length(z1)/2)+ex);
            bly = b1y(floor(length(x1)/2)-ex:end,floor(length(y1)/2)-ex:end,1:floor(length(z1)/2)+ex);
            blz = b1z(floor(length(x1)/2)-ex:end,floor(length(y1)/2)-ex:end,1:floor(length(z1)/2)+ex);
            xx = x1(floor(length(x1)/2)-ex:end);
            yy = y1(floor(length(y1)/2)-ex:end);
            zz = z1(1:floor(length(z1)/2)+ex);  
        elseif g == 7
            blx = b1x(floor(length(x1)/2)-ex:end,1:floor(length(y1)/2)+ex,floor(length(z1)/2)-ex:end);
            bly = b1y(floor(length(x1)/2)-ex:end,1:floor(length(y1)/2)+ex,floor(length(z1)/2)-ex:end);
            blz = b1z(floor(length(x1)/2)-ex:end,1:floor(length(y1)/2)+ex,floor(length(z1)/2)-ex:end);
            xx = x1(floor(length(x1)/2)-ex:end);
            yy = y1(1:floor(length(y1)/2)+ex);
            zz = z1(floor(length(z1)/2)-ex:end);  
        elseif g == 8
            blx = b1x(floor(length(x1)/2)-ex:end,floor(length(y1)/2)-ex:end,floor(length(z1)/2)-ex:end);
            bly = b1y(floor(length(x1)/2)-ex:end,floor(length(y1)/2)-ex:end,floor(length(z1)/2)-ex:end);
            blz = b1z(floor(length(x1)/2)-ex:end,floor(length(y1)/2)-ex:end,floor(length(z1)/2)-ex:end);
            xx = x1(floor(length(x1)/2)-ex:end);
            yy = y1(floor(length(y1)/2)-ex:end);
            zz = z1(floor(length(z1)/2)-ex:end); 
        end
        %figure
        %using 50% of available triangles
        p1 = patch(isosurface(yy,xx,zz,bly,0),'visible','off');
        p1 = reducepatch(p1,0.8);
        %p1.EdgeColor = 'red';
        p2 = patch(isosurface(yy,xx,zz,blx,0),'visible','off');
        p2 = reducepatch(p2,0.8);
        %p2.EdgeColor = 'blue';
        p3 = patch(isosurface(yy,xx,zz,blz,0),'visible','off');
        p3 = reducepatch(p3,0.8);
        %p3.EdgeColor = 'green';

    %     clf
    %     S=s1; trisurf(S.faces, S.vertices(:,1),S.vertices(:,2),S.vertices(:,3),'FaceAlpha', 0.5, 'FaceColor', 'r');
    %     S=s2; trisurf(S.faces, S.vertices(:,1),S.vertices(:,2),S.vertices(:,3),'FaceAlpha', 0.5, 'FaceColor', 'g');
    %     S=s3; trisurf(S.faces, S.vertices(:,1),S.vertices(:,2),S.vertices(:,3),'FaceAlpha', 0.5, 'FaceColor', 'b');
        'patches reduced, computing surface intersections'
        try
            if ~isempty(p1.vertices) || ~isempty(p2.vertices)
                [~, Surf12] = SurfaceIntersection(p1, p2);
                sv1 = Surf12.vertices;
            else
                sv1 = [];        
            end
        catch
            'S12 failed'
            sv1 = [];
        end
        try
            'intersection 1 done'
            if ~isempty(p1.vertices) || ~isempty(p3.vertices)
                [~, Surf13] = SurfaceIntersection(p1, p3);
                sv2 = Surf13.vertices;
            else
                sv2 = [];
            end
        catch
            'S13 failed'
            sv2 = [];
        end
        try
            'intersection 2 done'
            if ~isempty(p3.vertices) || ~isempty(p2.vertices)
                [~, Surf23] = SurfaceIntersection(p3, p2);
                sv3 = Surf23.vertices;
            else
                sv3 = [];
            end
            'intersection 3 done'
        catch
            'S23 failed'
            sv3 = [];
        end

    %     clf
    %     hold on
    %     if ~isempty(Surf12.vertices)
    %         trisurf(Surf12.faces, Surf12.vertices(:,1),Surf12.vertices(:,2),Surf12.vertices(:,3),'EdgeColor', 'r', 'FaceColor', 'r');
    %     end
    %     if ~isempty(Surf13.vertices)
    %         trisurf(Surf13.faces, Surf13.vertices(:,1),Surf13.vertices(:,2),Surf13.vertices(:,3),'EdgeColor', 'g', 'FaceColor', 'g');
    %     end
    %     if ~isempty(Surf23.vertices)
    %         trisurf(Surf23.faces, Surf23.vertices(:,1),Surf23.vertices(:,2),Surf23.vertices(:,3),'EdgeColor', 'b', 'FaceColor', 'b');
    %     end
    %     title ('Surface/Surface intersections')
    %     legend({'x-y', 'y-z', 'x-z'});

        if ~isempty(sv1)
            o1 = 0;
            c1 = interp3(yy,xx,zz,blz,sv1(1,1),sv1(1,2),sv1(1,3));
            if c1 > 0
                o1 = 1; 
            else
                o1 = -1;
            end
            for i = 2:length(sv1)
                c2 = interp3(yy,xx,zz,blz,sv1(i,1),sv1(i,2),sv1(i,3));
                if c2 > 0
                    o2 = 1;
                elseif c2 < 0
                    o2 = -1;
                else
                    o2 = 0;
                end
                if o1 + o2 == 0
                    dist = (sv1(i,1)-sv1(i-1,1))^2 + (sv1(i,2)-sv1(i-1,2))^2 + (sv1(i,3)-sv1(i-1,3))^2;
                    if sqrt(dist) < 0.5
                        bz_dif = c2-c1;
                        sl = abs(c1/bz_dif);
                        v = sv1(i,:) - sv1(i-1,:);
                        r = sv1(i-1,:) + sl*v;
                        bbbz = interp3(yy,xx,zz,blz,r(1),r(2),r(3));
                        %scatter3(r(1),r(2),r(3))
                        if abs(bbbz) < 0.1
                            nulls(end+1,:) = r;
                        end
                    end
                end
                c1 = c2;   
                o1 = o2;
            end
        end

        if ~isempty(sv2)
            o1 = 0;
            c1 = interp3(yy,xx,zz,blx,sv2(1,1),sv2(1,2),sv2(1,3));
            if c1 > 0
                o1 = 1; 
            else
                o1 = -1;
            end
            for i = 2:length(sv2)
                c2 = interp3(yy,xx,zz,blx,sv2(i,1),sv2(i,2),sv2(i,3));
                if c2 > 0
                    o2 = 1;
                elseif c2 < 0
                    o2 = -1;
                else
                    o2 = 0;
                end
                if o1 + o2 == 0
                    dist = (sv2(i,1)-sv2(i-1,1))^2 + (sv2(i,2)-sv2(i-1,2))^2 + (sv2(i,3)-sv2(i-1,3))^2;
                    if sqrt(dist) < 0.5
                        bx_dif = c2-c1;
                        sl = abs(c1/bx_dif);
                        v = sv2(i,:) - sv2(i-1,:);
                        r = sv2(i-1,:) + sl*v;
                        bbbx = interp3(yy,xx,zz,blx,r(1),r(2),r(3));
                        if abs(bbbx) < 0.1
                            nulls(end+1,:) = r;
                        end
                    end
                end
                c1 = c2;   
                o1 = o2;
            end
        end

        if ~isempty(sv3)
            o1 = 0;
            c1 = interp3(yy,xx,zz,bly,sv3(1,1),sv3(1,2),sv3(1,3));
            if c1 > 0
                o1 = 1; 
            else
                o1 = -1;
            end
            for i = 2:length(sv3)
                c2 = interp3(yy,xx,zz,bly,sv3(i,1),sv3(i,2),sv3(i,3));
                if c2 > 0
                    o2 = 1;
                elseif c2 < 0
                    o2 = -1;
                else
                    o2 = 0;
                end
                if o1 + o2 == 0
                    dist = (sv3(i,1)-sv3(i-1,1))^2 + (sv3(i,2)-sv3(i-1,2))^2 + (sv3(i,3)-sv3(i-1,3))^2;
                    if sqrt(dist) < 0.5
                        by_dif = c2-c1;
                        sl = abs(c1/by_dif);
                        v = sv3(i,:) - sv3(i-1,:);
                        r = sv3(i-1,:) + sl*v;
                        bbby = interp3(yy,xx,zz,bly,r(1),r(2),r(3));
                        if abs(bbby) < 0.1
                            nulls(end+1,:) = r;
                        end
                    end
                end
                c1 = c2;   
                o1 = o2;
            end
        end
    end
end
