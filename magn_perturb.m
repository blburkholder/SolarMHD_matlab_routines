function [dbx,dby,dbz,poynt] = magn_perturb(bx,by,bz,rho,xpo,ypo,zp1,xfin1,yfin1,zfin1,x,y,z,sx,sy,sz,ex,ey,ez)
    b = sqrt((bx.^2+by.^2+bz.^2));
    Va = b./sqrt((rho));

    xxfin = zeros(size(xfin1));
    yyfin = zeros(size(yfin1));

    n = 2;
    for i = 1:length(xfin1)
        for j = 1:length(xfin1)
            if xfin1(i,j) > x(end-1)
                xxfin(i,j) = x(end-1) - mod(xfin1(i,j),x(end-1));
            elseif xfin1(i,j) < x(end-1)
                xxfin(i,j) = abs(xfin1(i,j));
            end

            if yfin1(i,j) > y(end-1)
                yyfin(i,j) = y(end-1) - mod(yfin1(i,j),y(end-1));
            elseif yfin1(i,j) < y(end-1)
                yyfin(i,j) = abs(yfin1(i,j));
            end
        end
    end

    % bbz = interp3(y,x,z,bz,mod(yfin1,y(end-1)),mod(xfin1,x(end-1)),z(n)*ones(size(zfin1)));
    % bby = interp3(y,x,z,by,mod(yfin1,y(end-1)),mod(xfin1,x(end-1)),z(n)*ones(size(zfin1)));
    % bbx = interp3(y,x,z,bx,mod(yfin1,y(end-1)),mod(xfin1,x(end-1)),z(n)*ones(size(zfin1)));
    % VVa = interp3(y,x,z,Va,mod(yfin1,y(end-1)),mod(xfin1,x(end-1)),z(n)*ones(size(zfin1)));
    % vx = interp3(y,x,z,sx./rho,mod(yfin1,y(end-1)),mod(xfin1,x(end-1)),z(n)*ones(size(zfin1)));
    % vy = interp3(y,x,z,sy./rho,mod(yfin1,y(end-1)),mod(xfin1,x(end-1)),z(n)*ones(size(zfin1)));

    %IM PRETTY SURE THIS IS THE RIGHT WAY TO MOD THE COORDINATES
    bbz = interp3(y,x,z,bz,xxfin,yyfin,z(n)*ones(size(zfin1)));
    bby = interp3(y,x,z,by,xxfin,yyfin,z(n)*ones(size(zfin1)));
    bbx = interp3(y,x,z,bx,xxfin,yyfin,z(n)*ones(size(zfin1)));

    eex = interp3(y,x,z,ex,xxfin,yyfin,z(n)*ones(size(zfin1)));
    eey = interp3(y,x,z,ey,xxfin,yyfin,z(n)*ones(size(zfin1)));
    eez = interp3(y,x,z,ez,xxfin,yyfin,z(n)*ones(size(zfin1)));

    VVa = interp3(y,x,z,Va,xxfin,yyfin,z(n)*ones(size(zfin1)));
    vx = interp3(y,x,z,sx./rho,xxfin,yyfin,z(n)*ones(size(zfin1)));
    vy = interp3(y,x,z,sy./rho,xxfin,yyfin,z(n)*ones(size(zfin1)));
    vz = interp3(y,x,z,sz./rho,xxfin,yyfin,z(n)*ones(size(zfin1)));

    dbx = vx.*(bbx.^2 + bby.^2 + bbz.^2)./VVa;
    dby = vy.*(bbx.^2 + bby.^2 + bbz.^2)./VVa;
    dbz = vz.*(bbx.^2 + bby.^2 + bbz.^2)./VVa;

    poynt = ((eey.*dbz - eez.*dby).*bbx + (eez.*dbx - eex.*dbz).*bby + (eex.*dby - eey.*dbx).*bbz)./sqrt(bbx.^2+bby.^2+bbz.^2);
end
