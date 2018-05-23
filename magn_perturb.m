function [dbx,dby,dbz,poynt] = magn_perturb(bx,by,bz,rho,xpo,ypo,xfin1,yfin1,zfin1,x,y,z,sx,sy,sz,ex,ey,ez)
    Vax = abs(bx)./sqrt(rho);
    Vay = abs(by)./sqrt(rho);
    Vaz = abs(bz)./sqrt(rho);

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

    %IM PRETTY SURE THIS IS THE RIGHT WAY TO MOD THE COORDINATES
    bbz = interp3(x,y,z,bz,xxfin,yyfin,zfin1);
    bby = interp3(x,y,z,by,xxfin,yyfin,zfin1);
    bbx = interp3(x,y,z,bx,xxfin,yyfin,zfin1);
    bb = sqrt(bbx.^2+bby.^2+bbz.^2);

    eex = interp3(x,y,z,ex,xxfin,yyfin,zfin1);
    eey = interp3(x,y,z,ey,xxfin,yyfin,zfin1);
    eez = interp3(x,y,z,ez,xxfin,yyfin,zfin1);

    VVax = interp3(x,y,z,Vax,xxfin,yyfin,zfin1);
    VVay = interp3(x,y,z,Vay,xxfin,yyfin,zfin1);
    VVaz = interp3(x,y,z,Vaz,xxfin,yyfin,zfin1);
    VVa = sqrt(VVax.^2 + VVay.^2 + VVaz.^2);

    vx = interp3(x,y,z,sx./rho,xxfin,yyfin,zfin1);
    vy = interp3(x,y,z,sy./rho,xxfin,yyfin,zfin1);
    vz = interp3(x,y,z,sz./rho,xxfin,yyfin,zfin1); 
    
    dbx = -vx.*bb./VVa;
    dby = -vy.*bb./VVa;
    dbz = -vz.*bb./VVa;

    poynt = ((eey.*dbz - eez.*dby).*bbx + (eez.*dbx - eex.*dbz).*bby + (eex.*dby - eey.*dbx).*bbz)./sqrt(bbx.^2+bby.^2+bbz.^2);
end
