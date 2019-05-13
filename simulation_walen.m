function [inds1,inds2] = simulation_walen(x,bx,by,bz,vx,vy,vz,vax,vay,vaz,wll,wlh)  
    %wlh = walen limit high
    %wll = walen limit low

    max_windl = 1;
    vxd = zeros(2,1);
    vyd = zeros(2,1);
    vzd = zeros(2,1);
    vaxd = zeros(2,1);
    vayd = zeros(2,1);
    vazd = zeros(2,1);

    inds1 = [];
    inds2 = [];

    start = 2;
    first = 0;
    here = 0;

    ppp = zeros(length(vx),1);

    while start <= length(vx) - 2
        vxd(1) = abs(vx(start) - vx(start+1)); vxd(2) = abs(vx(start) - vx(start-1));
        vyd(1) = abs(vy(start) - vy(start+1)); vyd(2) = abs(vy(start) - vy(start-1));  
        vzd(1) = abs(vz(start) - vz(start+1)); vzd(2) = abs(vz(start) - vz(start-1));   
        vaxd(1) = abs(vax(start) - vax(start+1)); vaxd(2) = abs(vax(start) - vax(start-1));
        vayd(1) = abs(vay(start) - vay(start+1)); vayd(2) = abs(vay(start) - vay(start-1)); 
        vazd(1) = abs(vaz(start) - vaz(start+1)); vazd(2) = abs(vaz(start) - vaz(start-1));

        p = [vxd;vyd;vzd]\[vaxd;vayd;vazd];
        mslope = p(1)*[vxd;vyd;vzd];
        [goodness ,~] = rsquare([vaxd;vayd;vazd],mslope);

%        ppp(start) = p(1);

        if (p(1) >= wll) && (p(1) <= wlh) && goodness > 0.9
            first = first + 1;
            if first == 1
                here = start;
            end
        else
            if first > 0
                ind1 = here-max_windl;
                ind2 = here+(first-1) + max_windl;
                bjump = (bx(ind1)-bx(ind2))^2+(by(ind1)-by(ind2))^2+(bz(ind1)-bz(ind2))^2;
                if bjump > 0.001
                    inds1(end+1) = ind1;
                    inds2(end+1) = ind2;
                end
            end
            first = 0;
        end
        start = start + 1;
    end

%figure
%plot(ppp)
%axis([0 length(ppp) wll wlh])
end