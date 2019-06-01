function [inds1,inds2] = twoDsimulation_walen(x,bx,by,vx,vy,vax,vay,wll,wlh)  
    %wlh = walen limit high
    %wll = walen limit low

    max_windl = 2;
    vxd = zeros(4,1);
    vyd = zeros(4,1);
    vaxd = zeros(4,1);
    vayd = zeros(4,1);

    inds1 = [];
    inds2 = [];

    start = 3;
    first = 0;
    here = 0;

    ppp = zeros(length(vx),1);

    while start <= length(vx) - 3
        vxd(1) = abs(vx(start) - vx(start+1)); vxd(2) = abs(vx(start) - vx(start-1));
        vxd(3) = abs(vx(start) - vx(start+2)); vxd(4) = abs(vx(start) - vx(start-2));
        vyd(1) = abs(vy(start) - vy(start+1)); vyd(2) = abs(vy(start) - vy(start-1));
        vyd(3) = abs(vy(start) - vy(start+2)); vyd(4) = abs(vy(start) - vy(start-2));      
        vaxd(1) = abs(vax(start) - vax(start+1)); vaxd(2) = abs(vax(start) - vax(start-1));
        vaxd(3) = abs(vax(start) - vax(start+2)); vaxd(4) = abs(vax(start) - vax(start-2));
        vayd(1) = abs(vay(start) - vay(start+1)); vayd(2) = abs(vay(start) - vay(start-1));
        vayd(3) = abs(vay(start) - vay(start+2)); vayd(4) = abs(vay(start) - vay(start-2));    

        p = [vxd;vyd]\[vaxd;vayd];
        mslope = p(1)*[vxd;vyd];
        [goodness ,~] = rsquare([vaxd;vayd],mslope);

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
                bjump = (bx(ind1)-bx(ind2))^2+(by(ind1)-by(ind2))^2;
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