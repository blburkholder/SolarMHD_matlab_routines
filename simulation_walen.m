function [ind1,ind2] = simulation_walen(x,rho,bx,by,bz,vx,vy,vz,vax,vay,vaz,wll,wlh)  
    %wlh = walen limit high
    %wll = walen limit low
    ind1 = 0;
    ind2 = 0;

    px = zeros(length(vx),1);
    max_windl = 4;
    wind_inc = 1;
    vxd = zeros(max_windl/wind_inc,1);
    vyd = zeros(max_windl/wind_inc,1);
    vzd = zeros(max_windl/wind_inc,1);
    vaxd = zeros(max_windl/wind_inc,1);
    vayd = zeros(max_windl/wind_inc,1);
    vazd = zeros(max_windl/wind_inc,1);

    inds1 = [];
    inds2 = [];

    start = max_windl+1;
    first = 0;
    pmax = 0;
    here = 0;
    heremax = 0;
    rsqr = 0;
    while start <= length(vx) - max_windl
        for wind = wind_inc:wind_inc:max_windl
            vxd(wind/wind_inc) = abs(vx(start-wind) - vx(start+wind)); 
            vyd(wind/wind_inc) = abs(vy(start-wind) - vy(start+wind));     
            vzd(wind/wind_inc) = abs(vz(start-wind) - vz(start+wind));     
            vaxd(wind/wind_inc) = abs(vax(start-wind) - vax(start+wind));     
            vayd(wind/wind_inc) = abs(vay(start-wind) - vay(start+wind));     
            vazd(wind/wind_inc) = abs(vaz(start-wind) - vaz(start+wind));
        end
        p = polyfit([vxd;vyd;vzd],[vaxd;vayd;vazd],1);
        f = polyval(p,[vxd;vyd;vzd]);
        [goodness ,~] = rsquare([vaxd;vayd;vazd],f);

        px(start) = p(1);
        if (p(1) >= wll) && (p(1) <= wlh)
            first = first + 1;
            if first == 1
                here = start;
                rsqr = 0;
            end
            start = start+max_windl;
            rsqr = rsqr + goodness;
        else
            if first > 0 && rsqr/first >  0.5
                ind1 = here-max_windl/2;
                ind2 = here+(first+1)*max_windl/2;
                bjump = (bx(ind1)-bx(ind2))^2+(by(ind1)-by(ind2))^2+(bz(ind1)-bz(ind2))^2;
                if bjump > 0.01
                    inds1(end+1) = ind1;
                    inds2(end+1) = ind2;
                end
            end
            first = 0;
            start = start+1;
        end
    end

    figure
    subplot(5,1,1); hold on
    plot(x,bx,'r')
    plot(x,by,'g')
    plot(x,bz,'b')
    title('magnetic field')
    set(gca,'xticklabel',[])
    legend('B_x','B_y','B_z')
    axis([-2 2 -1.1 1.1])
    subplot(5,1,2);
    plot(x,rho)
    title('density')
    set(gca,'xticklabel',[])
    axis([-2 2 0.8 1.5])
    subplot(5,1,5); hold on
    for j = 1:2
        for i = inds1(j):inds2(j)
            plot([x(i) x(i)],[-2 2],'y')
        end
    end
    p1 = plot(x,vax,'r');
    p2 = plot(x,vay,'g');
    p3 = plot(x,vaz,'b');
    title('Alfven velocity')
    xlabel('x')
    legend([p1,p2,p3],'v_{Ax}','v_{Ay}','v_{Az}')
    axis([-2 2 -0.1 1.1])
    subplot(5,1,3); hold on
    for j = 1:2
        for i = inds1(j):inds2(j)
            plot([x(i) x(i)],[-2 2],'y')
        end
    end
    p1 = plot(x,vx,'r');
    p2 = plot(x,vy,'g');
    p3 = plot(x,vz,'b');
    title('velocity')
    legend([p1,p2,p3],'v_x','v_y','v_z')
    axis([-2 2 -1 1])
    set(gca,'xticklabel',[])
end