function [] = oneDsimulatorgator(fin)
    [x,rho,p,temp,vx,vy,vz,bx,by,bz] = oneD_read(fin);     

    figure
    subplot(2,2,1); hold on
    plot(x,bx(:,1),'r')
    plot(x,by(:,1),'g')
    plot(x,bz(:,1),'b')
    plot(x,p(:,1)./(bx(:,1).^2+by(:,1).^2+bz(:,1).^2),'k')
    title('magnetic field')
    legend('B_x','B_y','B_z','\beta')
    axis([-2 2 -1.5 2])
    set(gca,'xticklabel',[])

    subplot(2,2,3); hold on
    plot(x,p(:,1))
    plot(x,p(:,1)./rho(:,1).^(5/3))
    legend('p','S')
    title('Thermal Pressure and Specific Entropy')
    xlabel('x')
    axis([-2 2 0 3])
    
%     subplot(3,2,2); hold on
%     plot(x,abs(bx(:,1))./sqrt(rho(:,1)),'r')
%     plot(x,abs(by(:,1))./sqrt(rho(:,1)),'g')
%     plot(x,abs(bz(:,1))./sqrt(rho(:,1)),'b')
%     title('Alfven velocity')
%     legend('v_{Ax}','v_{Ay}','v_{Az}')
%     %axis([-2 2 0.0 2.5])
%     set(gca,'xticklabel',[])

    subplot(2,2,2); hold on
    plot(x,vx(:,1),'r')
    plot(x,vy(:,1),'g')
    plot(x,vz(:,1),'b')
    title('velocity')
    set(gca,'xticklabel',[])
    legend('v_x','v_y','v_z','location','southeast')
    axis([-2 2 -0.5 1.5])

    subplot(2,2,4)
    plot(x,rho(:,1))
    title('density')
    axis([-2 2 0.7 1.1])
    xlabel('x')

    [inds1,inds2] = simulation_walen(x,bx(:,end),by(:,end),bz(:,end),vx(:,end),vy(:,end),vz(:,end),...
            abs(bx(:,end))./sqrt(rho(:,end)),abs(by(:,end))./sqrt(rho(:,end)),abs(bz(:,end))./sqrt(rho(:,end)),0.9,1.1);

    figure
    subplot(2,2,1); hold on
    plot(x,bx(:,end),'r')
    plot(x,by(:,end),'g')
    plot(x,bz(:,end),'b')
    plot(x,p(:,end)./(bx(:,end).^2+by(:,end).^2+bz(:,end).^2),'k')
    title('magnetic field')
    set(gca,'xticklabel',[])
    legend('B_x','B_y','B_z','\beta')
    axis([-2 2 -1.5 2])
    subplot(2,2,3); hold on
    plot(x,p(:,end))
    plot(x,p(:,end)./rho(:,end).^(5/3))
    legend('p','S')
    title('Thermal Pressure and Specific Entropy')
    xlabel('x')
    axis([-2 2 0.0 3])
%     subplot(2,2,2); hold on
% %     for j = 1:1:length(inds1)
% %         for i = inds1(j):inds2(j)
% %             plot([x(i) x(i)],[0 2.5],'y')
% %         end
% %     end
%     p1 = plot(x,abs(bx(:,end))./sqrt(rho(:,end)),'r');
%     p2 = plot(x,abs(by(:,end))./sqrt(rho(:,end)),'g');
%     p3 = plot(x,abs(bz(:,end))./sqrt(rho(:,end)),'b');
%     title('Alfven velocity')
%     legend([p1,p2,p3],'v_{Ax}','v_{Ay}','v_{Az}')
%     %axis([-2 2 0.0 2.5])
%     set(gca,'xticklabel',[])
     subplot(2,2,2); hold on
    for j = 1:length(inds1)
        for i = inds1(j):inds2(j)
            plot([x(i) x(i)],[-1.0 0.2],'y')
        end
    end
    p1 = plot(x,vx(:,end),'r');
    p2 = plot(x,vy(:,end),'g');
    p3 = plot(x,vz(:,end),'b');
    title('velocity')
    legend([p1,p2,p3],'v_x','v_y','v_z','location','southeast')
    axis([-2 2 -0.5 1.5])
    set(gca,'xticklabel',[])
    subplot(2,2,4)
    plot(x,rho(:,end))
    title('density')
    axis([-2 2 0.7 1.1])
    xlabel('x')

%     subplot(3,2,5)
%     plot(x,p(:,end))
%     title('thermal pressure')
%     %axis([-2 2 0 2])
%     xlabel('x')
end