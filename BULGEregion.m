%fname = {'2dsims/magtap0080p1n2','2dsims/magtap0065p1b2','2dsims/magtap0065p1','2dsims/magtap0080p.1'};
fname = {'2dsims/gf_runs/magtap0070_p1b2by','2dsims/gf_runs/magtap0070_p1n2by','2dsims/gf_runs/magtap0070_p4by','2dsims/gf_runs/magtap0075_p1by','2dsims/gf_runs/magtap0080_p.1by'};
figure; hold on
for kkk = 1
    [x,y,rho,p,vx,vy,vz,bx,by,bz] = twoD_read(fname{kkk});

    upperbulgecut = 700;
    lowerbulgecut = 550;
    jetcut = 300;

    jx = zeros(size(bx));
    jy = zeros(size(bx));
    jz = zeros(size(bx));

    jx(2:end-1,:) = (bz(1:end-2,:) - bz(3:end,:))/(y(2)-y(1))/2;
    jy(:,2:end-1) = -(bz(:,1:end-2) - bz(:,3:end))/(y(2)-y(1))/2;
    jz(:,2:end-1) = (by(:,1:end-2) - by(:,3:end))/(x(2)-x(1))/2;
    jz(2:end-1,:) =  jz(2:end-1,:)-(bx(1:end-2,:) - bx(3:end,:))/(y(2)-y(1))/2;
% 
    figure; hold on
    pcolor(x,y,sqrt(jx.^2+jy.^2+jz.^2))
    plot([x(127) x(627)],[y(upperbulgecut) y(upperbulgecut)],'w--','LineWidth',2)
    plot([x(127) x(627)],[y(lowerbulgecut) y(lowerbulgecut)],'w-','LineWidth',2)
    plot([x(127) x(627)],[y(jetcut) y(jetcut)],'w:','LineWidth',2)
    shading interp
    %[xx,yy] = meshgrid(x,y);
    h = streamslice(x,y,bx,by);
    set(h,'Color','k')
    daspect([1 1 1])
    title('|j| and magnetic field')
    xlabel('x')
    ylabel('y')
    colorbar
    axis([x(1) x(end) y(1) y(end)])
    caxis([0 2])
% 
    figure; hold on
    pcolor(x,y,p)
    shading interp
    colormap(hot)
    res = 50;
    [xx,yy] = meshgrid(x,y);
    quiver(xx(1:res:end,1:res:end),yy(1:res:end,1:res:end),vx(1:res:end,1:res:end),vy(1:res:end,1:res:end))
    daspect([1 1 1])
    title('thermal pressure and flow')
    xlabel('x')
    ylabel('y')
    colorbar
    axis([x(1) x(end) y(1) y(end)])

    skinterval = 127:627;
   figure
   linespeccers = {'k:','k-','k--'};
%     if kkk == 4
%         linespeccers = {'r-','r-','r-'};
%     elseif kkk == 5
%         linespeccers = {'b-','b-','b-'};
%     else
%         linespeccers = {'g-','g-','g-'};
%     end

    ccutline = [jetcut,lowerbulgecut,upperbulgecut];
    tittys = [{'jet'},{'lower bulge'},{'middle bulge'}];

    for hh = 1:3
        cutline = ccutline(hh);

%         subplot(5,3,hh); hold on;
%         plot(x(skinterval),by(cutline,skinterval),linespeccers{hh},'LineWidth',2)
%         title(tittys(hh))
%         if hh == 1
%             ylabel('B_y');
%         else
%             set(gca,'YTickLabel',[])     
%         end
%         xlim([-3*hh 0])
%         ylim([-1.5 1.5])
%         set(gca,'XTickLabel',[])
%         set(gca,'FontSize',16)
%         grid on
% 
%         subplot(5,3,hh+3); hold on;
%         plot(x(skinterval),bz(cutline,skinterval),linespeccers{hh},'LineWidth',2)
%         if hh == 1
%             ylabel('B_z');
%         else
%             set(gca,'YTickLabel',[])    
%         end
%         xlim([-3*hh 0])
%         ylim([0.4 1])
%         set(gca,'XTickLabel',[])
%         set(gca,'FontSize',16)
%         grid on
% 
%         subplot(5,3,hh+6); hold on;
%         rotation = 180*acos((bz(cutline,1)*bz(cutline,:)+by(cutline,1)*by(cutline,:))./(sqrt(bz(cutline,1)^2+by(cutline,1)^2)*sqrt(bz(cutline,:).^2+by(cutline,:).^2)))/pi;
%         plot(x(skinterval),rotation(skinterval),linespeccers{hh},'LineWidth',2)
%         if hh == 1
%             ylabel('rotation');
%             set(gca,'YTick',[0 20 40 60])
%         else
%             set(gca,'YTick',[0 20 40 60])
%             set(gca,'YTickLabel',[])    
%         end
%         xlim([-3*hh 0])
%         ylim([0 60])
%         set(gca,'XTickLabel',[])
%         set(gca,'FontSize',16)
%         grid on
% 
%         subplot(5,3,hh+9); hold on;
%         plot(x(skinterval),rho(cutline,skinterval),linespeccers{hh},'LineWidth',2)
%         if hh == 1
%             ylabel('density');
%         else
%             set(gca,'YTickLabel',[])    
%         end
%         xlim([-3*hh 0])
%         set(gca,'XTickLabel',[])
%         set(gca,'FontSize',16)
%         ylim([0.5 4])
%         grid on
% 
%         subplot(5,3,hh+12); hold on;
%         entrop = p./(rho.^(5/3));
%         plot(x(skinterval),entrop(cutline,skinterval)/entrop(cutline,1),linespeccers{hh},'LineWidth',2)
%         if hh == 1
%             ylabel('S/S_0');
%         else
%             set(gca,'YTickLabel',[])    
%         end
%         xlabel('x')
%         xlim([-3*hh 0])
%         set(gca,'FontSize',16)
%         ylim([0.9 2])
%         grid on


        subplot(3,2,1); hold on
        plot(x(skinterval),by(cutline,skinterval),linespeccers{hh})
        title('B_y')

         subplot(3,2,2); hold on
%         plot(x(skinterval),abs(by(cutline,skinterval))./sqrt(rho(cutline,skinterval)),linespeccers{hh})
%         title('v_{Ay}')
        rotation = 180*acos((bz(cutline,1)*bz(cutline,:)+by(cutline,1)*by(cutline,:))./(sqrt(bz(cutline,1)^2+by(cutline,1)^2)*sqrt(bz(cutline,:).^2+by(cutline,:).^2)))/pi;
        subplot(3,2,2); hold on
        plot(x(skinterval),rotation(skinterval),linespeccers{hh})
        title('rotation')

        subplot(3,2,3); hold on
        plot(x(skinterval),vy(cutline,skinterval),linespeccers{hh})
        title('v_y')

        subplot(3,2,4); hold on
        plot(x(skinterval),rho(cutline,skinterval),linespeccers{hh})
        title('density')

        subplot(3,2,5); hold on
        plot(x(skinterval),p(cutline,skinterval),linespeccers{hh})
        title('thermal pressure')
        xlabel('x')

        subplot(3,2,6); hold on
        plot(x(skinterval),p(cutline,skinterval)./rho(cutline,skinterval).^(5/3),linespeccers{hh})
        title('entropy')
        xlabel('x')

        [inds1,inds2] = twoDsimulation_walen(x(skinterval),bx(cutline,skinterval),by(cutline,skinterval),...
            vx(cutline,skinterval),vy(cutline,skinterval),...
            abs(bx(cutline,skinterval))./sqrt(rho(cutline,skinterval)),...
            abs(by(cutline,skinterval))./sqrt(rho(cutline,skinterval)),0.9,1.1);

        subplot(3,2,3)
        for j = 1:length(inds1)
            plot(x((126+inds1(j)):(126+inds2(j))),vy(cutline,(126+inds1(j)):(126+inds2(j))),'y','LineWidth',2)
        end    
    end
end

