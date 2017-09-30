bbx = bx(2:end-1,2:end-1,2:end-1);
bby = by(2:end-1,2:end-1,2:end-1);
bbz = bz(2:end-1,2:end-1,2:end-1);
may_null = zeros(size(bbx)-1);

t_8x = zeros(8,1);
t_8y = zeros(8,1);
t_8z = zeros(8,1);

xx = x(2:end-1);
yy = y(2:end-1);
zz = z(2:end-1);

ddx = xx(2) - xx(1);
ddy = yy(2) - yy(1);

bvec = [0;0;0];
b2 = 0;
jays = 0;

for ix = 1:nx-3
    for iy = 1:ny-3
        for iz = 1:nz-3
            t_8y = [bbx(ix,iy,iz) bbx(ix+1,iy,iz) bbx(ix,iy+1,iz)...
                bbx(ix,iy,iz+1) bbx(ix+1,iy+1,iz) bbx(ix+1,iy,iz+1)...
                bbx(ix,iy+1,iz+1) bbx(ix+1,iy+1,iz+1)];
            t_8x = [bby(ix,iy,iz) bby(ix+1,iy,iz) bby(ix,iy+1,iz)...
                bby(ix,iy,iz+1) bby(ix+1,iy+1,iz) bby(ix+1,iy,iz+1)...
                bby(ix,iy+1,iz+1) bby(ix+1,iy+1,iz+1)];
            t_8z = [bbz(ix,iy,iz) bbz(ix+1,iy,iz) bbz(ix,iy+1,iz)...
                bbz(ix,iy,iz+1) bbz(ix+1,iy+1,iz) bbz(ix+1,iy,iz+1)...
                bbz(ix,iy+1,iz+1) bbz(ix+1,iy+1,iz+1)];

            hx = t_8x;
            hy = t_8y;
            hz = t_8z;

            hx(hx < 0) = -1;
            hx(hx > 0) = 1;
            hy(hy < 0) = -1;
            hy(hy > 0) = 1;
            hz(hz < 0) = -1;
            hz(hz > 0) = 1;
 
            if abs(sum(hx)) < 8 && abs(sum(hy)) < 8 && abs(sum(hz)) < 8
                %doesnt always return an integer so something is wrong
                %td = topo_deg(t_8x,t_8y,t_8z);
                td = 1;
                if td == 1 || td == -1                
                    may_null(ix,iy,iz) = 1;   
                    %if we have a null do newton's method!
                    ddz = zz(iz+1) - zz(iz);

                    % coefficients for trilinear interpolation
                    ax = t_8x(1);
                    b_x = t_8x(2) - t_8x(1);
                    cx = t_8x(3) - t_8x(1);
                    dx = t_8x(5) - t_8x(2) - t_8x(3) + t_8x(1);
                    ex = t_8x(4) - t_8x(1);
                    fx = t_8x(6) - t_8x(2) - t_8x(4) + t_8x(1);
                    gx = t_8x(7) - t_8x(3) - t_8x(4) + t_8x(1);
                    hx = t_8x(8) - t_8x(5) - t_8x(6) - t_8x(7) + t_8x(2) + t_8x(3) + t_8x(4) - t_8x(1);

                    ay = t_8y(1);
                    b_y = t_8y(2) - t_8y(1);
                    cy = t_8y(3) - t_8y(1);
                    dy = t_8y(5) - t_8y(2) - t_8y(3) + t_8y(1);
                    ey = t_8y(4) - t_8y(1);
                    fy = t_8y(6) - t_8y(2) - t_8y(4) + t_8y(1);
                    gy = t_8y(7) - t_8y(3) - t_8y(4) + t_8y(1);
                    hy = t_8y(8) - t_8y(5) - t_8y(6) - t_8y(7) + t_8y(2) + t_8y(3) + t_8y(4) - t_8y(1);

                    az = t_8z(1);
                    b_z = t_8z(2) - t_8z(1);
                    cz = t_8z(3) - t_8z(1);
                    dz = t_8z(5) - t_8z(2) - t_8z(3) + t_8z(1);
                    ez = t_8z(4) - t_8z(1);
                    fz = t_8z(6) - t_8z(2) - t_8z(4) + t_8z(1);
                    gz = t_8z(7) - t_8z(3) - t_8z(4) + t_8z(1);
                    hz = t_8z(8) - t_8z(5) - t_8z(6) - t_8z(7) + t_8z(2) + t_8z(3) + t_8z(4) - t_8z(1);

                    steps = 20;
                    b2l = 1000;
                    jl = 0;
                    %9 start points, middle of cell, middle of 8 sub cells    
                    for k = 1:9
                        if k == 1
                            stx = ddx/2;
                            sty = ddy/2;    
                            stz = ddz/2;
                        else
                            binar = dec2bin(k-2,3);
                            stx = ddx/4 + str2num(binar(1))*ddx/2;
                            sty = ddy/4 + str2num(binar(2))*ddy/2;
                            stz = ddz/4 + str2num(binar(3))*ddz/2;
                        end

                        dyz = ddy*ddz;
                        dxz = ddx*ddz;
                        dxy = ddx*ddy;
                        dxyz = ddx*ddy*ddz;

                        j = 0;                    
                        while j < steps
                            %jacobian matrix
                            jac = [...
                                b_x/ddx + dx*sty/dxy + fx*stz/dxz + hx*sty*stz/dxyz,...
                                b_y/ddx + dy*sty/dxy + fy*stz/dxz + hy*sty*stz/dxyz,...   
                                b_z/ddx + dz*sty/dxy + fz*stz/dxz + hz*sty*stz/dxyz;...

                                cx/ddy + dx*stx/dxy + gx*stz/dyz + hx*stx*stz/dxyz,...
                                cy/ddy + dy*stx/dxy + gy*stz/dyz + hy*stx*stz/dxyz,...
                                cz/ddy + dz*stx/dxy + gz*stz/dyz + hz*stx*stz/dxyz;...

                                ex/ddz + fx*stx/dxz + gx*sty/dyz + hx*stx*sty/dxyz,...
                                ey/ddz + fy*stx/dxz + gy*sty/dyz + hy*stx*sty/dxyz,...
                                ez/ddz + fz*stx/dxz + gz*sty/dyz + hz*stx*sty/dxyz];

                            bvec = [ax + b_x*stx/ddx + cx*sty/ddy + dx*stx*sty/(ddx*ddy) + ex*stz/ddz + fx*stx*stz/(ddx*ddz) + gx*sty*stz/(ddy*ddz) + hx*stx*sty*stz/(ddx*ddy*ddz);...
                                ay + b_y*stx/ddx + cy*sty/ddy + dy*stx*sty/(ddx*ddy) + ey*stz/ddz + fy*stx*stz/(ddx*ddz) + gy*sty*stz/(ddy*ddz) + hy*stx*sty*stz/(ddx*ddy*ddz);...
                                az + b_z*stx/ddx + cz*sty/ddy + dz*stx*sty/(ddx*ddy) + ez*stz/ddz + fz*stx*stz/(ddx*ddz) + gz*sty*stz/(ddy*ddz) + hz*stx*sty*stz/(ddx*ddy*ddz)];
                            if sum(bvec.^2) < 0.00000001
                                break;
                            else
                                j = j + 1;
                            end

                            %constrained line search that I made up
%                             alph = 1;
%                             while (true)
%                                 s = [stx;sty;stz] - alph*(jac\bvec);
%                                 if s(1) <= ddx && s(1) >= 0 && s(2) <= ddy && s(2) >= 0 && s(3) <= ddz && s(3) >= 0
%                                     bvec2 = [ax + b_x*s(1)/ddx + cx*s(2)/ddy + dx*s(1)*s(2)/(ddx*ddy) + ex*s(3)/ddz + fx*s(1)*s(3)/(ddx*ddz) + gx*s(2)*s(3)/(ddy*ddz) + hx*s(1)*s(2)*s(3)/(ddx*ddy*ddz);...
%                                         ay + b_y*s(1)/ddx + cy*s(2)/ddy + dy*s(1)*s(2)/(ddx*ddy) + ey*s(3)/ddz + fy*s(1)*s(3)/(ddx*ddz) + gy*s(2)*s(3)/(ddy*ddz) + hy*s(1)*s(2)*s(3)/(ddx*ddy*ddz);...
%                                         az + b_z*s(1)/ddx + cz*s(2)/ddy + dz*s(1)*s(2)/(ddx*ddy) + ez*s(3)/ddz + fz*s(1)*s(3)/(ddx*ddz) + gz*s(2)*s(3)/(ddy*ddz) + hz*s(1)*s(2)*s(3)/(ddx*ddy*ddz)];
%                                     
%                                     f1 = sum(bvec2.^2);
%                                     if sum(bvec.^2) < f1
%                                         %this means not descent direction
%                                         s(1) = ddx;
%                                         break
%                                     end
% 
%                                     alph = alph*0.9;
%                                     s = [stx;sty;stz] - alph*(jac\bvec);
% 
%                                     bvec2 = [ax + b_x*s(1)/ddx + cx*s(2)/ddy + dx*s(1)*s(2)/(ddx*ddy) + ex*s(3)/ddz + fx*s(1)*s(3)/(ddx*ddz) + gx*s(2)*s(3)/(ddy*ddz) + hx*s(1)*s(2)*s(3)/(ddx*ddy*ddz);...
%                                         ay + b_y*s(1)/ddx + cy*s(2)/ddy + dy*s(1)*s(2)/(ddx*ddy) + ey*s(3)/ddz + fy*s(1)*s(3)/(ddx*ddz) + gy*s(2)*s(3)/(ddy*ddz) + hy*s(1)*s(2)*s(3)/(ddx*ddy*ddz);...
%                                         az + b_z*s(1)/ddx + cz*s(2)/ddy + dz*s(1)*s(2)/(ddx*ddy) + ez*s(3)/ddz + fz*s(1)*s(3)/(ddx*ddz) + gz*s(2)*s(3)/(ddy*ddz) + hz*s(1)*s(2)*s(3)/(ddx*ddy*ddz)];
%                                     f2 = sum(bvec2.^2);
% 
%                                     if f1 < f2 || f1 == f2 || alph < 0.0000001
%                                         break;
%                                     end                                    
%                                 else
%                                     alph = alph*0.9;
%                                 end
%                             end
%                                         

                            %or take "relaxed" newton step to satisfy wolfe    
                            s = [stx;sty;stz] - 0.1*(jac\bvec);

                            %does step remain within cell?
                            if s(1) > ddx || s(1) < 0 || s(2) > ddy || s(2) < 0 || s(3) > ddz || s(3) < 0
                                %dont accept if first step goes outside cell
                                %if j == 1
                                    bvec = [nan,nan,nan];
                                %end
                                    
                                break;
                            else
                                stx = s(1);
                                sty = s(2);
                                stz = s(3);
                            end 
                        end
                        %keep the starting point with smallest min(B^2)
                        if sum(bvec.^2) < b2l
                            b2l = sum(bvec.^2);
                            jl = j;
                        elseif b2l == 1000
                            b2l = nan;
                        end
                    end
                    b2(end+1) = (b2l);
                    jays(end+1) = jl;
                end
            end
        end
    end
end

figure
scatter(1:length(b2),b2,[],jays,'filled')
set(gca,'yscale','log')
colorbar

% figure
% pcolor(sum(there_is_a_null_here,3))
% shading interp

% f = figure;
% ax = axes('Parent',f);
% iz = 1;
% h = pcolor(ax,x(1:end-3),y(1:end-3),there_is_a_null_here(:,:,iz)');
% shading interp
% %colorbar
% b = uicontrol('Parent',f,'Style','slider','Position',[81,0,419,23],...
%               'value',iz, 'min',1, 'max',nz-3);
% set(b, 'SliderStep', [1/(nz-3) , 10/(nz-3) ]);
% b.Callback  = @(objHandle,~)  update_plot3(x(1:end-3),y(1:end-3),there_is_a_null_here,get(objHandle,'Value'))



