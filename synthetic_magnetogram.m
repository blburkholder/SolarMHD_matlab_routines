function [a,frac,cmn1,cmn2] = synthetic_magnetogram(x,y,cmn1,cmn2,xmodes,ymodes,bz,flux_control)
[cx,cy] = meshgrid(x,y-y(end-1)/2);

a = zeros(length(x),length(y));

if flux_control == 0
    for m = 1:xmodes
        for n = 1:ymodes
            am = pi*(m-1)/x(end-1);
            bn = pi*(n-1)/y(end-1);

            modd = mod(n+m-2,2);
            if modd == 0
                a = a + cmn1(m,n)*sin(cx*am).*sin(cy*bn) +...
                    cmn2(m,n)*cos(cx*am).*cos(cy*bn);
            end
        end
    end
end

%great plain of gaussians
if flux_control == 1
    nspots = 100 + floor(100*rand(1));
    randos = zeros(nspots,4);
    randos(:,1) = 10*randn(nspots,1);
    randos(:,2) = max(x)*rand(nspots,1);
    randos(:,3) = max(y)*rand(nspots,1)-max(y)/2;
    randos(:,4) = 10*rand(nspots,1);
    for i = 1:nspots
        a = a + randos(i,1)*exp((-(cx-randos(i,2)).^2-(cy-randos(i,3)).^2)/randos(i,4));
    end
end


%%%%for normalizing total flux to observation
[~,tot] = flux(a,x(2)-x(1),y(2)-y(1));
[~,tot1] = flux(bz,x(2)-x(1),y(2)-y(1));
frac = tot/tot1;

