wx = [-flipud(x(4:end));x];
wwx = [-flipud(xpo(4:end));xpo];
wftv = [fliplr(flipud(ftvol(:,3:end))),ftvol(:,2:end)];
wvx = [-fliplr(flipud(sx(:,3:end,:)./rho(:,3:end,:))),sx(:,2:end,:)./rho(:,2:end,:)];
wvy = [-fliplr(flipud(sy(:,3:end,:)./rho(:,3:end,:))),sy(:,2:end,:)./rho(:,2:end,:)];
wvz = [-fliplr(flipud(sz(:,3:end,:)./rho(:,3:end,:))),sz(:,2:end,:)./rho(:,2:end,:)];
wbx = [-fliplr(flipud(bx(:,3:end,:))),bx(:,2:end,:)];
wby = [-fliplr(flipud(by(:,3:end,:))),by(:,2:end,:)];
wbz = [fliplr(flipud(bz(:,3:end,:))),bz(:,2:end,:)];
wrho = [fliplr(flipud(rho(:,3:end,:))),rho(:,2:end,:)];

wftv(wftv > 150) = 150;

nnnn = nulls(end-5:end-2,:);
nnn = [mean(nnnn(:,1)),mean(nnnn(:,2)),mean(nnnn(:,3))];

[bsx,bsy,bsz] = meshgrid(nnn(2),y,z); 

nflines = 3;
ss = 0.9;
vvvec = (-ss:2*ss/nflines:ss);
[xx1,yy1,zz1] = meshgrid(nnn(2)+vvvec,nnn(1)+vvvec,nnn(3)+vvvec);

figure
% pcolor(wx,y,wvx(:,:,2).^2+wvy(:,:,2).^2)
% shading interp
daspect([1 1 1])
view([0 0])
title('Perturbation at fan footpoint')
%axis([-15 50 10 70 0 35])

plines1 = stream3(wx,y,z,wbx,wby,wbz,xx1,yy1,zz1); %spine
plines2 = stream3(wx,y,z,-wbx,-wby,-wbz,xx1,yy1,zz1); %fan

ep2 = zeros(length(plines1),1);
for i = 1:length(plines1)
    xyzfin = plines2{i}(end,:);
    ep2(i) = interp3(wx,y,z,(wvx.^2+wvy.^2+wvz.^2)./wrho,xyzfin(1),xyzfin(2),xyzfin(3)); 
end

colormap(parula)
mep2 = max(ep2);
cm = colormap;

for i = 1:length(plines1)
    jgz = floor(63*ep2(i)/mep2)+1;
    h = streamline(plines1(i));
    set(h,'Color',cm(jgz,:));
%     h = streamline(plines2(i));
%     set(h,'Color',cm(jgz,:))
end

colormap(parula)
axis off

% xyzfin = zeros(length(plines1),3);
% for i = 1:length(plines1)
%     xyzfin(i,:) = plines1{i}(end,:); 
% end
% 
% figure
% scatter(xyzfin(:,1),xyzfin(:,2),[],ep2/mep2,'.')

figure
ep1 = zeros(length(plines1),1);
for i = 1:length(plines1)
    xyzfin = plines1{i}(end,:);
    ep1(i) = interp3(wx,y,z,(wvx.^2+wvy.^2+wvz.^2)./wrho,xyzfin(1),xyzfin(2),xyzfin(3)); 
end

colormap(parula)
mep1 = max(ep1);
cm = colormap;

for i = 1:length(plines1)
    jgz = floor(63*ep1(i)/mep1)+1;
     h = streamline(plines1(i));
    set(h,'Color',cm(jgz,:)) 
    h = streamline(plines2(i));
    set(h,'Color',cm(jgz,:))
end

daspect([1 1 1])
view(3)
title('Perturbation at spine footpoint')
