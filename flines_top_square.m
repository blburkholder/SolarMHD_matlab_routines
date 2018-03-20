load('mentrop0.mat')

%property which will paint bottom/top boundary
prop = ftvol;

%prop quantity maximum in case of extraneous values
maxi = 400;
%where to start flines
n = nz-1;

count = 75;
[row,col] = size(prop);
bs = zeros(count,2);

center1 = [50,30];
center2 = [7,35];
center3 = [40,78];
center4 = [15,10];

bs1(1,:) = center1;
bs2(1,:) = center2;
bs3(1,:) = center3;
bs4(1,:) = center4;

xspread1 = 10;
yspread1 = 10;
xspread2 = 2;
yspread2 = 2;
xspread3 = 5;
yspread3 = 10;
xspread4 = 8;
yspread4 = 8;

for i = 2:count
    bs1(i,:) = center1 + [2*xspread1*rand(1,1)-xspread1,2*yspread1*rand(1,1)-yspread1];
    bs2(i,:) = center2 + [2*xspread2*rand(1,1)-xspread2,2*yspread2*rand(1,1)-yspread2];
    bs3(i,:) = center3 + [2*xspread3*rand(1,1)-xspread3,2*yspread3*rand(1,1)-yspread3];
    bs4(i,:) = center4 + [2*xspread4*rand(1,1)-xspread4,2*yspread4*rand(1,1)-yspread4];
end

prop(prop > maxi) = maxi;
prop(prop < -maxi) = -maxi;
%prop = prop/max(max(prop));

%figure
%hold on
%pcolor(xpo,ypo,prop)
%shading interp
%colormap(gray)
%colorbar
%scatter(bs1(:,1),bs1(:,2),[20],'m.')
%scatter(bs2(:,1),bs2(:,2),[30],'r.')
% scatter(bs3(:,1),bs3(:,2),[20],'g.')
% scatter(bs4(:,1),bs4(:,2),[20],'b.')
%daspect([1 1 1])
%axis([0,x(end-1),0,y(end-1)])

load('mentrop0_b.mat')
prop = zfin1;
maxi = 10
prop(prop > maxi) = maxi;
%prop = prop/max(max(prop));
figure
hold on
pcolor(xpo,ypo,prop)
shading interp
%colormap(gray)
colorbar
streamslice(x,y,sx(:,:,2)./rho(:,:,2),sy(:,:,2)./rho(:,:,2),'y')

plines1 = perio_stream3(x(2:end-1),y(2:end-1),z(2:end-1),...
    bx(2:end-1,2:end-1,2:end-1),by(2:end-1,2:end-1,2:end-1),...
    bz(2:end-1,2:end-1,2:end-1),bs1(:,1),bs1(:,2),...
    z(n)*ones(length(bs1(:,1)),1));
% plines2 = perio_stream3(x(2:end-1),y(2:end-1),z(2:end-1),...
%     bx(2:end-1,2:end-1,2:end-1),by(2:end-1,2:end-1,2:end-1),...
%     bz(2:end-1,2:end-1,2:end-1),bs2(:,1),bs2(:,2),...
%     z(n)*ones(length(bs2(:,1)),1));
plines3 = perio_stream3(x(2:end-1),y(2:end-1),z(2:end-1),...
    bx(2:end-1,2:end-1,2:end-1),by(2:end-1,2:end-1,2:end-1),...
    bz(2:end-1,2:end-1,2:end-1),bs3(:,1),bs3(:,2),...
    z(n)*ones(length(bs3(:,1)),1));
plines4 = stream3(x(2:end-1),y(2:end-1),z(2:end-1),...
    -bx(2:end-1,2:end-1,2:end-1),-by(2:end-1,2:end-1,2:end-1),...
    -bz(2:end-1,2:end-1,2:end-1),bs4(:,1),bs4(:,2),...
    z(n)*ones(length(bs4(:,1)),1));
[~,plines44] = perio_stream32(x(2:end-1),y(2:end-1),z(2:end-1),...
    bx(2:end-1,2:end-1,2:end-1),by(2:end-1,2:end-1,2:end-1),...
    bz(2:end-1,2:end-1,2:end-1),bs4(:,1),bs4(:,2),...
    z(n)*ones(length(bs4(:,1)),1));

%this will put wittle red stars at the conjugate footpoints
% conj1 = zeros(length(plines1),3);
% conj2 = zeros(length(plines1),3);
% conj3 = zeros(length(plines1),3);
% conj4 = zeros(length(plines1),3);

% for l = 1:length(plines1)
%     if ~isempty(plines1{l})
%         conj1(l,1) = plines1{l}(end,1);
%         conj1(l,2) = plines1{l}(end,2);
%         conj1(l,3) = plines1{l}(end,3);
% 
%         conj2(l,1) = plines2{l}(end,1);
%         conj2(l,2) = plines2{l}(end,2);
%         conj2(l,3) = plines2{l}(end,3);
% 
%         conj3(l,1) = plines3{l}(end,1);
%         conj3(l,2) = plines3{l}(end,2);
%         conj3(l,3) = plines3{l}(end,3);
%     
%         conj4(l,1) = plines4{l}(end,1);
%         conj4(l,2) = plines4{l}(end,2);
%         conj4(l,3) = plines4{l}(end,3);
%     end
% end
%scatter(conj1(:,1),conj1(:,2),'m.')
%streamline(plines1)
%scatter(conj2(:,1),conj2(:,2),'r.')
%scatter(conj3(:,1),conj3(:,2),'g.')
%scatter(conj4(:,1),conj4(:,2),'y.')

daspect([1 1 1])

colormap(bone)
axis off
colorbar off
%axis([0,x(end-1),0,y(end-1)])

h = streamline(plines1);
set(h,'Color',[1 0 1])
h = streamline(plines3);
set(h,'Color',[0 1 0])
h = streamline(plines4);
set(h,'Color',[1 1 0])
h = streamline(plines44);
set(h,'Color',[1 1 0])
view(3)
