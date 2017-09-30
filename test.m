
n = nz-1;
%bx_n = bx./(sqrt(bx.^2+by.^2+bz.^2));
%by_n = by./(sqrt(bx.^2+by.^2+bz.^2));
bz_n = bz./(sqrt(bx.^2+by.^2+bz.^2));
theta = acos(bz_n(2:end-1,2:end-1,n));

%prop = bz(2:end-1,2:end-1,n);
prop = theta;

figure
hold on
colormap('jet')
pcolor(y(2:end-1),x(2:end-1),prop)
shading interp
streamslice(y,x,by(:,:,n),bx(:,:,n))
title('s3 tf')
colorbar

figure
xxxxx=[-2:.001:2],yyyyy=(sqrt(cos(xxxxx)).*cos(200*xxxxx)+sqrt(abs(xxxxx))-0.7).*(4-xxxxx.*xxxxx).^0.01
plot(xxxxx,yyyyy);
figure
penny
why
why
why
figure
spy