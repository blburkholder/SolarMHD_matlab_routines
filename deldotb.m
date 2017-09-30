where = 5;
ddb = zeros(nx,ny,nz);

%integrating del dot b over a plane
for ix = 2:nx-1
  for iy = 2:ny-1
    ddb(ix,iy) = 1/(z(where-1)+z(where+1))*(bz(ix,iy,where+1) - bz(ix,iy,where-1)) +...
      1/(x(where-1)+x(where+1))*(bx(ix+1,iy,where) - bx(ix-1,iy,where)) +...
      1/(y(where-1)+y(where+1))*(by(ix,iy+1,where) - by(ix,iy-1,where));
  end
end
sum(sum(ddb))
figure
pcolor(ddb)
shading interp
colorbar