r = 0.:1.1276372445109879E-002:11.547005383792516;
z = 0.:9.7656250000000000E-003:10.;
fileI = fopen('gd_3-psi.oa7.txt','r');
out = fscanf(fileI,'%f');
fclose(fileI);
psi = reshape(out,1025,1025);
[C,h] = contourf(r,z,psi,20);
savefig('gd_3-psi.oa7.fig');
fileI = fopen('gd_3-aphi.oa7.txt','r');
out = fscanf(fileI,'%f');
fclose(fileI);
aphi = reshape(out,1025,1025);
[C,h] = contourf(r,z,aphi,20);
savefig('gd_3-aphi.oa7.fig');
