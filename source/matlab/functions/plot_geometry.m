function plot_geometry(gm)%slit_width,pix_width,num_spectral,num_spatial,x_slit,y_slit,tau_slit,x_fp,y_fp,z_fp,tau_fp,tip_fp,tilt_fp,xy_p,slit);
%% angles are in degrees
%% rewritten 8/31/2020 to pass structure
txt = sprintf('slit width %d slit dx=%d, dy=%d tau=%d\n dx=%0.1f,dy=%0.1f,dz=%0.1f,tau=%0.1f,tip=%0.1f,tilt=%0.1f\n',gm.slit_width,gm.x_slit,gm.y_slit,gm.tau_slit,gm.x_fp,gm.y_fp,gm.z_fp,gm.tau_fp,gm.tip_fp,gm.tilt_fp);
figure()
red = [1,0,0];orange = [1,0.5,0];yellow = [1,1,0.1];green = [0,1,0];blue = [0.2 0.2 0.4];violet =[0.6,0,1];
colorOrder =[red;orange;yellow;green;blue;violet];
set(gca,'ColorOrder',colorOrder);

for k=1:4:gm.num_spectral %only show every 4th pixel so the out-of-focus pixels dont always overlap
    for j=1:gm.num_spatial
        X = [gm.xy_p(1,k,j,2) gm.xy_p(1,k,j,3) gm.xy_p(1,k,j,5) gm.xy_p(1,k,j,4)];
        Y = [gm.xy_p(2,k,j,2) gm.xy_p(2,k,j,3) gm.xy_p(2,k,j,5) gm.xy_p(2,k,j,4)];
        pgon = polyshape( X , Y ); 
        nextColor = colorOrder(mod(k,6)+1,:);
        plot(pgon,'FaceColor',nextColor);
        hold on
    end
end
X = [gm.slit(1,2) gm.slit(1,3) gm.slit(1,5) gm.slit(1,4)];
Y = [gm.slit(2,2) gm.slit(2,3) gm.slit(2,5) gm.slit(2,4)];
pgon = polyshape( X , Y ); 
plot(pgon);
%whitebg('black');
xlabel('physical microns');
ylabel('physical microns');
axis([-6000 6000 -6000 6000]);
axis equal %sets equal aspect ratio
title(txt);
hold off
txt = sprintf('');

end
