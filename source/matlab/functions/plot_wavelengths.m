function plot_wavelengths(wl)%slit_width,pix_width,num_spectral,num_spatial,x_slit,y_slit,tau_slit,x_fp,y_fp,z_fp,tau_fp,tip_fp,tilt_fp,dd,wavelengths);
%% rewritten 8/31/2020 to pass structure
red = [1,0,0];orange = [1,0.5,0];yellow = [1,1,0.1];green = [0,1,0];blue = [0.2 0.2 0.4];violet =[0.6,0,1];
ideal_wavelength = wl.dd.*[0:wl.num_spectral-1];
ch=[1:64];
txt = sprintf('slit width %d slit dx=%d, dy=%d tau=%d\n FP dx=%0.1f,dy=%0.1f,dz=%0.1f,tau=%0.1f,tip=%0.1f,tilt=%0.1f\n',wl.slit_width,wl.x_slit,wl.y_slit,wl.tau_slit,wl.x_fp,wl.y_fp,wl.z_fp,wl.tau_fp,wl.tip_fp,wl.tilt_fp);
figure()
colorOrder =[red;orange;yellow;green;blue;violet];
set(gca,'ColorOrder',colorOrder);
hold on
for j=1:wl.num_spatial
    c = colorOrder(mod(j,6)+1,:);
    %wavelengths are stored in units of 10nm but we want to plot in um
    y1=wl.wavelengths(:,j,2)/100-ideal_wavelength';
    y2=wl.wavelengths(:,j,1)/100-ideal_wavelength';
    y3=wl.wavelengths(:,j,3)/100-ideal_wavelength';
    plot(ch,y1,'Color',c,'LineWidth',2.0);%,ch,y2,'Color',c,ch,y3,'Color',c)
    plot(ch,y2,'Color',c,'LineWidth',2.0);
    plot(ch,y3,'Color',c,'LineWidth',2.0);
end
%whitebg('black');
xlabel('channel');
ylabel('spectral sampling (microns)');
axis([0 65 -5 5]);
title(txt);
txt = sprintf('');
hold off

end
