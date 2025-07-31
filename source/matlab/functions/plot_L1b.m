function plot_L1b(filepath,latitudes,longitudes,time,spec_rad);
num_images = size(latitudes);
%txt = sprintf('slit width %d slit dx=%d, dy=%d tau=%d\n dx=%0.1f,dy=%0.1f,dz=%0.1f,tau=%0.1f,tip=%0.1f,tilt=%0.1f\n',slit_width,x_slit,y_slit,tau_slit,x_fp,y_fp,z_fp,tau_fp,tip_fp,tilt_fp);
txt = sprintf(erase(filepath,'\'));
figure()
plot(longitudes,latitudes);
xlabel('longitude');
ylabel('latitude');
title(txt);

figure()
plot(latitudes,reshape(spec_rad(10,:,:),num_images));
xlabel('latitude');
ylabel('channel 10');
title(txt);

figure()
for i=1:num_images(1)
    plot(time,latitudes(i,:));
    hold on;
end
xlabel('time (seconds)');
ylabel('latitude');
title(txt);
hold off;

% below is leftover from geometry plotter, will maybe be useful for
% plotting footprints

red = [1,0,0];orange = [1,0.5,0];yellow = [1,1,0.1];green = [0,1,0];blue = [0.2 0.2 0.4];violet =[0.6,0,1];
colorOrder =[red;orange;yellow;green;blue;violet];
set(gca,'ColorOrder',colorOrder);

for k=1:num_images(1)
    for m=1:num_images(2)
        %save this for corners when L1b gets more populated
        %X = [xy_p(1,k,j,2) xy_p(1,k,j,3) xy_p(1,k,j,5) xy_p(1,k,j,4)];
        %Y = [xy_p(2,k,j,2) xy_p(2,k,j,3) xy_p(2,k,j,5) xy_p(2,k,j,4)];
        %pgon = polyshape( X , Y ); 
        %nextColor = colorOrder(mod(k,6)+1,:);
        %plot(pgon,'FaceColor',nextColor);
        %hold on
    end
end
%whitebg('black');
%xlabel('longitude');
%ylabel('latitude');
%axis([-6000 6000 -6000 6000]);
%axis equal %sets equal aspect ratio
%title(txt);
%hold off
txt = sprintf('');

end
