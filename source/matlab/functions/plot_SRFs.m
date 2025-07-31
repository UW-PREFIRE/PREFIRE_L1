function plot_SRFs(wl,k1,k2,j1,j2,wavelen,SRF,SRFbyOrder)

red = [1,0,0];orange = [1,0.5,0];yellow = [1,1,0.1];green = [0,1,0];
blue = [0.2 0.2 0.4];violet =[0.6,0,1];grey=[0.5,0.5,0.5];black = [0,0,0];
colorOrder =[red;orange;yellow;green;blue;violet;grey;black];
set(gca,'ColorOrder',colorOrder);

figure();
for k=k1:k2
    for j=j1:j2
        plot(wavelen,SRF(:,k,j));
        hold on;
    end
end
xlabel('wavelength (um)');
ylabel('throughput');    
title(['Total SRF=']);
hold off;

figure();
for k=k1:k2
    for j=j1:j2
        plot(wavelen,SRFbyOrder(:,k,j,3)); 
        hold on;
    end
end
xlabel('wavelength (um)');
ylabel('throughput');   
title(['SRF m=1']);
hold off;

figure();
for k=k1:k2
    for j=j1:j2
        plot(wavelen,SRFbyOrder(:,k,j,2)); 
        hold on;
    end    
end
xlabel('wavelength (um)');
ylabel('throughput');
title(['SRF m=2']);
hold off;

figure();
for k=k1:k2
    for j=j1:j2
        plot(wavelen,SRFbyOrder(:,k,j,1)); 
        hold on;
    end
end
xlabel('wavelength (um)');
ylabel('throughput');
title(['SRF m=3']);
hold off;

%figure();
%plot(wavelen,SRF(:,8,j),wavelen,SRF(:,35,j)); %plotting only jth scene
%xlabel('wavelength (um)');
%ylabel('throughput');
%title(['Total SRF Channels 7 and 34']); %v10 includes channel 0 now

end
