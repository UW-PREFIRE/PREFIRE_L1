function forward_model_plots(spectral_channels, num_images, channel_center_wavelens, latitudes, sig, radfile);

num_xtrack = num_images(1);
num_intrack = num_images(2);

figure();
x=channel_center_wavelens;
ax = gca;
ax.YAxis.Exponent = 0;
xlabel('wavelength (um)');
ylabel('counts');
hold off;
%plot(channel_center_wavelens, sig(2:64,1,2));
plot(channel_center_wavelens, sig(1:64,1,12));

figure();
%sig2=reshape(sig,spectral_channels+1,num_intrack*num_xtrack);
sig2=reshape(sig,spectral_channels,num_intrack*num_xtrack);
lats = reshape(latitudes,1,num_intrack*num_xtrack);
plot(lats, sig2(2,:),lats, sig2(12,:), lats, sig2(22,:), lats, sig2(32,:),lats, sig2(42,:),lats, sig2(1,:));
xlabel("Latitude");
ylabel("Counts");
title(radfile);

end
