function [fig] = quicklook_L1A_data( ...
    obs, obs_time, calseq, gain_at_obs, space_at_obs, ...
    cal_radiance, channel_number, scene_number)

fig = figure('Position', [0 750 600 750]);

% assuming input channel number starts at 0.
c = channel_number + 1;
s = scene_number;
obsDN = squeeze(obs(c,s,:));
spaceDN = squeeze(space_at_obs(c,s,:));
cal_rad = squeeze(cal_radiance(c,s,:));
gain = squeeze(gain_at_obs(c,s,:));

space_seq = squeeze(calseq.space_seqmean_v(:,c,s));
bbcal_seq = squeeze(calseq.cal_seqmean_v(:,c,s));

% we will plot all the DN time series with floating y-axis ranges,
% but we do this after offseting to be zero mean. The required
% offsets are then displayed in the annotations.
spaceDN_offset = mean(space_seq, 'all');
rounded_spaceDN_offset = round(spaceDN_offset, -2);

obsDN_offset = mean(obsDN, 'all');
rounded_obsDN_offset = round(obsDN_offset, -2);

bbcalDN_offset = mean(bbcal_seq, 'all');
rounded_bbcalDN_offset = round(bbcalDN_offset, -2);

% this looks like extra work - but hopefully makes it clearer as to
% how things map to the annotations (ylabels and legend)
rounded_deltaDN_obs = rounded_obsDN_offset - rounded_spaceDN_offset;
rounded_deltaDN_bbcal = rounded_bbcalDN_offset - rounded_spaceDN_offset;

obsDN_offset     =     obsDN - rounded_spaceDN_offset - rounded_deltaDN_obs;
spaceDN_offset   =   spaceDN - rounded_spaceDN_offset;
space_seq_offset = space_seq - rounded_spaceDN_offset;
bbcal_seq_offset = bbcal_seq - rounded_spaceDN_offset - rounded_deltaDN_bbcal;

% relative times
obs_rt = obs_time - obs_time(1);
space_rt = calseq.space_seqmean_ts - obs_time(1);
bbcal_rt = calseq.contrast_seqmean_ts - obs_time(1);

% 15 minute padding on either side - not sure this is the best
% method, we want to see "some" of the calseq around the obs data.
time_range = [obs_rt(1) - 900.0, obs_rt(end) + 900.0];

if rounded_spaceDN_offset < 0
    ylabel_string = ['Raw DN + ' num2str(-rounded_spaceDN_offset)];
else
    ylabel_string = ['Raw DN - ' num2str(rounded_spaceDN_offset)];
end

% raw DN in the Earth observation frames
subplot(5,1,1);
plot(obs_rt, obsDN_offset, '.g');
if rounded_deltaDN_obs < 0
    legend({['Obs DN + ' num2str(-rounded_deltaDN_obs)]})
else
    legend({['Obs DN - ' num2str(rounded_deltaDN_obs)]})
end
ylabel(ylabel_string);
xlim(time_range);

% raw DN in the blackbody cal target frames
subplot(5,1,2);
plot(bbcal_rt, bbcal_seq_offset, 'sr', 'MarkerSize', 12);
if rounded_deltaDN_bbcal < 0
    legend({['Cal DN + ' num2str(-rounded_deltaDN_bbcal)]})
else
    legend({['Cal DN - ' num2str(rounded_deltaDN_bbcal)]})
end
ylabel(ylabel_string);
xlim(time_range);

% raw DN in the space view frames
subplot(5,1,3);
plot(obs_rt, spaceDN_offset, '.b');
hold on;
plot(space_rt, space_seq_offset, 'sb', 'MarkerSize', 12);
hold off;
ylabel(ylabel_string);
legend({'Space DN', 'Space Seq DN'}, 'NumColumns', 3);
ylabel(ylabel_string);
xlim(time_range);

% gain derived from cal seq
subplot(5,1,4);
plot(obs_rt, gain, '.c');
legend({'Gain'});
ylabel('Gain [DN / RU]')
xlim(time_range);

% final calibrated radiance
subplot(5,1,5);
plot(obs_rt, cal_rad, '.m');
legend({'Calibrated Radiance'});
ylabel('RU [W / (m^2 sr um)]');
xlabel('Time since start of observations [s]')
xlim(time_range);

end
