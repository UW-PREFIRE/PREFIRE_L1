function [cal_only, cal_indices_range, cal_starts, cal_stops] = ...
                               get_cal_data(image_data, encoder_data, varargin)
if(length(varargin)>0)
    cal_en = varargin{1};
else
    cal_en = 3962; %default value
end

if(length(varargin)>1)
    cal_range=varargin{2};
else
    cal_range = 5; %slop for encoder position
end

count = (1:size(image_data,3));
cal_indices = count(and(encoder_data > cal_en - cal_range, encoder_data < cal_en + cal_range));
edges = cal_indices(1);%[]; %try to capture first cal
for i=2:length(cal_indices)
    step = cal_indices(i)-cal_indices(i-1);
    if (step > 2) %find non-consecutive points
        edges = cat(1,edges,cal_indices(i));
    end
end    
cal_indices_range = [edges+1 edges+2 edges+3 edges+4 edges+5 edges+6]; %
cal_starts = edges+1;
cal_stops = edges+6;
cal_indices_range = reshape(cal_indices_range',[1,6*size(edges,1)]);
cal_only = image_data(:,:,cal_indices_range);  

end
