%% Metabolic state M1, low-adapted
%
%  Goal: quantify instantaneous growth rate, mu, in response to single nutrient downshift
%
%        In this data, single E. coli are adapted to steady high nutrient
%        before a shift to low nutrient.


%  Output: 
%
%        i.  Two plots of growth rate response to upshifts in real time.
%             (1) mean of each experimental replicate, one curve each
%             (2) mean + standard dev across experimental replicates

%       ii.  Also outputs data as two vectors, growth rate and time,
%            in a .mat workspace titled, quantified_m2.mat


%  General strategy:
%
%         Part 0. initialize folder with stored meta data
%         Part 1. plot and align replicate single upshift data
%         Part 2. quantify M1_low and M1_high


%  last updated: jen, 2020 Oct 2
%  commit: quantify pre-shift (high) and post-shift (low) growth rate for M2


% OK let's go!

%% Part 0. initialize

clc
clear

% 0. initialize complete meta data
source_data = '/Users/jen/Documents/StockerLab/Source_data';
cd(source_data)
load('storedMetaData.mat')


% 0. define shift type, growth rate and time bin of interest
shiftType = 'downshift';
specificGrowthRate = 'log2';
specificColumn = 3; % log2 growth rate
xmin = -0.5;
xmax = 3.5;
timePerBin = 75;

counter = 0; % because timescales will differ between experiments
color_single = 'DarkCyan'; % color designations


%% Part 1A. curves for single shift data


% 1. create array of experiments of interest, then loop through each
exptArray = [26,27]; % downshift data


for e_shift = 1:length(exptArray)
    
    counter = counter + 1;
    
    % 2. initialize experiment meta data
    index = exptArray(e_shift); 
    date = storedMetaData{index}.date;
    
    % define which frames to ignore (noisy tracking)
    if strcmp(date,'2018-08-09') == 1
        ignoredFrames = [115,116,117];
    elseif strcmp(date,'2018-08-08') == 1
        ignoredFrames = [112,113,114];
    end
    
    timescale = storedMetaData{index}.timescale;
    timescale_vector(counter) = timescale;
    
    bubbletime = storedMetaData{index}.bubbletime;
    expType = storedMetaData{index}.experimentType;
    shiftTime = storedMetaData{index}.shiftTime;        % sec
    disp(strcat(date, ': analyze!'))
    
    
    
    % 3. load measured data
    source_data = '/Users/jen/Documents/StockerLab/Source_data';
    cd(source_data)
    filename = strcat('lb-fluc-',date,'-width1p7-jiggle-0p5.mat');
    load(filename,'D5','T')
    
    
    
    % 4. specify condition of interest (fluc) and build data matrix
    condition = 1;                                      % 1 = fluctuating
    xy_start = storedMetaData{index}.xys(condition,1);
    xy_end = storedMetaData{index}.xys(condition,end);
    conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
    
    
    
    % 5. isolate volume (Va), timestamp, mu, drop and curveID data
    volumes = getGrowthParameter(conditionData,'volume');             % calculated va_vals (cubic um)
    timestamps_sec = getGrowthParameter(conditionData,'timestamp');   % timestamp in seconds
    isDrop = getGrowthParameter(conditionData,'isDrop');              % isDrop, 1 marks a birth event
    curveFinder = getGrowthParameter(conditionData,'curveFinder');    % curve finder (ID of curve in condition)
    trackNum = getGrowthParameter(conditionData,'trackNum');          % track number (not ID from particle tracking)
    clear xy_start xy_end
    
    
    
    % 6. calculate growth rate
    growthRates = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,trackNum);
    clear curveFinder trackNum isDrop volumes
    
    
    
    % 7. isolate data to stabilized regions of growth
    %    NOTE: errors (excessive negative growth rates) occur at trimming
    %          point if growth rate calculation occurs AFTER time trim.
    
    minTime = 2.5; % for single shift data, unlike original fluc
    maxTime = bubbletime(condition);
    timestamps_sec = conditionData(:,2); % time in seconds converted to hours
    timestamps_hr = timestamps_sec / 3600;
    clear condition
    
    % trim to minumum
    times_trim1 = timestamps_hr(timestamps_hr >= minTime);
    conditionData_trim1 = conditionData(timestamps_hr >= minTime,:);
    growthRates_trim1 = growthRates(timestamps_hr >= minTime,:);
    
    % trim to maximum
    if maxTime > 0
        conditionData_trim2 = conditionData_trim1(times_trim1 <= maxTime,:);
        growthRates_trim2 = growthRates_trim1(times_trim1 <= maxTime,:);
    else
        conditionData_trim2 = conditionData_trim1;
        growthRates_trim2 = growthRates_trim1;
    end
    clear growthRates conditionData maxTime minTime timestamps_sec timestamps_hr

    
     
    % 8. isolate selected specific growth rate
    growthRt = growthRates_trim2(:,specificColumn);
    % specificColumn is already defined in Part B.
    % not re-defining it here ensures that we use the same metric between both
    clear growthRates_trim1 growthRates_trim2
    

    
    % 9. isolate corrected timestamp
    correctedTime = conditionData_trim2(:,22); % col 22 = timestamps corrected for signal lag
    clear D5 T isDrop conditionData_trim1
    
    
    
    % 10. assign NaN to all growth rates associated with frames to ignore
    frameNum = conditionData_trim2(:,16); % col 16 = original frame number
    growthRt_ignorant = growthRt;
    for fr = 1:length(ignoredFrames)
        growthRt_ignorant(frameNum == ignoredFrames(fr),1) = NaN;
    end
    clear fr
    
    
    
    % 11. remove nans from data analysis
    growthRt_noNaNs = growthRt_ignorant(~isnan(growthRt_ignorant),:);
    correctedTime_noNans = correctedTime(~isnan(growthRt_ignorant),:);
    clear growthRt growthRt_ignorant correctedTime frameNum
    
   
    
    % 12. assign corrected timestamps to bins, by which to accumulate growth data
    bins = ceil(correctedTime_noNans/timePerBin);      % bin 1 = first timePerBin sec of experiment
    bins_unique = (min(bins):max(bins))';              % avoids missing bins due to lack of data
    
    
    % find bins after shift
    % generalized for single shift experiments
    first_postshiftBin_single = ceil(shiftTime/timePerBin) + 1; % shift occurs at shiftTime/timePerBin, so first full bin with shifted data is the next one
    postshiftBins_single{counter} = (first_postshiftBin_single:max(bins))';
    
    
    
    % 13. choose which pre-shift data bins to plot
    % single shift experiments don't have high/low phase interruptions
    % however, they do have bins missing data!
    numPreshiftBins = 10;
    preShift_bins{counter} = numPreshiftBins;

    % determine pre-shift bins
    index_single_shift = find(bins_unique == first_postshiftBin_single);
    pre_downshiftBins{counter} = bins_unique(index_single_shift-numPreshiftBins : index_single_shift-1);
    

    
    % 14. collect growth rate data into bins and calculate stats
    %     WARNING: bins variable does not i
    binned_growthRate{counter} = accumarray(bins,growthRt_noNaNs,[],@(x) {x});
    binned_mean{counter} = accumarray(bins,growthRt_noNaNs,[],@mean);
    clear bins
    
    
    
    % 15. plot response in growth rate for all timescales over time
    timePerBin_min = timePerBin/60; % time per bin in minutes
    color = rgb(color_single);
    
    
    % in single shift data, not all time bins have data!
    % plot accordingly to avoid sharp drops
    preDownshift_times = ((numPreshiftBins-1)*-1:0)*timePerBin_min;
    postDownshift_times = (1:length(binned_mean{counter}(postshiftBins_single{counter})))*timePerBin_min;
    downshift_times_gapped = [preDownshift_times, postDownshift_times];
    
    preDownshift_growth_gapped = binned_mean{counter}(pre_downshiftBins{counter});
    postDownshift_growth_single = binned_mean{counter}(postshiftBins_single{counter});
    downshift_growth_gapped = [preDownshift_growth_gapped; postDownshift_growth_single];
    
    
    % don't plot zeros that are place holders for gaps in data
    downshift_growth = downshift_growth_gapped(downshift_growth_gapped > 0);
    downshift_times = downshift_times_gapped(downshift_growth_gapped > 0);
    
    
    figure(1) % mean of each replicate
    plot(downshift_times,downshift_growth,'Color',color,'LineWidth',1)
    grid on
    hold on
    title(strcat('response to downshift, binned every (',num2str(timePerBin),') sec'))
    
    
    % store data for shaded error bars
    binned_singles{counter} = downshift_growth_gapped;
    binned_singles_times{counter} = downshift_times_gapped;
    

    % save mean signal (of replicate that reaches G_low)
    if strcmp(date,'2018-08-09') == 1
        source_data = '/Users/jen/Documents/StockerLab/Source_data';
        cd(source_data)
        save('response_singleDownshift.mat','downshift_growth','downshift_times')
    end

     
end
xlabel('time (min)')
ylabel(strcat('growth rate: (', specificGrowthRate,')'))
axis([numPreshiftBins*-1*timePerBin_min,160,xmin,xmax])

clear e_shift date
clear replicate_means replicate_means_mean replicate_means_std


%% Part 1B. find mean and standard deviation between single-shift replicates


% extract 2d matrix of mean replicate data, with columns replicate and rows being time (period bin)
for rep = 1:2 % two replicates for upshift and downshift as of 2018-12-04
    col = rep;
    
    curr_rep_means = binned_singles{col}(1:149);
    curr_rep_means(curr_rep_means == 0) = NaN; % make zeros NaN prior to averaging between replicates
    replicate_single_means(rep,:) = curr_rep_means;
    replicate_single_times(rep,:) = binned_singles_times{col}(1:149); % demonstrates indeces are of the same timeline
    
end
clear rep    


% calculate mean and stdev across replicates (rows)
replicate_single_means_mean = mean(replicate_single_means);
replicate_single_means_std = std(replicate_single_means);



% don't plot zeros that are place holders for gaps in data
downshift_means_single = replicate_single_means_mean(replicate_single_means_mean > 0);
downshift_stds = replicate_single_means_std(replicate_single_means_mean > 0);
downshift_times_single = downshift_times_gapped(replicate_single_means_mean > 0);

downshift_replicate_means = replicate_single_means(:,replicate_single_means_mean > 0);


% plot
figure(2) % mean of replicate means
hold on
plot(downshift_times_single,downshift_means_single,'Color',color,'LineWidth',1,'Marker','.')
title('reponse to downshift: mean of replicate means')
ylabel('growth rate: log2 (1/hr)')
xlabel('time (min)')

figure(3) % mean and std of replicate means
hold on
ss = shadedErrorBar(downshift_times_single,downshift_replicate_means,{@nanmean,@nanstd},'lineprops',{'Color',color},'patchSaturation',0.3);
title('response to downshift, mean and std of replicates')

% set face and edge properties
ss.mainLine.LineWidth = 3;
ss.patch.FaceColor = color;
axis([-5,120,-0.5,3.4])

% save only single shift data
figure(3)
ylabel('instaneous growth rate (1/h)')
xlabel('time (min)')


% save mean signal (across replicates)
M2_growthrate = downshift_means_single;
M2_time = downshift_times_single;

repo = '/Users/jen/randoness';
cd(repo)
save('quantified_m2.mat','M2_growthrate','M2_time')


%% Part 2. quantify M1_low and M1_high

clc
clear
load('quantified_m2.mat')

% 1. M1_high is steady-state low growth rate
%    average data points before shift
preshift_t = find(M2_time < 0);
preshift_G = M2_growthrate(preshift_t);
M2_high = mean(preshift_G);

% 2. M1_low is first growth rate after shift (t=0)
postshift_t = M2_time(M2_time > 0);
postshift_t = postshift_t(postshift_t < 35);
postshift_G = M2_growthrate(M2_time > 0);
postshift_G = postshift_G(postshift_t < 35);

M2_low = mean(postshift_G);

