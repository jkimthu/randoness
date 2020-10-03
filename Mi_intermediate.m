%%  Metabolic state Mi, intermediate, 60 min fluctuations



%  Goal: quantify instantaneous growth rate, mu, in response to each pair
%        of nutrient shifts (high, low) as cells are shifted from metabolic
%        state M1 (low-adapted)
%  
%        In this data, single E. coli are adapted to steady low nutrient
%        before a shift to T=60 min fluctuations between high and low.


%  Outputs: 
%        i.  a plot of growth rate vs nutrient phase
%            - growth rates are binned by period fraction (20th of a period)
%              plotting each period in succession, colored by time

%       ii.  a plot of growth rate vs time
%            - growth rates are binned every 75 sec, to match single shift data


%  Strategy:
%
%  Part 1:
%     0. initialize analysis parameters
%     0. initialize complete meta data

%  Part 2:
%     1. for all experiments in dataset:
%           2. initialize experiment meta data
%           3. load measured data
%           4. gather specified condition data
%           5. isolate parameters for growth rate calculations
%           6. calculate growth rate
%           7. generate a period vector by which to isolate data
%           8. isolate growth rate of interest
%           9. for each period, bin growth rates into 20ths of period timescale
%                   10. isolate data from current period,
%                   11. bin growth rates by 20th of period,
%                       assigning growth rate to time value of the middle of two timepoints (not end value)!
%                   12. calculate average growth rate per timebin
%                   13. plot
%    14. repeat for all experiments



%  Last edit: jen, 2020 Oct 2
%  commit: quantify post-shift (high) and post-shift (low) growth rate for Mi


% Okie, go go let's go!

%% Part 0. initialize analysis

clc
clear

% 0. initialize analysis parameters
binsPerPeriod = 30;


% 0. define growth rate of interest
%specificGrowthRate = 'log2';
specificColumn = 3;


% 0. initialize complete meta data
source_data = '/Users/jen/Documents/StockerLab/Source_data';
cd(source_data)
load('storedMetaData.mat')
exptArray = [37,38,39]; % list experiments by index


% 0. initialize colors per successive period
colorSpectrum = {'SlateGray','DarkCyan','CadetBlue','DeepSkyBlue','DodgerBlue','Navy','Indigo',};


%% Part 1. loop through each period and isolate growth rates

% 1. for all experiments in dataset
exptCounter = 0;

for e = 1:length(exptArray)
    
    
    % 2. initialize experiment meta data
    index = exptArray(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    expType = storedMetaData{index}.experimentType;
    bubbletime = storedMetaData{index}.bubbletime;
    xys = storedMetaData{index}.xys;
    shiftTime = storedMetaData{index}.shiftTime;
    
    disp(strcat(date, ': analyze!'))
    exptCounter = exptCounter + 1;

    
    
    % 3. load measured data
    source_data = '/Users/jen/Documents/StockerLab/Source_data';
    cd(source_data)
    filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    load(filename,'D5','T');
    
        
    
    % 4. gather specified condition data
    condition = 1; % 1 = fluctuating; 3 = ave nutrient condition
    if strcmp(expType,'steady2fluc') == 1
            xy_start = xys{condition}(1);
            xy_end = xys{condition}(end);
        else
            xy_start = xys(condition,1);
            xy_end = xys(condition,end);
    end
    conditionData = buildDM(D5, T, xy_start, xy_end, index, expType);
    clear D5 T xy_start xy_end xys
    

    
    % 5. isolate parameters for growth rate calculations
    volumes = getGrowthParameter(conditionData,'volume');             % volume = calculated va_vals (cubic um)
    timestamps_sec = getGrowthParameter(conditionData,'timestamp');   % ND2 file timestamp in seconds
    isDrop = getGrowthParameter(conditionData,'isDrop');              % isDrop == 1 marks a birth event
    curveFinder = getGrowthParameter(conditionData,'curveFinder');    % col 5  = curve finder (ID of curve in condition)
    trackNum = getGrowthParameter(conditionData,'trackNum');          % track number, not ID from particle tracking
    
    
        
    % 6. calculate growth rate
    growthRates = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,trackNum);
    clear trackNum curveFinder isDrop volumes
    
    
    
    % 7. trim timestamps and growth rate by bubble time
    %     NOTE: errors (excessive negative growth rates) occur at trimming
    %           point if growth rate calculation occurs AFTER time trim.
    timestamps_corrected = getGrowthParameter(conditionData,'correctedTime'); % signal corrected time
    if strcmp(date,'2017-10-10') == 1
        timestamps_corrected = timestamps_sec;
    end
        
    timestamps_hr = timestamps_sec/3600;
    maxTime = bubbletime(condition);
    
    if maxTime ~= 0
        timestamps_trimmed = timestamps_corrected(timestamps_hr <= maxTime);
        growthRates_trimmed = growthRates(timestamps_hr <= maxTime,:);
    else
        timestamps_trimmed = timestamps_corrected;
        growthRates_trimmed = growthRates;
    end
    
    
    
    % 8. generate a bin vector for entire time series
    %    assigning growth rate to time value of the middle of two timepoints (not end value)!
    
    % note: replicate three is not at an even timestamp
    % solution: adjust shiftTime and timestamps to help bin correctly

    if strcmp(date,'2019-07-18') == 1
        timestamps_trimmed = timestamps_trimmed + 5*60; % 5 min off of round 3.5 shift time * 60 sec
        shiftTime = shiftTime + 5*60;
    end
    timestep_sec = 60+57;
    timeInSeconds_middle = timestamps_trimmed - (timestep_sec/2);
    timeInPeriods = timeInSeconds_middle/timescale;
    binVector = ceil(timeInPeriods*binsPerPeriod);

    clear timestamps_hr timestamps_sec maxTime growthRates timestamps_corrected
    clear timeInPeriods timeInSeconds_middle timeInPeriods timestamps_trimmed timestep_sec
    
    
  
    % 9. isolate growth rate of interest and bin by timestamp
    growthRt = growthRates_trimmed(:,specificColumn); % log2
    growthRt(binVector == 0) = [];
    binVector(binVector == 0) = [];
    growthRt_binned = accumarray(binVector,growthRt,[],@(x) {x});
    clear growthRates_trimmed

    
    
    % 10. initialize first period, based on shift time (time of fluctuation onset)
    shiftBin = (shiftTime/timescale)*binsPerPeriod;
    period_shift = [shiftBin-1, shiftBin:shiftBin+binsPerPeriod]; % includes two bins before upshift
    period_first = period_shift - binsPerPeriod;
    clear period_shift shiftBin
    
    
    % 11. for each period, bin growth rates into period fraction
    %numPeriods = max(ceil(binVector/binsPerPeriod));
    for pp = 1:7 %numPeriods
        
        period_current = period_first + (30*(pp-1));
        
        
        % 12. isolate data from current period
        isEnd = period_current > length(growthRt_binned);
        if sum(isEnd) > 0
            period_current = period_current(isEnd == 0);
        end
        
        period_mus = growthRt_binned(period_current);

        
   
        % 13. calculate average growth rate per timebin
        signal_current = cellfun(@nanmean,period_mus);
        signals_compiled{pp,:} = signal_current;

        
        
        % 13. plot
        color = rgb(colorSpectrum{pp});
        figure(e)
        plot(signal_current,'Color',color,'LineWidth',1)
        hold on
        grid on
        title('growth rate: log2 mean')
        xlabel(strcat('period bin (',num2str(binsPerPeriod),' binsPerPeriod)'))
        ylabel('growth rate (1/hr)')
        axis([0,binsPerPeriod+3,-1,3])
        title(date)
        
        clear signal_current
    
        
    end
    
    signals_all{exptCounter} = signals_compiled;
    clear signals_compiled
    
    
% 14. repeat for all experiments
end
clear binVector bubbletime color condition conditionData
clear growthRt_binned growthRt isEnd period_current period_first period_mus pp
clear shape shiftTime specificColumn date e filename exptCounter exptArray


%% Part 2. compile replicate data with errorbars

% for each period
time = -1:30;
for period = 1:7
    
    
    % 1. isolate current period signal from both replicate
    rep1 = signals_all{1}{period};
    rep2 = signals_all{2}{period};
    rep3 = signals_all{3}{period};
    
    
    % 2. ensure replicate data are of same length, trim if needed
    %    compile replicate data with columns being time and rows being replicate (period bin)
    if length(rep1) ~= length(rep2)
        rl = [length(rep1); length(rep2)];
        rep_compiled(1,:) = rep1(1:min(rl));
        rep_compiled(2,:) = rep2(1:min(rl));
        rep_compiled(3,:) = rep3(1:min(rl));
    else
        rep_compiled(1,:) = rep1;
        rep_compiled(2,:) = rep2;
        rep_compiled(3,:) = rep3;
    end
    
    
    % 3. calculate mean and stdev across replicates (rows)
    rep_mean = mean(rep_compiled);
    rep_std = std(rep_compiled);
    
    
    % 4. plot
    color = rgb(colorSpectrum{period});
    
    figure(1) % mean of replicate means
    hold on
    plot(time(1:length(rep_mean)),rep_mean,'Color',color,'LineWidth',1,'Marker','.')
    ylabel('instantaneous growth rate (1/h)')
    xlabel('time (min)')
    
    figure(1) % mean and std of replicate means
    hold on
    ss = shadedErrorBar(time(1:length(rep_mean)),rep_compiled,{@nanmean,@nanstd},'lineprops',{'Color',color},'patchSaturation',0.3);

    % set face and edge properties
    ss.mainLine.LineWidth = 3;
    ss.patch.FaceColor = color;
    axis([-2 31 -0.6 2.8])
    
    clear rep_compiled
end


%% Part 3. quantify Mi_low and Mi_high for each of 6 periods
%


%  strategy: from each period signal 
%        
%        Mi_high = growth rate value at t=6 min into period (index 5)
%        Mi_low = growth rate value at t=36 min into period (index 18)
%        

clear index timescale expType specificGrowthRate
clc


% 0. initialize number of replicate experiments, periods and vectors for integrals
numReps = length(signals_all);
numPeriods = 7;
Mi_high = nan(numReps,numPeriods);
Mi_low = nan(numReps,numPeriods);


% 0. initialize bins associated with each integral
bin_high = 5; % t = 6 min
bin_low = 20; % t = 36 min


% 1. loop through each period of each replicate, calculating integrals
for rep = 1:numReps
    
    currentReplicate = signals_all{rep};
    
    for per = 1:numPeriods
        
        currentPeriod = currentReplicate{per};
        Mi_high(rep,per) = currentPeriod(bin_high);  % current Mi_high
        Mi_low(rep,per) = currentPeriod(bin_low); % current Mi_low
       
        clear currentPeriod
        
    end
    
end
clear rep per currentReplicate currentPeriod


% 2. calculate mean and sem of each period
mean_high = mean(Mi_high);
mean_low = mean(Mi_low);

std_high = std(Mi_high);
std_low = std(Mi_low);

sem_high = std_high./sqrt(numReps);
sem_low = std_low./sqrt(numReps);

clear std_high std_low std_total

