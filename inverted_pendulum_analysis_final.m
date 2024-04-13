%% Load Data
clear all; 
close all; 
clc; 
% List for filepaths
filePaths1 = {
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S1_D1_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S2_D1_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S3_D1_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S4_D1_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S5_D1_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S6_D1_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S7_D1_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S8_D1_data.mat';
};

% Initialize scoreMat1
scoreMat1 = zeros(length(filePaths1), 5);

for i = 1:length(filePaths1)
    % Load data from the file
    load(filePaths1{i}, 'data');

    % Extract MaxPoints from the loaded data
    maxPoints1 = data.MaxPoints;

    % Sort MaxPoints based on FBLocation
    FBLocation = data.FBLocation;
    [~, sortIdx] = sortrows([FBLocation, (1:length(FBLocation))', maxPoints1]);
    maxPoints1 = maxPoints1(sortIdx);

    % Calculate mean for each pendulum length and store in scoreMat1
    col = 1;
    for j = 1:80:400
        scoreMat1(i, col) = mean(maxPoints1(j:j + 79));
        col = col + 1;
    end
end


% Calculate mean and standard deviation across columns of scoreMat
MeanScore1 = mean(scoreMat1); 
STDScore1 = std(scoreMat1);
PendLengths = [0.75, 1, 1.5, 2, 4];
display(scoreMat1)
% Create figure for shaded area
figure;
X = [PendLengths fliplr(PendLengths)];
Y = [MeanScore1+STDScore1 fliplr(MeanScore1-STDScore1)];
fill(X, Y, [0, 0.5, 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
hold on;
plot(PendLengths, MeanScore1, '-', 'Color', [0, 0.5, 0], 'MarkerSize', 8);

xlabel('Pendulum Lengths');
ylabel('Mean Score');
title('Mean Score & Standard Deviation ~ Pendulum Length');
grid on;

hold off;

%% % List for filepaths
filePaths2 = {
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S1_D2_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S2_D2_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S3_D2_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S4_D2_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S5_D2_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S6_D2_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S7_D2_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S8_D2_data.mat';
};

% Initialize scoreMat1
scoreMat2 = zeros(length(filePaths2), 5);

for i = 1:length(filePaths2)
    % Load data from the file
    load(filePaths2{i}, 'data');

    % Extract MaxPoints from the loaded data
    maxPoints2 = data.MaxPoints;

    % Sort MaxPoints based on FBLocation
    FBLocation = data.FBLocation;
    [~, sortIdx] = sortrows([FBLocation, (1:length(FBLocation))', maxPoints2]);
    maxPoints2 = maxPoints2(sortIdx);

    % Calculate mean for each pendulum length and store in scoreMat1
    col = 1;
    for j = 1:80:400
        scoreMat2(i, col) = mean(maxPoints1(j:j + 79));
        col = col + 1;
    end
end


% Calculate mean and standard deviation across columns of scoreMat
MeanScore2 = mean(scoreMat2); 
STDScore2 = std(scoreMat2);
PendLengths = [0.75, 1, 1.5, 2, 4];
display(scoreMat2)
% Create figure for shaded area
figure;
X2 = [PendLengths fliplr(PendLengths)];
Y2 = [MeanScore2+STDScore2 fliplr(MeanScore2-STDScore2)];
fill(X2, Y2, [0, 0.5, 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
hold on;
plot(PendLengths, MeanScore2, '-', 'Color', [0, 0.5, 0], 'MarkerSize', 8);

xlabel('Pendulum Lengths');
ylabel('Mean Score');
title('Mean Score & Standard Deviation ~ Pendulum Length');
grid on;

hold off;

%% %% Load Data
clear all; 
close all; 
clc; 

% List for filepaths Day 1
filePaths1 = {
    'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S1_D1_data.mat';
    'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S2_D1_data.mat';
    'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S3_D1_data.mat';
    'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S4_D1_data.mat';
    'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S5_D1_data.mat';
    'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S6_D1_data.mat';
    'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S7_D1_data.mat';
    'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S8_D1_data.mat';
};

% Initialize scoreMat1
scoreMat1 = zeros(length(filePaths1), 5);

for i = 1:length(filePaths1)
    % Load data from the file
    load(filePaths1{i}, 'data');

    % Extract MaxPoints from the loaded data
    maxPoints1 = data.MaxPoints;

    % Sort MaxPoints based on FBLocation
    FBLocation = data.FBLocation;
    [~, sortIdx] = sortrows([FBLocation, (1:length(FBLocation))', maxPoints1]);
    maxPoints1 = maxPoints1(sortIdx);

    % Calculate mean for each pendulum length and store in scoreMat1
    col = 1;
    for j = 1:80:400
        scoreMat1(i, col) = mean(maxPoints1(j:j + 79));
        col = col + 1;
    end
end

% Calculate mean and standard deviation across columns of scoreMat for Day 1
MeanScore1 = mean(scoreMat1); 
STDScore1 = std(scoreMat1);
PendLengths = [0.75, 1, 1.5, 2, 4];

% List for filepaths Day 2
filePaths2 = {
    'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S1_D2_data.mat';
    'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S2_D2_data.mat';
    'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S3_D2_data.mat';
    'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S4_D2_data.mat';
    'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S5_D2_data.mat';
    'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S6_D2_data.mat';
    'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S7_D2_data.mat';
    'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S8_D2_data.mat';
};

% Initialize scoreMat2
scoreMat2 = zeros(length(filePaths2), 5);

for i = 1:length(filePaths2)
    % Load data from the file
    load(filePaths2{i}, 'data');

    % Extract MaxPoints from the loaded data
    maxPoints2 = data.MaxPoints;

    % Sort MaxPoints based on FBLocation
    FBLocation = data.FBLocation;
    [~, sortIdx] = sortrows([FBLocation, (1:length(FBLocation))', maxPoints2]);
    maxPoints2 = maxPoints2(sortIdx);

    % Calculate mean for each pendulum length and store in scoreMat2
    col = 1;
    for j = 1:80:400
        scoreMat2(i, col) = mean(maxPoints2(j:j + 79));
        col = col + 1;
    end
end

% Calculate mean and standard deviation across columns of scoreMat for Day 2
MeanScore2 = mean(scoreMat2); 
STDScore2 = std(scoreMat2);

% Create figure for shaded area
figure;
hold on;

% Plot for Day 1
X = [PendLengths fliplr(PendLengths)];
Y = [MeanScore1+STDScore1 fliplr(MeanScore1-STDScore1)];
plot(PendLengths, MeanScore1, '-', 'Color', [0.4940 0.1840 0.5560], 'MarkerSize', 8);
fill(X, Y, [0.4940 0.1840 0.5560], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
% Plot for Day 2
X2 = [PendLengths fliplr(PendLengths)];
Y2 = [MeanScore2+STDScore2 fliplr(MeanScore2-STDScore2)];
plot(PendLengths, MeanScore2, '-', 'Color', 'r', 'MarkerSize', 8);
fill(X2, Y2, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
xlabel('Pendulum Lengths');
ylabel('Mean Score');
title('Mean Score & Standard Deviation ~ Pendulum Length');
grid on;
legend('Day 1','Day 1 STD' ,'Day 2', 'Day 2 STD');

hold off;


