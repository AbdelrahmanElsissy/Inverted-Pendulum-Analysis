% Specify the folder path
folderPath = 'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData'; 
% Get a list of MAT files in the folder
matFiles = dir(fullfile(folderPath, '*.mat'));

% Initialize variables to store data
allMaxPoints = [];

% Loop through each MAT file
for fileIndex = 1:length(matFiles)
    filePath = fullfile(folderPath, matFiles(fileIndex).name);

    load(filePath, 'data');

    % Extract FBlocation and Maxpoints from the loaded data
    maxPoints = data.MaxPoints;
    FBLocation = data.FBLocation;
    [~, sortIdx] = sortrows([FBLocation, (1:length(FBLocation))', maxPoints]);
    maxPoints = maxPoints(sortIdx);

    allMaxPoints = [allMaxPoints, maxPoints];

    col = 1;
    for i = 1:80:400
        scoreMat(fileIndex, col) = mean(maxPoints(i:i + 79, 1));
        col = col + 1;
    end
end
display(scoreMat)
pendLengths = [0.75 1 1.5 2 4]
PendLengths = [0.75, 1, 1.5, 2, 4];
figure;
plot(PendLengths, mean(scoreMat))

% Calculate mean and standard deviation
% meanMaxPoints = mean(allMaxPoints, 2);
% stdMaxPoints = std(allMaxPoints, 0, 2);

% Calculate mean and standard deviation across columns of scoreMat
MeanScore = mean(scoreMat); 
STDScore = std(scoreMat);

% Create figure for shaded area
figure;
X = [PendLengths fliplr(PendLengths)];
Y = [MeanScore+STDScore fliplr(MeanScore-STDScore)];
fill(X, Y, [0.4940 0.1840 0.5560], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
hold on;
plot(PendLengths, MeanScore, '-', 'Color', [0.4940 0.1840 0.5560], 'MarkerSize', 8);

xlabel('Pendulum Lengths');
ylabel('Mean Score');
title('Mean Score & Standard Deviation ~ Pendulum Length');
grid on;


hold off;
%% 

for fileIndex = 1:length(matFiles)
    filePath = fullfile(folderPath, matFiles(fileIndex).name);

    load(filePath, 'data');
    allFallflag = [];
    % Extract FallFlag and FBLocation from the data structure
    FallFlag = data.FallFlag; % Ensure that the field name is correct
    FBLocation = data.FBLocation; % Ensure that the field name is correct
    
    % Sort the data based on FBLocation, maintaining the order of FallFlag
    [~, sortIdx] = sortrows([FBLocation, (1:length(FBLocation))', FallFlag]);
    FallFlag = FallFlag(sortIdx);
    
    % Store the MaxPoints data in the allFallflag matrix
    allFallFlag(:, fileIndex) = FallFlag;
    % allFallflag = [allFallFlag, FallFlag];

    % Block of 80 = one pendulum length and calculate one mean of it
     col = 1;
        for i = 1:80:400
            FallMat(fileIndex, col) = sum(FallFlag(i:i+79,1));
            col=col+1;
        end
end

PendLengths = [0.75, 1, 1.5, 2, 4];
figure;
plot(PendLengths, mean(FallMat))

% Calculate mean and standard deviation
meanFallFlag = mean(allFallFlag, 2);
stdFallFlag = std(allFallFlag, 0, 2);
% Calculate mean and standard deviation across columns of fallflags
MeanFall = mean(FallMat); 
STDFall = std(FallMat);

% Create figure for shaded area
figure;
X = [PendLengths fliplr(PendLengths)];
Y = [MeanFall+STDFall fliplr(MeanFall-STDFall)];
fill(X, Y, [0.4940 0.1840 0.5560], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
hold on;
plot(PendLengths, MeanFall, '-', 'Color', [0.4940 0.1840 0.5560], 'MarkerSize', 8);
xlabel('Pendulum Lengths');
ylabel('Number of Falls');
title('Number of Falls ~ Pendulum Length');
grid on;

% Are subjects improving through blocks?

% Initialize variables to store data
allMaxPoints2 = [];
% Loop through each MAT file
for fileIndex = 1:length(matFiles)
    filePath = fullfile(folderPath, matFiles(fileIndex).name);

    load(filePath, 'data');

    % Extract FBlocation and Maxpoints from the loaded data
    maxPoints = data.MaxPoints;
    FBLocation = data.FBLocation;
    [~, sortIdx] = sortrows([FBLocation, (1:length(FBLocation))', maxPoints]);
    maxPoints = maxPoints(sortIdx);

    allMaxPoints2 = [allMaxPoints2, maxPoints];

    col = 1;
    for i = 1:40:400
        % Calculate mean for the current block
        scoreMat2(fileIndex, col) = mean(maxPoints(i:i + 39, 1));

        col = col + 1;
    end
end
display(scoreMat2)
% Separate the data into two blocks based on odd and even indices
block1Scores = scoreMat2(:, 1:2:end);  % Odd indices
block2Scores = scoreMat2(:, 2:2:end);  % Even indices

% Calculate mean and standard deviation for each block
meanBlock1 = mean(block1Scores);
stdBlock1 = std(block1Scores);

meanBlock2 = mean(block2Scores);
stdBlock2 = std(block2Scores);

% Create figure for shaded area
figure;
% Plot for Block 1
X1 = [PendLengths fliplr(PendLengths)];
Y1 = [meanBlock1 + stdBlock1 fliplr(meanBlock1 - stdBlock1)];
fill(X1, Y1, [0.4940 0.1840 0.5560], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
hold on;
plot(PendLengths, meanBlock1, '-', 'Color', [0.4940 0.1840 0.5560], 'MarkerSize', 8);

% Plot for Block 2
X2 = [PendLengths fliplr(PendLengths)];
Y2 = [meanBlock2 + stdBlock2 fliplr(meanBlock2 - stdBlock2)];
fill(X2, Y2, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(PendLengths, meanBlock2, '-', 'Color', 'r', 'MarkerSize', 8);

xlabel('Pendulum Lengths');
ylabel('Mean Score');
title('Mean Score & Standard Deviation ~ Pendulum Length');
legend('Block 1 ± Std','Block 1 Mean', 'Block 2 ± Std', 'Block 2 Mean');
grid on;

hold off;

% Loop through each MAT file
for fileIndex = 1:length(matFiles)
    filePath = fullfile(folderPath, matFiles(fileIndex).name);

    load(filePath, 'data');

    % Extract velocities from the loaded data
    velocities = data.FrameData.Velocity;  
    FBLocation = data.FBLocation;
    % Sort FBLocation and apply the sorting to velocities
    [~, sortIdx] = sort(FBLocation);
    FBLocation = FBLocation(sortIdx);
    velocities = velocities(sortIdx, :, :);

    % Take the absolute value of velocities
    velocities = abs(velocities);
    % Extract forces in the x-dimension
    velocitiesX = velocities(:, :, 1);  % Assuming x-dimension is the first dimension
    % Calculate mean force in the x-dimension for every 80 trials
     col = 1;
     for i = 1:80:400
         avgvelocitiesX(fileIndex, col) = mean(velocitiesX(i:i + 79, :), 'all');
         col = col + 1;
     end
     
end

% Calculate mean and standard deviation across columns of velocities
MeanvelocitiesX = mean(avgvelocitiesX); 
STDvelocitiesX = std(avgvelocitiesX);

% Create figure for shaded area
figure;
X4 = [PendLengths fliplr(PendLengths)];
Y4 = [MeanvelocitiesX+STDvelocitiesX fliplr(MeanvelocitiesX-STDvelocitiesX)];
fill(X4, Y4, [0.4940 0.1840 0.5560], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
hold on;
plot(PendLengths, MeanvelocitiesX, '-', 'Color', [0.4940 0.1840 0.5560], 'MarkerSize', 8);
xlabel('Pendulum Lengths');
ylabel('Mean Velocities');
title('Mean Velocities & Standard Deviation ~ Pendulum Length');
grid on;

%Lenght vs handle position
avghandlepositionsX = [];
% Loop through each MAT file
for fileIndex = 1:length(matFiles)
    filePath = fullfile(folderPath, matFiles(fileIndex).name);

    load(filePath, 'data');

    % Extract velocities from the loaded data
    handlepositions = data.FrameData.Position;  
    FBLocation = data.FBLocation;
    % Sort FBLocation and apply the sorting to velocities
    [~, sortIdx] = sort(FBLocation);
    FBLocation = FBLocation(sortIdx);
    handlepositions = handlepositions(sortIdx, :, :);

    % Take the absolute value of velocities
    handlepositions = abs(handlepositions);
    % Extract forces in the x-dimension
    handlepositionsX = handlepositions(:, :, 1);  % Assuming x-dimension is the first dimension
    % Calculate mean force in the x-dimension for every 80 trials
     col = 1;
     for i = 1:80:400
         avghandlepositionsX(fileIndex, col) = mean(handlepositionsX(i:i + 79, :), 'all');
         col = col + 1;
     end
     
end

% Calculate mean and standard deviation across columns of scoreMat
MeanhandlepositionsX = mean(avghandlepositionsX); 
STDhandlepositionsX = std(avghandlepositionsX);

% Create figure for shaded area
figure;
plot(PendLengths, MeanhandlepositionsX, '-', 'Color', [0.4940 0.1840 0.5560], 'MarkerSize', 8);
hold on;
% Create shaded area for standard deviation
X7 = [PendLengths, fliplr(PendLengths)];
Y_upper = MeanhandlepositionsX + STDhandlepositionsX;
Y_lower = MeanhandlepositionsX - STDhandlepositionsX;
fill(X7, [Y_upper, fliplr(Y_lower)], [0.4940 0.1840 0.5560], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

xlabel('Pendulum Lengths');
ylabel('Handle Position');
title('Handle Position & Standard Deviation ~ Pendulum Length');
grid on;

hold off;

