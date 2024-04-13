%% Load Data
clear all; 
close all; 
clc; 
% Initialize the data structure
allData = struct('subject', cell(1, 8), 'day1', cell(1, 8), 'day2', cell(1, 8));

% Define the folder path
folderPath = 'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\';

% Load data of 8 subjects
for i = 1:8
    % Construct the file paths for each day
    day1_file = fullfile(folderPath, sprintf('S%d_D1_data.mat', i));
    day2_file = fullfile(folderPath, sprintf('S%d_D2_data.mat', i));

    % Check if the files exist before loading
    if exist(day1_file, 'file') == 2 && exist(day2_file, 'file') == 2
        % Load data for each day
        day1_data = load(day1_file);
        day2_data = load(day2_file);

        % Save data in structure
        allData(i).subject = i;
        allData(i).day1 = day1_data;
        allData(i).day2 = day2_data;
    else
        disp(['Error: Files not found for subject ', num2str(i)]);
    end
end


%% Find indices where visual perturbation starts 

% Initialize vis_pert_start matrix
vis_pert_start = zeros(400, 8);

for subject = 1:8
    for i = 1:400
        % Find the index where visual perturbation starts
        vis_start_index = find(allData(subject).day1.data.FrameData.VisState(i,:) == 1, 1, 'first');
      
        % Check if the result is empty (no visual perturbation)
        if ~isempty(vis_start_index)
            vis_pert_start(i, subject) = vis_start_index;


        else
            % Assign a default value (e.g., 0) when no perturbation is found
            vis_pert_start(i, subject) = 0;
        end
    end 
end


%% Set to 0 

for subject = 1:8 
    allData(subject).day1.data.FrameData.TrialTimeZero = cell(400, 1);
    %allData(subject).day2.data.FrameData.TrialTimeZero = cell(400, 1);
end 

for subject = 1:8 
    for trials = 1:400
        pert_start_index = vis_pert_start(trials, subject);
        
        % Check if perturbation start index is not 0
        if pert_start_index ~= 0
            allData(subject).day1.data.FrameData.TrialTimeZero{trials} = ... 
                allData(subject).day1.data.FrameData.TrialTime(trials, 1:allData(subject).day1.data.FrameData.Frames(trials, 1)) - ...
                allData(subject).day1.data.FrameData.TrialTime(trials, pert_start_index);
        else
            % fprintf('Subject: %d, Trial: %d excluded (Perturbation Start Index is 0)\n', subject, trials);
            % Optionally mark or exclude these trials
            % For example: allData(subject).day1.data.FrameData.TrialTimeZero{trials} = NaN;
            continue; % Skip to the next iteration
        end
    end
end


%% Interpolation 

for subject = 1:8 
    for trials = 1:400 
        % Check if TrialTimeZero is empty for the current trial
        if isempty(allData(subject).day1.data.FrameData.TrialTimeZero{trials})
            % Skip this trial
            %fprintf('Subject: %d, Trial: %d skipped (TrialTimeZero is empty)\n', subject, trials);
            continue;
        end
        
        frames = allData(subject).day1.data.FrameData.Frames(trials, 1);
        force_X = allData(subject).day1.data.FrameData.Forces(trials, 1:frames, 1);

        %fprintf('Subject: %d, Trial: %d, frames: %d, size(TrialTimeZero): %s\n', subject, trials, frames, mat2str(size(allData(subject).day1.data.FrameData.TrialTimeZero{trials})));

        % Interpolate only if TrialTimeZero is not empty
        allData(subject).day1.data.FrameData.Force_X_Interp{trials} = ...
            interp1(allData(subject).day1.data.FrameData.TrialTimeZero{trials}(1:frames), force_X, -0.1:0.001:0.3); 

        % Additional processing (subtracting value at 230)
        allData(subject).day1.data.FrameData.Force_X_Interp{trials} = ...
            allData(subject).day1.data.FrameData.Force_X_Interp{trials} - ...
            allData(subject).day1.data.FrameData.Force_X_Interp{trials}(230);
    end
end



%%

% Matrix for Conditions
% condition +3      PertCondition = 3 or 6
% condition -3      PertCondition = 4 or 5
% 3 matrices with different trials and the same condition

% Initialize the ForceMatrix
PendLengthVec = [0.75, 1, 1.5, 2, 4];

for subject = 1:8
    for i = 1:5
        for j = 1:3
            allData(subject).day1.data.FrameData.ForceMatrix{j, i} = [];
        end
    end
end

for subject = 1:8
    for trials = 1:400
        Flag = 1;
        PendLength = allData(subject).day1.data.FBLocation(trials) * 2;
        Col = find(PendLengthVec == PendLength);
        
        if (allData(subject).day1.data.PertCondition(trials, 1) == 3 || allData(subject).day1.data.PertCondition(trials, 1) == 6)
            Row = 1;
        elseif (allData(subject).day1.data.PertCondition(trials, 1) == 7 || allData(subject).day1.data.PertCondition(trials, 1) == 8)
            Row = 2;
        elseif (allData(subject).day1.data.PertCondition(trials, 1) == 4 || allData(subject).day1.data.PertCondition(trials, 1) == 5)
            Row = 3;
        else
            Flag = 0;
        end
        
        if Flag == 1
            allData(subject).day1.data.FrameData.ForceMatrix{Row, Col} = [allData(subject).day1.data.FrameData.ForceMatrix{Row, Col}; allData(subject).day1.data.FrameData.Force_X_Interp{trials}];
        end
    end
end

%% Find indices where visual perturbation starts 

% Initialize vis_pert_start matrix
vis_pert_start = zeros(400, 8);

for subject = 1:8
    for i = 1:400
        % Find the index where visual perturbation starts
        vis_start_index = find(allData(subject).day2.data.FrameData.VisState(i,:) == 1, 1, 'first');
      
        % Check if the result is empty (no visual perturbation)
        if ~isempty(vis_start_index)
            vis_pert_start(i, subject) = vis_start_index;


        else
            % Assign a default value (e.g., 0) when no perturbation is found
            vis_pert_start(i, subject) = 0;
        end
    end 
end


%% Set to 0 

for subject = 1:8 
    %allData(subject).day1.data.FrameData.TrialTimeZero = cell(400, 1);
    allData(subject).day2.data.FrameData.TrialTimeZero = cell(400, 1);
end 

for subject = 1:8 
    for trials = 1:400
        pert_start_index = vis_pert_start(trials, subject);
        
        % Check if perturbation start index is not 0
        if pert_start_index ~= 0

                allData(subject).day2.data.FrameData.TrialTimeZero{trials} = ... 
                allData(subject).day2.data.FrameData.TrialTime(trials, 1:allData(subject).day2.data.FrameData.Frames(trials, 1)) - ...
                allData(subject).day2.data.FrameData.TrialTime(trials, pert_start_index);
        else
            % fprintf('Subject: %d, Trial: %d excluded (Perturbation Start Index is 0)\n', subject, trials);
            % Optionally mark or exclude these trials
            % For example: allData(subject).day1.data.FrameData.TrialTimeZero{trials} = NaN;
            continue; % Skip to the next iteration
        end
    end
end


%% Interpolation 

for subject = 1:8 
    for trials = 1:400 
        % Check if TrialTimeZero is empty for the current trial
        if isempty(allData(subject).day2.data.FrameData.TrialTimeZero{trials})
            % Skip this trial
            %fprintf('Subject: %d, Trial: %d skipped (TrialTimeZero is empty)\n', subject, trials);
            continue;
        end
        
        frames = allData(subject).day2.data.FrameData.Frames(trials, 1);
        force_X = allData(subject).day2.data.FrameData.Forces(trials, 1:frames, 1);

        %fprintf('Subject: %d, Trial: %d, frames: %d, size(TrialTimeZero): %s\n', subject, trials, frames, mat2str(size(allData(subject).day1.data.FrameData.TrialTimeZero{trials})));

        % Interpolate only if TrialTimeZero is not empty
        allData(subject).day2.data.FrameData.Force_X_Interp{trials} = ...
            interp1(allData(subject).day2.data.FrameData.TrialTimeZero{trials}(1:frames), force_X, -0.1:0.001:0.3); 

        % Additional processing (subtracting value at 230)
        allData(subject).day2.data.FrameData.Force_X_Interp{trials} = ...
            allData(subject).day2.data.FrameData.Force_X_Interp{trials} - ...
            allData(subject).day2.data.FrameData.Force_X_Interp{trials}(260);
    end
end



%%

% Matrix for Conditions
% condition +3      PertCondition = 3 or 6
% condition -3      PertCondition = 4 or 5
% 3 matrices with different trials and the same condition

% Initialize the ForceMatrix
PendLengthVec = [0.75, 1, 1.5, 2, 4];

for subject = 1:8
    for i = 1:5
        for j = 1:3
            allData(subject).day2.data.FrameData.ForceMatrix{j, i} = [];
        end
    end
end

for subject = 1:8
    for trials = 1:400
        Flag = 1;
        PendLength = allData(subject).day2.data.FBLocation(trials) * 2;
        Col = find(PendLengthVec == PendLength);
        
        if (allData(subject).day2.data.PertCondition(trials, 1) == 3 || allData(subject).day2.data.PertCondition(trials, 1) == 6)
            Row = 1;
        elseif (allData(subject).day2.data.PertCondition(trials, 1) == 7 || allData(subject).day2.data.PertCondition(trials, 1) == 8)
            Row = 2;
        elseif (allData(subject).day2.data.PertCondition(trials, 1) == 4 || allData(subject).day2.data.PertCondition(trials, 1) == 5)
            Row = 3;
        else
            Flag = 0;
        end
        
        if Flag == 1
            allData(subject).day2.data.FrameData.ForceMatrix{Row, Col} = [allData(subject).day2.data.FrameData.ForceMatrix{Row, Col}; allData(subject).day2.data.FrameData.Force_X_Interp{trials}];
        end
    end
end

%%

% Loop through the columns (P1 to P5) for each subject
for subject = 1:8
    for i = 1:5
        % Calculate mean for each condition (x, y, z) and then combine them
        x_d1 = mean(abs(allData(subject).day1.data.FrameData.ForceMatrix{1, i}(:, 280:330)));
        %y_d1 = mean(abs(allData(subject).day1.data.FrameData.ForceMatrix{2, i}(:, 280:330)));
        z_d1 = mean(abs(allData(subject).day1.data.FrameData.ForceMatrix{3, i}(:, 280:330)));

        x_d2 = mean(abs(allData(subject).day2.data.FrameData.ForceMatrix{1, i}(:, 280:330)));
        %y_d2 = mean(abs(allData(subject).day2.data.FrameData.ForceMatrix{2, i}(:, 280:330)));
        z_d2 = mean(abs(allData(subject).day2.data.FrameData.ForceMatrix{3, i}(:, 280:330)));

        % Calculate the final mean for the current column (P1 to P5)
        P_means(subject, i) = mean([mean(x_d1), mean(z_d1), mean(x_d2), mean(z_d2)]);
        P_sd(subject, i) = std([mean(x_d1), mean(z_d1), mean(x_d2), mean(z_d2)]);
    end
end 

P_means_final = mean(P_means);
P_sd_final = std(P_means);

% Define custom RGB values for the colors
custom_colors = [
   0, 0.0667, 1.0000;   % Dark purple (0.75 m)
    0.4667, 0, 1.0000;   % purple (1.5 m)
    0.7176, 0.2745, 1.0000;   % pink (2.25 m)
    1, 0, 1;       % Yellow (3 m)
    1.0000, 0.0745, 0.6510   % Light Yellow (4 m) 
];

% Plot the bar plot with custom colors
figure;

bar(P_means, 'FaceColor', 'flat');

% Loop through each bar to set face color and add error bars
for i = 1:5
    % Set the face color for each bar
    h = bar(i, P_means_final(i), 'FaceColor', custom_colors(i, :), 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
    hold on;

    % Calculate error bars (standard deviation)
    errorBarColor = custom_colors(i, :);
    line([i, i], [P_means_final(i) - P_sd_final(i), P_means_final(i) + P_sd_final(i)], 'Color', errorBarColor, 'LineWidth', 4);

    % Store the handle to the bar for setting legend later
    bars(i) = h;
end


% Customize the plot
xlabel('Pendulum Length [m]');
str = sprintf('Visuomotor Feedback Response \n180-230 ms [N]');
ylabel(str);
xticks(1:5);
yticks(0:0.2:2);
xticklabels({'0.75', '1', '1.5', '2', '4'});
grid off; 
box off; 
set(gca);

%%

% Loop through the columns (P1 to P5) for each subject
for subject = 1:8
    for i = 1:5
        % Calculate mean for each condition (x, y, z) and then combine them
        x_d1 = mean(abs(allData(subject).day1.data.FrameData.ForceMatrix{1, i}(:, 330:401)));
        %y_d1 = mean(abs(allData(subject).day1.data.FrameData.ForceMatrix{2, i}(:, 330:401)));
        z_d1 = mean(abs(allData(subject).day1.data.FrameData.ForceMatrix{3, i}(:, 330:401)));

        x_d2 = mean(abs(allData(subject).day2.data.FrameData.ForceMatrix{1, i}(:, 330:401)));
        %y_d2 = mean(abs(allData(subject).day2.data.FrameData.ForceMatrix{2, i}(:, 330:401)));
        z_d2 = mean(abs(allData(subject).day2.data.FrameData.ForceMatrix{3, i}(:, 330:401)));

        % Calculate the final mean for the current column (P1 to P5)
        P_means(subject, i) = mean([mean(x_d1), mean(z_d1), mean(x_d2), mean(z_d2)]);
        P_sd(subject, i) = std([mean(x_d1), mean(z_d1), mean(x_d2), mean(z_d2)]);
    end
end 

P_means_final = mean(P_means);
P_sd_final = std(P_means);

% Define custom RGB values for the colors
custom_colors = [
   0, 0.0667, 1.0000;   % Dark purple (0.75 m)
    0.4667, 0, 1.0000;   % purple (1.5 m)
    0.7176, 0.2745, 1.0000;   % pink (2.25 m)
    1, 0, 1;       % Yellow (3 m)
    1.0000, 0.0745, 0.6510   % Light Yellow (4 m)  
];
% Plot the bar plot with custom colors

bar(P_means, 'FaceColor', 'flat');

% Loop through each bar to set face color and add error bars
for i = 1:5
    % Set the face color for each bar
    h = bar(i, P_means_final(i), 'FaceColor', custom_colors(i, :), 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
    hold on;

    % Calculate error bars (standard deviation)
    errorBarColor = custom_colors(i, :);
    line([i, i], [P_means_final(i) - P_sd_final(i), P_means_final(i) + P_sd_final(i)], 'Color', errorBarColor, 'LineWidth', 4);

    % Store the handle to the bar for setting legend later
    bars(i) = h;
end


% Customize the plot
% Customize the plot
xlabel('Pendulum Length [m]');
str = sprintf('Visuomotor Feedback Response \n230-300 ms [N]');
ylabel(str);
xticks(1:5);
yticks(0:0.2:2);
xticklabels({'0.75', '1', '1.5', '2', '4'});
grid off; 
box off; 
set(gca);

