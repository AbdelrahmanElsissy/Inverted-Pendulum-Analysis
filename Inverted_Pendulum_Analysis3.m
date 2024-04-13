%% homework part 2

% List for filepaths
filePaths = {
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S1_D1_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S1_D2_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S2_D1_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S2_D2_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S3_D1_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S3_D2_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S4_D1_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S4_D2_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S5_D1_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S5_D2_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S6_D1_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S6_D2_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S7_D1_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S7_D2_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S8_D1_data.mat';
'C:\Users\PC\Desktop\Human Robotics\Inverted Pendulum\ExtractedData\S8_D2_data.mat';
};
% Schleife durch jede Datei
for i = 1:length(filePaths)
    % Lade die Daten aus der Datei
    load(filePaths{i}, 'data');

% Annahme: data.FrameData.VisState ist eine Matrix mit 400 Zeilen
[rowCount, ~] = size(data.FrameData.VisState);

% Initialisiere das Feld Pert_Start_Index mit Nullen
data.Pert_Start_Index = zeros(rowCount, 1);

for row = 1:rowCount
    % Suche nach dem Index, an dem data.FrameData.VisState(row, :) >= 1 ist
    index = find(data.FrameData.VisState(row, :) >= 1, 1, 'first');

     % keine Perturbation => setzte für den index 1 ein
    if isempty(index)
    index=1;
    end

    % Aktualisiere das Feld Pert_Start_Index mit dem gefundenen Index
    data.Pert_Start_Index(row, 1) = index;
end


    % Speichere die aktualisierten Daten zurück in die Datei
    save(filePaths{i}, 'data');
end
%% %% cut für schnellere Laufzeit (nicht kompletten Code laufen lassen)
% Zeitpunkte auf Null gleichsetzen


% Schleife durch jede Datei und indiziere data.FrameData.TrialTimeZero
for i = 1:length(filePaths)
    % Lade die Daten aus der Datei
    load(filePaths{i}, 'data');

    % Initialisiere TrialTimeZero als Matrix mit der gleichen Größe wie TrialTime
    %data.FrameData.TrialTimeZero = zeros(size(data.FrameData.TrialTime));

        %für jede Zeile - Loop sprint immer eine Zeile weiter - Zeilindex ist j
    %for j=1:1:400
    %data.FrameData.TrialTimeZero(j,:)= data.FrameData.TrialTime(i,:)- data.FrameData.TrialTime(i,data.Pert_Start_Index(i,1));
    %end
    
    % Anzahl der Trials
    numTrials = size(data.FrameData.TrialTime, 1);
    
    % Initialisiere TrialTimeZero als Zellarray, jedes Element ist ein Zeitvektor
    data.FrameData.TrialTimeZero = cell(numTrials, 1);
    
    % TrialTime wird zu TrialTimeZero, wo der Zeitpunkt 0 gleich dem Start
    % der Purturbation ist
    for trial = 1:numTrials
        data.FrameData.TrialTimeZero{trial} = data.FrameData.TrialTime(trial,1:data.FrameData.Frames(trial,1))- data.FrameData.TrialTime(trial,data.Pert_Start_Index(trial,1));
    end


    % Speichere die aktualisierten Daten zurück in die Datei
    save(filePaths{i}, 'data');
end
%% Cut 
% Interpolation

% Schleife durch jede Datei und speichere die interpoliere die Kraft in data.FrameData.ForceInterpolated
for i = 1:length(filePaths)
    % Lade die Daten aus der Datei
    load(filePaths{i}, 'data');

    %für jede Zeile - Loop sprint immer eine Zeile weiter - Zeilindex ist trial
    for trial=1:1:400
           Frame=data.FrameData.Frames(trial,1);
           forceDataX = data.FrameData.Forces(trial,1:Frame,1);
           data.FrameData.ForceInterpolated2(trial,:) = interp1(data.FrameData.TrialTimeZero{trial}(1:Frame),forceDataX,0.180:0.001:0.230);
           %we start at minus 100 and therefor take force=0 at 100msec =>
           %200
           data.FrameData.ForceInterpolated2(trial,:)= data.FrameData.ForceInterpolated2(trial,:)-data.FrameData.ForceInterpolated2(trial,230);
    end       

    % Speichere die aktualisierten Daten zurück in die Datei
    save(filePaths{i}, 'data');
end
%% %% Cut 
%Matrix für Grafik -3 0 +3


%Matrix for Conditions
%condition +3      PertCondition = 3 oder 6
%condtiton -3      PertCondition = 4 oder 5
% 3 matrixes with different trials and th same condition

%matrix mit condition +3
% Initialisiere die Matrix matrix_plus3

PendLengthVec=[0.75,1,1.5,2,4];
for i=1:5
    for j=1:3
        ForceMatrix2{j,i}=[];
    end
end
% Schleife durch jede Datei
for i = 1:length(filePaths)
    % Lade die Daten aus der Datei
    load(filePaths{i}, 'data');
    
    for trial=1:400
        Flag=1;
        PendLength=data.FBLocation(trial)*2;
        Col=find(PendLengthVec==PendLength);
        if (data.PertCondNum(trial,1) == 3 || data.PertCondNum(trial,1) == 6)
            Row=1;
        elseif (data.PertCondNum(trial,1) == 7 || data.PertCondNum(trial,1) == 8)
            Row=2;    
        elseif (data.PertCondNum(trial,1) == 4 || data.PertCondNum(trial,1) == 5)
            Row=3;
        else
           Flag=0;
        end
        if Flag==1
            ForceMatrix{Row,Col}=[ForceMatrix{Row,Col}(:,:);data.FrameData.ForceInterpolated2(trial,:)];
        end
    end

end
%% %% %%
% Create figure
figure;

% Define custom colors based on pendulum length with greater intensity differences
colors = [
    0.7, 0.1, 0.1;   % Dark Red (0.75 m)
    0.9, 0.2, 0.1;   % Red-Orange (1.5 m)
    1, 0.5, 0;       % Orange-Yellow (2.25 m)
    1, 0.8, 0;       % Yellow (3 m)
    1, 0.95, 0.6     % Light Yellow (4 m)
];

% Loop over pendulum lengths
for pendLength = 1:5
    % Initialize matrix to store forces for both conditions
    combinedForces2 = [];
    
    % Loop over selected conditions (1 and 3)
    for condition = [1, 3]
        % Extract forces for the current condition and pendulum length
        forcesMatrix2 = ForceMatrix2{condition, pendLength};
        
        % For condition 1, change the values to their negation
        if condition == 1
            forcesMatrix2 = -forcesMatrix2;
        end
        
        % Concatenate forces for both conditions
        combinedForces2 = [combinedForces2; forcesMatrix2];
    end
    
    % Calculate mean force for combined conditions
    meanForce2 = mean(combinedForces2);
    
    % Calculate standard error for each pendulum length
    stdError2 = std(combinedForces2) / sqrt(size(combinedForces2, 1));
    
    % Plot mean force for each pendulum length with color based on custom colors
    plot(-0.1:0.001:0.3, meanForce, 'LineWidth', 1.5, 'Color', colors(pendLength, :));
    hold on;
    
    % Plot shaded area for standard error
    X = [0.180:0.001:0.230, fliplr(-0.1:0.001:0.3)];
    Y = [meanForce + stdError, fliplr(meanForce - stdError)];
    fill(X, Y, colors(pendLength, :), 'EdgeColor', 'none', 'FaceAlpha', 0.3);


    % Plot shaded area from 0.180 to 0.230 in light grey
    if pendLength == 2 % Only for condition 0
        X_lightgrey = [0.180, 0.230, 0.230, 0.180];
        Y_lightgrey = [-0.4, -0.4, 1.6, 1.6];
        fill(X_lightgrey, Y_lightgrey, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
    end
    % Plot shaded area from 0.230 to 0.300 in grey
    if pendLength == 2 % Only for condition 0
        X_grey = [0.230, 0.300, 0.300, 0.230];
        Y_grey = [-0.4, -0.4, 1.6, 1.6];
        fill(X_grey, Y_grey, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    end
end

% Customize plot
xlabel('Time [ms]');
ylabel('Hand Lateral Force [N]');
title('Mean Force ~ Time for Different Pendulum Lengths (Conditions -3 and 3)');
grid off;

% Add dotted vertical line at time = 0
xline(0, '--k', 'LineWidth', 1.2);

% Add a horizontal dashed line at y = 0
yline(0, 'k--', 'LineWidth', 1.2);
   

% Create a custom legend with the specified colors
legend('0.75 m', '1.5 m','2.25 m','3 m', '3 m', '4 m', 'Location', 'Best');

