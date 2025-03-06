% Clear workspace and figures
clear;

% -----------------------------
% 1. Define Grid Centers
% -----------------------------

% Define x and y grid centers
x_centers = -1.0 : 0.1 : 1.0;      % x centers from -0.7 to 0.7 with step 0.05
y_centers = -0.1 : 0.005 : -0.005;  % y centers from -0.1 to -0.005 with step 0.005

% Initialize an empty data matrix with NaN (to represent missing values)
data_grid = NaN(length(y_centers), length(x_centers));

% -----------------------------
% 2. Define Known Data Points
% -----------------------------

% Example data points with known x and y indices
% Each row: [x_index, y_index, value]
known_data = [
    -0.5, -0.005, -14.14; 
    -0.6, -0.005, -8.785;   
    0.5, -0.005, -14.14;  
    0.6, -0.005, -8.785;

    -0.5, -0.01, -15.81;
    -0.6, -0.01, -9.38;
    0.5, -0.01, -15.81;
    0.6, -0.01, -9.38;

    -0.5, -0.02, -20.755;
    -0.6, -0.02, -10.93;
    0.5, -0.02, -20.755;
    0.6, -0.02, -10.93;

    -0.7, -0.05, -10.325;
    0.7, -0.05, -10.325;
];

for I0 = (-0.4:0.1:0.4)
    known_data = [known_data; I0, -0.005, -24];
    known_data = [known_data; I0, -0.01, -24];
    known_data = [known_data; I0, -0.02, -24];
end

for I0 = (-0.5:0.2:0.5)
    known_data = [known_data; I0, -0.05, -24];
end

% Populate the data_grid with the known values
for i = 1:size(known_data, 1)
    x_idx = known_data(i, 1); % x index in x_centers
    y_idx = known_data(i, 2); % y index in y_centers
    value = known_data(i, 3); % data value
    
    % Set the data value in the grid
    data_grid(round((y_idx+0.1)*200)+1, round((x_idx+1)*10)+1) = value;
end

% -----------------------------
% 3. Plotting the Data
% -----------------------------

% Create a new figure window with equal width and height
%figure('Position', [100, 100, 600, 600]);  % Adjust figure size as needed

% Use imagesc to display the data with NaN values left as white
imagesc(x_centers, y_centers, data_grid);

% Set colormap (choose any MATLAB colormap you prefer)
colormap(jet); % Other options: 'parula', 'hot', 'cool', etc.

% Add colorbar to indicate the scale of values
colorbar;

% Set NaN values to white by customizing the colormap
colormap_with_nan = [1, 1, 1; colormap];  % Add white at the start of the colormap
colormap(colormap_with_nan);
clim([-25,-6]);

% -----------------------------
% 4. Configuring the Axes
% -----------------------------

% Get current axes handle
ax = gca;

% Set X and Y ticks to correspond to grid centers
ax.XTick = x_centers;
ax.YTick = y_centers;

% Enable grid lines for clarity
grid on;

% Customize grid line properties
ax.GridColor = [0 0 0]; % Black grid lines
ax.GridAlpha = 0.5;     % Semi-transparent grid lines
ax.LineWidth = 1;       % Thickness of grid lines

% Set aspect ratio to ensure square cells
axis square;

% Reverse the Y-axis to match typical matrix display format
set(ax, 'YDir', 'normal');

% Set axis limits as specified
xlim([-1.1, 1.1]);
ylim([-0.12, 0]);

% -----------------------------
% 5. Adding Labels and Title
% -----------------------------

xlabel('I0');
ylabel('J0');
title('K0=-0.3');

% -----------------------------
% 6. Save the Figure (Optional)
% -----------------------------

% Save the figure as a PNG image
% Uncomment the following line to save
% saveas(gcf, 'Discrete2DGrid_MissingData.png');
