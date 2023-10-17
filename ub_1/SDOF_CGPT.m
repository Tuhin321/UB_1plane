clear all
clc

%% Example System Parameters
mass = 1.0;        % Mass of the rotor (kg)
stiffness = 100.0; % Stiffness of the bearings (N/m)
damping = 2.0;     % Damping coefficient of the bearings (Ns/m)
external_force = 5.0; % External force applied to the rotor (N)
displacement_amplitude = 0.1; % Amplitude of the rotor displacement (m)
velocity_amplitude = 0.2;    % Amplitude of the rotor velocity (m/s)

%% Simulate the system response
num_samples = 100; % Number of data samples
time = linspace(0, 1, num_samples); % Time vector (1 second simulation)
displacements = displacement_amplitude * sin(2 * pi * time); % Sinusoidal displacement
velocities = velocity_amplitude * cos(2 * pi * time);      % Cosine velocity

% Calculate forces based on the equation of motion for SDOF system
forces = -(stiffness * displacements + damping * velocities) + external_force;

% Generate training data
X_train = [displacements; velocities]; % Input features (displacements and velocities)
Y_train = forces; % Output labels (forces)

%% NEURAL NETWORK
% Define the percentage for validation data (e.g., 20%)
validation_percentage = 20;

% Split the data into training and validation sets
num_samples = size(X_train, 2); % Total number of samples
num_validation = round(validation_percentage / 100 * num_samples); % Number of validation samples
num_training = num_samples - num_validation; % Number of training samples

% Shuffle the data (if needed)
shuffled_indices = randperm(num_samples);
X_shuffled = X_train(:, shuffled_indices);
Y_shuffled = Y_train(shuffled_indices);

% Split the data
X_training = X_shuffled(:, 1:num_training);
Y_training = Y_shuffled(1:num_training);

X_validation = X_shuffled(:, (num_training + 1):end);
Y_validation = Y_shuffled((num_training + 1):end);

% Define a custom feedforward neural network model with ReLU activation
hiddenLayerSizes = [128, 64, 32]; % Three hidden layers with 128, 64, and 32 neurons
model = patternnet(hiddenLayerSizes);

% Configure training options with early stopping
options = trainingOptions('adam', ...
    'MaxEpochs', 1000, ... % Increase maximum epochs
    'Verbose', true, ...
    'ValidationData', {X_validation, Y_validation}, ... % Add validation data
    'Plots', 'training-progress', ... % Plot training progress
    'ValidationPatience', 20, ... % Number of epochs with no improvement before stopping
    'MiniBatchSize', 16); % Adjust mini-batch size as needed

% Train the model
[model1, ~] = train(model, X_train, Y_train); % Use transposed data


% Simulate the system response for a new input
new_displacement = displacement_amplitude * sin(2 * pi * 1.5); % Example new displacement
new_velocity = velocity_amplitude * cos(2 * pi * 1.5);      % Example new velocity
predicted_force = model1([new_displacement; new_velocity]);   % Use the trained model

% Calculate actual force based on the equation of motion
actual_force = -(stiffness * new_displacement + damping * new_velocity) + external_force;

% Display predicted and actual forces
fprintf('Predicted Force: %.2f N\n', predicted_force);
fprintf('Actual Force: %.2f N\n', actual_force);

%% Plot predicted force vs. actual force
figure;
plot(predicted_force, actual_force, 'o');
hold on;
plot(actual_force, actual_force, '-b');
xlabel('Predicted Force (N)');
ylabel('Actual Force (N)');
title('Predicted Force vs. Actual Force');
grid on;
legend('Predicted vs. Actual', 'Perfect Fit');
