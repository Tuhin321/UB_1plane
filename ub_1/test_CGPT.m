% Generate synthetic data for a 44-degree of freedom system with equations of motion:
% F = -K * X - C * V
% Where F is a vector of forces, K is a stiffness matrix, X is a vector of displacements,
% C is a damping matrix, and V is a vector of velocities.

% Generate synthetic data (simplified example)
rng(0); % Set random seed for reproducibility
num_samples = 100;
num_dofs = 10;

% Generate random stiffness and damping matrices (simplified)
K_true = rand(num_dofs, num_dofs); % True stiffness matrix
C_true = rand(num_dofs, num_dofs); % True damping matrix

X_train = rand(num_samples, num_dofs); % Input features (displacements)

% Calculate ground truth forces using simplified equations
Y_train = -(K_true * X_train' + C_true * randn(num_dofs, num_samples));

% Define a custom feedforward neural network model with ReLU activation
hiddenLayerSizes = [128, 64, 32]; % Three hidden layers with 128, 64, and 32 neurons
model = patternnet(hiddenLayerSizes);

% Configure training options
options = trainingOptions('adam', ...
    'MaxEpochs', 500, ...
    'Verbose', true);

% Transpose the input and target data
X_train = X_train'; % Transpose the input data
Y_train = Y_train;  % Transpose is not needed for the target data

% Train the model
[model, ~] = train(model, X_train, Y_train); % Use transposed data

% Predict forces for new data (adjust input size accordingly)
X_new = rand(1, num_dofs); % Example input data for a 44-DOF system
V_new = rand(1, num_dofs); % Example input data for a 44-DOF system
predicted_forces = model(X_new'); % Use the trained model

% Calculate actual forces based on the equations of motion
actual_forces = -K_true * X_new' - C_true * V_new';

% Display predicted and actual forces
fprintf('Predicted Forces:\n');
disp(predicted_forces);
fprintf('Actual Forces:\n');
disp(actual_forces);

% Plot predicted forces vs. actual forces
figure;
plot(predicted_forces, actual_forces, 'o');
hold on
plot(actual_forces, actual_forces, '-b');
xlabel('Predicted Forces');
ylabel('Actual Forces');
title('Predicted Forces vs. Actual Forces');
grid on;

