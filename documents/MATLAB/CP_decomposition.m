% Generate synthetic data
rng(0); % Set random seed for reproducibility
students = 100;
subjects = 5;
semesters = 4;

% True factors
student_factor = rand(students, 2); % Two latent factors for students
subject_factor = rand(subjects, 2); % Two latent factors for subjects
semester_factor = rand(semesters, 2); % Two latent factors for semesters

% Generate tensor
data = zeros(students, subjects, semesters);
for i = 1:students
    for j = 1:subjects
        for k = 1:semesters
            data(i, j, k) = sum(student_factor(i, :) .* subject_factor(j, :) .* semester_factor(k, :));
        end
    end
end

% Add noise
noise_level = 0.1;
data = data + noise_level * randn(size(data));

% Reshape data tensor into a matrix
reshaped_data = reshape(data, students, []);

% Perform singular value decomposition (SVD) on reshaped data
[U, S, V] = svd(reshaped_data, 'econ');

% Extract factors from U and V matrices
rank = 3;
student_factor_hat = U(:, 1:rank) * sqrt(S(1:rank, 1:rank)); % Adjust for singular values
subject_factor_hat = V(:, 1:rank) * sqrt(S(1:rank, 1:rank)); % Adjust for singular values
semester_factor_hat = randn(semesters, rank); % Random factors for semester (just for demonstration)

% Reshape factors into tensors
student_factor_hat = reshape(student_factor_hat, students, []);
subject_factor_hat = reshape(subject_factor_hat, subjects, []);
semester_factor_hat = reshape(semester_factor_hat, semesters, []);

% Calculate reconstruction
% Initialize core tensor
core_tensor = rand(rank, rank, rank);
reconstructed_data = zeros(size(data));
for i = 1:rank
    for j = 1:rank
        for k = 1:rank
            reconstructed_data = reconstructed_data + ...%core_tensor(i, j, k) * ...
                tensorprod(tensorprod(student_factor_hat(:, i), subject_factor_hat(:, j), 2, 2), semester_factor_hat(:, k), 3, 2);
        end
    end
end
%for i = 1:rank
%    reconstructed_data = reconstructed_data + tensorprod(tensorprod(student_factor_hat(:, 1),subject_factor_hat(:, 1),2,2),semester_factor_hat(:, 1),3,2);
%end

% Calculate reconstruction error
mse = norm(data(:) - reconstructed_data(:))^2 / numel(data);

disp(['Original Data Size: ', num2str(size(data))]);
disp(['CP Decomposition Rank: ', num2str(rank)]);
disp(['Reconstruction Error (MSE): ', num2str(mse)]);
