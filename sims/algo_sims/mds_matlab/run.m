close all;clear all;

%%%SET PARAMS
% n = input("Num nodes: ");
% field_length = input("R x R field, in meters. R:");
% radio_range = input("radio range (0-1.0, fraction of R)");
n=9;
field_length=3;
radio_range=1;
% anchors = [0 -.9*field_length; ...
%            0 0];
anchors = [0 -.9*field_length];
num_anchors = size(anchors); num_anchors = num_anchors(1);
noise_stdev = 0.5;
labels={};
for i = 1:(n+num_anchors)
    labels = [labels; num2str(i)];
end
%add min_dist?

%%%GENERATE FAKE DATA
[temp_x, temp_y] = ndgrid(linspace(-(field_length/2),(field_length/2),sqrt(n)), linspace(-(field_length/2),(field_length/2),sqrt(n)));
square_grid = zeros(n,2);
for i = 1:n
    square_grid(i,1) = temp_x(i);
    square_grid(i,2) = temp_y(i);
end
orig_data_pts = square_grid;
% orig_data_pts = -(field_length/2) + field_length*rand(n,2); % [x y]
orig_data_pts = [anchors; orig_data_pts];

%%%GENERATE DISTANCES
disparity_matrix = zeros(n+num_anchors);
for i = 1:(n+num_anchors)
    for j = 1:(n+num_anchors)
        disparity_matrix(i,j) = norm(orig_data_pts(i,1:end)-orig_data_pts(j,1:end));
    end
end

%%%ADD IN NOISE
noise_matrix = zeros(n+num_anchors);
for i = 1:(n+num_anchors)
    for j = (i+1):(n+num_anchors)
        noise_matrix(i,j) = randn*noise_stdev; %normal distribution
        noise_matrix(j,i) = noise_matrix(i,j);
    end
end
disparity_matrix = disparity_matrix + noise_matrix;
disparity_matrix(disparity_matrix<0) = 0; %distances cannot be < 0

%%%BLACK MAGIC
tic
[guess_pts,eigvals] = matlab_cmdscale(disparity_matrix);
toc

% initial_guess = -(field_length/2) + field_length*rand(n+1,2); % [x y]
% tic
% guess_pts = smacof(disparity_matrix,initial_guess);
% toc

guess_pts = guess_pts(:,1:2); %in case of extradimensionality

%%%DEBUGGING
orig_data_pts
guess_pts
maxerr = max(abs(fake_pdist(orig_data_pts)-fake_pdist(guess_pts)))

%%%VISUALIZE
figure(1);
set(gcf, 'Position',  [2000, 300, 1350, 500])
subplot(1,2,1);
scatter(orig_data_pts(1:end,1),orig_data_pts(1:end,2));
text(orig_data_pts(:,1),orig_data_pts(:,2),labels);
title("Original");
axis([-(field_length) (field_length) -(field_length) (field_length)])
subplot(1,2,2);
scatter(guess_pts(1:end,1),guess_pts(1:end,2));
text(guess_pts(:,1),guess_pts(:,2),labels);
title("MDS first guess");
axis([-(field_length) (field_length) -(field_length) (field_length)])

%%%ROTATION CHECK
cmd_input = "";
altered_pts = guess_pts;
f="f";r="r";s="s";q="q";
while string(cmd_input) ~= 'q'
    cmd_input=input("operation: f, r, s? ");
    if cmd_input == "r" %rotate
        input_angle = input("rotate by how many degrees cw? ");
        theta = deg2rad(input_angle);
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        altered_pts = altered_pts * R;
    elseif cmd_input == "f" %flip
        altered_pts(1:end,1) = -1*altered_pts(1:end,1);
    elseif cmd_input == "s" %shift
        x_input = input("x shift? ");
        y_input = input("y shift? ");
        altered_pts(:,1) = altered_pts(:,1) + x_input;
        altered_pts(:,2) = altered_pts(:,2) + y_input;
    else 
        cmd_input = "q";
    end
    scatter(orig_data_pts(1:end,1),orig_data_pts(1:end,2));
    hold on
    scatter(altered_pts(1:end,1),altered_pts(1:end,2));
    text(altered_pts(:,1),altered_pts(:,2),labels);
    hold off
    title("After operations");
    axis([-(field_length) (field_length) -(field_length) (field_length)])
    sgtitle(strcat("Noise stdev = ",string(noise_stdev)," m"));
end