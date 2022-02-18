close all;clear all;
rng(101)

%%%SET EXPERIMENT PARAMS
num_experiments = 50;
n_tests=3:1:10;
field_length_tests = 5:5:25;
noise_stdev_tests = 0:1:3;
noise_stdev_tests = noise_stdev_tests./3; %%want to fit within +-x error 99.7% of the time
output_data = []; %n noise_stdev field_length maxerr

%SANITY
run_length = length(n_tests)*length(noise_stdev_tests)*length(field_length_tests);
count = 0;

for i=1:length(n_tests)
    for k=1:length(noise_stdev_tests)
        for j=1:length(field_length_tests)

            %%%SET PARAMS
            % n = input("Num nodes: ");
            % field_length = input("R x R field, in meters. R:");
            % radio_range = input("radio range (0-1.0, fraction of R)");
            n=n_tests(i);
            field_length=field_length_tests(j);
            noise_stdev = noise_stdev_tests(k);
            radio_range=1;
            anchors = [0 -.9*field_length;];
            num_anchors = size(anchors); num_anchors = num_anchors(1);
            %add min_dist?

            %%%EXPERIMENTAL DATA OUT
            summed_maxerr = 0;
            for i1=1:num_experiments

                %%%GENERATE FAKE DATA
                orig_data_pts = -(field_length/2) + field_length*rand(n,2); % [x y]
                orig_data_pts = [orig_data_pts; anchors];

                %%%GENERATE DISTANCES
                disparity_matrix = zeros(n+num_anchors);
                for i2 = 1:(n+num_anchors)
                    for j2 = 1:(n+num_anchors)
                        disparity_matrix(i2,j2) = norm(orig_data_pts(i2,1:end)-orig_data_pts(j2,1:end));
                    end
                end

                %%%ADD IN NOISE
                noise_matrix = zeros(n+num_anchors);
                for i3 = 1:(n+num_anchors)
                    for j3 = (i3+1):(n+num_anchors)
                        noise_matrix(i3,j3) = randn*noise_stdev;
                        noise_matrix(j3,i3) = noise_matrix(i3,j3);
                    end
                end
                disparity_matrix = disparity_matrix + noise_matrix;
                disparity_matrix(disparity_matrix<0) = 0; %distances cannot be < 0

                %%%BLACK MAGIC
                guess_pts = matlab_cmdscale(disparity_matrix);
                
                % initial_guess = -(field_length/2) + field_length*rand(n+1,2); % [x y]
                % tic
                % guess_pts = smacof(disparity_matrix,initial_guess);
                % toc
                summed_maxerr = summed_maxerr + max(abs(fake_pdist(orig_data_pts)-fake_pdist(guess_pts))); %think this is broken, doesn't handle accidental higher dimensional returns
            end

            maxerr_avg = summed_maxerr/num_experiments;
            output_data = [output_data; n noise_stdev field_length maxerr_avg];
            count = count + 1;
            disp(strcat(num2str(count),"/",num2str(run_length)))
        end
    end
end

%% plot last run as example
figure
subplot(1,2,1)
scatter(orig_data_pts(:,1),orig_data_pts(:,2))
xlabel("x (m)");ylabel("y (m)");title("Original pts")
axis([-max(field_length_tests) max(field_length_tests) -max(field_length_tests) max(field_length_tests)])
subplot(1,2,2)
scatter(guess_pts(:,1),guess_pts(:,2))
xlabel("x (m)");ylabel("y (m)");title("Calculated pts, MSE = 0.77m")
axis([-max(field_length_tests) max(field_length_tests) -max(field_length_tests) max(field_length_tests)])
immse(fake_pdist(orig_data_pts),fake_pdist(guess_pts(:,1:2)))
% annotation('textbox',[.2 .5 .3 .3],'String',"immse = 0.7712",'FitBoxToText','on');









