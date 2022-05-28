clearvars -except output_data n_tests field_length_tests noise_stdev_tests
close all
% n_tests=10:10:50;
% field_length_tests = 1:0.5:5;
% noise_stdev_tests = 0:0.5:2.5;
%[n noise_stdev field_length success_rate]

%%process for graphs of field_length and noise stdev vs max error 
for i=1:length(n_tests) %diff plots
    figure(i);
    hold on;
    num_trials_per_n_trial = length(field_length_tests)*length(noise_stdev_tests);
    start_index = num_trials_per_n_trial *(i-1)+1;
    num_trials_per_noise_trial = length(field_length_tests);
    for j=1:length(noise_stdev_tests) %diff traces
        start_row = start_index+((j-1)*num_trials_per_noise_trial);
        plot( output_data( start_row:(start_row+num_trials_per_noise_trial-1),3), ...
                output_data( start_row:(start_row+num_trials_per_noise_trial-1),4));
    end
%     plot3( output_data(start_index:(start_index+num_trials_per_trial-1),2), ...
%             output_data(start_index:(start_index+num_trials_per_trial-1),3), ...
%             output_data(start_index:(start_index+num_trials_per_trial-1),4), 'o' );
    
    xlabel("field length");ylabel("max error");
    title(strcat(string(n_tests(i))," nodes"));
% % legend({'0.0m','0.2m','0.4m','0.6m','0.8m','1.0m'},'Location','southeast');
    legend({'0.0m','1.0m','2.0m','3.0m'},'Location','northeast');
    axis([min(field_length_tests) max(field_length_tests) 0 5])
end

% %%%restructuring for node density
% %[noise_stdev node_density success_rate]
% restructured_data = zeros(length(output_data),3);
% for i=1:length(output_data)
%     restructured_data(i,1) = output_data(i,2);
%     restructured_data(i,2) = output_data(i,1)/(output_data(i,3))^2;
%     restructured_data(i,3) = output_data(i,4);
% end
% restructured_data = sortrows(restructured_data);
% num_trials_per_noise_trial = length(output_data)/length(noise_stdev_tests);
% figure(1)
% hold on
% for j=1:length(noise_stdev_tests)
%     start_row = 1+(j-1)*num_trials_per_noise_trial;
%     plot( restructured_data(start_row:(start_row+num_trials_per_noise_trial-1),2), ...
%           restructured_data(start_row:(start_row+num_trials_per_noise_trial-1),3));
% end
% 
% xlabel("Node density");ylabel("Success Rate");
% % legend({'0.0m','0.2m','0.4m','0.6m','0.8m','1.0m'},'Location','southeast');
% legend({'0.0m','0.25m','0.50m','0.75m','1.0m'},'Location','southeast');
% axis([0 1 0 1])
    