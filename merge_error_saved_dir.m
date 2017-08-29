function merge_error_saved_dir

% perform solo.merge_error_saved in the entire data in the current
% directory. if the one has already '*x.mat', then it skips that one
% because it's already done

% Assume at most 1 session in a day.


full_list = ls;
full_list = full_list(3:end,:);
data_list = []; % files that are not autosaved
for i = 1 : size(full_list,1)
    if contains(full_list(i,:),'autosave')
    elseif contains(full_list(i,:),'data_@pole_')
        data_list = [data_list; full_list(i,:)];
    end
end

for i = 1 : size(data_list,1)
    date = data_list(i,end-10:end-5);
    merge_error_saved(date);
end