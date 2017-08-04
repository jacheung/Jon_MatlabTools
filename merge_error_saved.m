function merge_error_saved(date)
% Sometimes, because of bControl error, there is 2 different saved files,
% usually having different alphabet at the end of the file name (e.g.,
% 1111111a.mat, 11111111b.mat, ...)
% This function is designed to merge those saved files into one for further
% analysis. The resulting file will have alphabet 'x' at the end (e.g.,
% 1111111x.mat). Practically, there will be 2-3 files max. 

% This also fixes another error sometimes occuring: n_done_trials mismatch
% Therefore, it is recommended to run this function for every saved files
% at least once.
% SidesSection_previous_sides are always having +1 or +2 indices, which
% makes it erroneous for merging, and this function makes the number of
% components same.
%
% Input argument date defines saved date, shown at the file name.
% 



% After this, previous_sides and previous_dstrs from SidesSection mean current sides and
% dstrs. Actually, previous_pole_distances, ap_positions, and pole_angles
% from MotorsSection are alreay meaning current ones.




% Assume at most 1 session in a day.

% Currently only for 2port_angdist
% 2016/06/22 JK

if isnumeric(date)
    date = num2str(date);
end
    
full_list = ls;
list = []; % files matching to the date information
for i = 1 : length(full_list)
    if strfind(full_list(i,:),date)
        if strfind(full_list(i,:),'autosave')
        else
            list = [list; full_list(i,:)];
        end
    end
end
% the resulting list will be sorted already.
% assume that the saved fields will be same across different files
% need to merge only saved and saved_history structs
savefn = [list(1,1:end-5), 'x.mat'];
for i = 1 : size(list,1)
    if strfind(list(i,:),savefn);
        sprintf('%s is already post-processed',savefn)
        return
    end
end

load(list(1,:))
savedfields = fieldnames(saved);
for idx_fields = 1 : length(savedfields)
    if strfind(savedfields{idx_fields,1}, 'hit_history')
        fieldname_hit_history = savedfields{idx_fields,1};
        n_done_trials = length(getfield(saved,fieldname_hit_history));
        idx_fields = length(savedfields);
    end
end

% compare every hit_history with it's one previous hit_history to see if
% they are saved after previous one is closed (flushed), or if there is
% redundancy in the saved files
num_file = size(list,1);
all_hit_history = cell(1,num_file);
for i = 1 : num_file
    load(list(i,:))
    all_hit_history{1,i} = getfield(saved,fieldname_hit_history);
end
keep_line = [];
for i = 1 : num_file - 1
    if (length(all_hit_history{i+1}) >= length(all_hit_history{i})) && ...
            isequal(all_hit_history{i+1}(1:length(all_hit_history{i})), all_hit_history{i})
    else
        keep_line = [keep_line; i];
    end
end
keep_line = [keep_line;num_file];
list = list(keep_line,:);


load(list(1,:))
saved.SidesSection_previous_sides = saved.SidesSection_previous_sides(2:end);
saved.SidesSection_previous_dstrs = saved.SidesSection_previous_dstrs(2:end);
if length(saved_history.AnalysisSection_NumTrials) ~= n_done_trials
    n_done_trials = length(saved_history.AnalysisSection_NumTrials);
    hit_history = getfield(saved,fieldname_hit_history);
    prompt_str = ['saved.',fieldname_hit_history,' = hit_history(1:n_done_trials);'];
    eval(prompt_str)
    saved.SidesSection_previous_sides = saved.SidesSection_previous_sides(1,1:n_done_trials);
    saved_history.RewardsSection_LastTrialEvents = saved_history.RewardsSection_LastTrialEvents(1:n_done_trials);
    if length(saved.SidesSection_previous_dstrs) > n_done_trials
        saved.SidesSection_previous_dstrs = saved.SidesSection_previous_dstrs(1,1:n_done_trials);
    end
    if isfield(saved,'MotorsSection_previous_pole_distances')
        saved.MotorsSection_previous_pole_distances = saved.MotorsSection_previous_pole_distances(1,1:n_done_trials);
    end
    if isfield(saved,'MotorsSection_previous_pole_ap_positions')
        saved.MotorsSection_previous_pole_ap_positions = saved.MotorsSection_previous_pole_ap_positions(1,1:n_done_trials);
    end
    if isfield(saved,'MotorsSection_previous_pole_angles')
        saved.MotorsSection_previous_pole_angles = saved.MotorsSection_previous_pole_angles(1,1:n_done_trials);
    end
end
% save(list(1,:),'saved', 'saved_history','saved_autoset'); % Don't change
% the original saved file. It's too risky.

result_saved = saved;
result_saved_history = saved_history;

if size(list,1) > 1
    for i = 2 : size(list,1)  
        load(list(i,:))
        savedfields = fieldnames(saved);
        savedhistoryfields = fieldnames(saved_history);
        n_done_trials = length(getfield(saved,fieldname_hit_history));
        
        saved.SidesSection_previous_sides = saved.SidesSection_previous_sides(2:end);
        saved.SidesSection_previous_dstrs = saved.SidesSection_previous_dstrs(2:end);

        if length(saved_history.AnalysisSection_NumTrials) ~= n_done_trials
            n_done_trials = length(saved_history.AnalysisSection_NumTrials);
            hit_history = getfield(saved,fieldname_hit_history);
            prompt_str = ['saved.',fieldname_hit_history,' = hit_history(1:n_done_trials);'];
            eval(prompt_str)
            saved.SidesSection_previous_sides = saved.SidesSection_previous_sides(1:n_done_trials);
            saved_history.RewardsSection_LastTrialEvents = saved_history.RewardsSection_LastTrialEvents(1:n_done_trials);
            if length(saved.SidesSection_previous_dstrs) > n_done_trials
                saved.SidesSection_previous_dstrs = saved.SidesSection_previous_dstrs(1:n_done_trials);
            end
            if isfield(saved,'MotorsSection_previous_pole_distances')
                saved.MotorsSection_previous_pole_distances = saved.MotorsSection_previous_pole_distances(1,1:n_done_trials);
            end
            if isfield(saved,'MotorsSection_previous_pole_ap_positions')
                saved.MotorsSection_previous_pole_ap_positions = saved.MotorsSection_previous_pole_ap_positions(1,1:n_done_trials);
            end
            if isfield(saved,'MotorsSection_previous_pole_angles')
                saved.MotorsSection_previous_pole_angles = saved.MotorsSection_previous_pole_angles(1,1:n_done_trials);
            end
        end

        for idx_fields = 1 : length(savedfields)
            temp = getfield(saved,savedfields{idx_fields,1});
            if length(temp) == n_done_trials;
                if size(temp,1) == n_done_trials;
                    prompt_str = ['result_saved.',savedfields{idx_fields,1},' = [getfield(result_saved,savedfields{idx_fields,1}); getfield(saved, savedfields{idx_fields,1})];'];
                    eval(prompt_str);
                elseif size(temp,2) == n_done_trials;
                    prompt_str = ['result_saved.',savedfields{idx_fields,1},' = [getfield(result_saved,savedfields{idx_fields,1}), getfield(saved, savedfields{idx_fields,1})];'];
                    eval(prompt_str);
                end
            end
        end
        for idx_fields = 1 : length(savedhistoryfields)
            temp = getfield(saved_history,savedhistoryfields{idx_fields,1});
            if length(temp) == n_done_trials;
                if strfind(savedhistoryfields{idx_fields,1},'AnalysisSection_Num')
                    oldnum = getfield(result_saved_history,savedhistoryfields{idx_fields,1});
                    oldnum1 = []; oldnum2 = []; oldnum3 = [];
                    for ii = 1 : size(oldnum,1)
                        oldnum1 = [oldnum1; str2num(oldnum{ii,1}(:,1:4))];
                        oldnum2 = [oldnum2; str2num(oldnum{ii,1}(:,7:10))];
                        oldnum3 = [oldnum3; str2num(oldnum{ii,1}(:,13:16))];
                    end
                    tempnum = getfield(saved_history,savedhistoryfields{idx_fields,1});
                    tempnum1 = []; tempnum2 = []; tempnum3 = [];
                    for ii = 1 : size(tempnum,1)
                        tempnum1 = [tempnum1; str2num(tempnum{ii,1}(:,1:4))];
                        tempnum2 = [tempnum2; str2num(tempnum{ii,1}(:,7:10))];
                        tempnum3 = [tempnum3; str2num(tempnum{ii,1}(:,13:16))];
                    end
                    tempnum1 = tempnum1 + max(oldnum1);
                    tempnum2 = tempnum2 + max(oldnum2);
                    tempnum3 = tempnum3 + max(oldnum3);
                    tempnum1 = [oldnum1;tempnum1];
                    tempnum2 = [oldnum2;tempnum2];
                    tempnum3 = [oldnum3;tempnum3];
                    tempresult = cell(size(tempnum1,1),1);
                    for ii = 1 : size(tempnum1,1)
                        tempresult1 = num2str(tempnum1(ii));
                        if length(tempresult1) == 1; tempresult1 = ['000', tempresult1];
                        elseif length(tempresult1) == 2; tempresult1 = ['00', tempresult1];
                        elseif length(tempresult1) == 3; tempresult1 = ['0', tempresult1];
                        end
                        tempresult2 = num2str(tempnum2(ii));
                        if length(tempresult2) == 1; tempresult2 = ['000', tempresult2];
                        elseif length(tempresult2) == 2; tempresult2 = ['00', tempresult2];
                        elseif length(tempresult2) == 3; tempresult2 = ['0', tempresult2];
                        end
                        tempresult3 = num2str(tempnum3(ii));
                        if length(tempresult3) == 1; tempresult3 = ['000', tempresult3];
                        elseif length(tempresult3) == 2; tempresult3 = ['00', tempresult3];
                        elseif length(tempresult3) == 3; tempresult3 = ['0', tempresult3];
                        end
                        tempresult{ii,1} = [tempresult1 ,'  ', tempresult2, '  ',tempresult3];
                    end
                    prompt_str = ['result_saved_history.',savedhistoryfields{idx_fields,1}, ' = tempresult;'];
                    eval(prompt_str);
                else
                    prompt_str = ['result_saved_history.',savedhistoryfields{idx_fields,1},' = [getfield(result_saved_history,savedfields{idx_fields,1}); getfield(saved_history, savedfields{idx_fields,1})];'];
                    eval(prompt_str);
                end
            end
        end
    end
end
saved = result_saved;
saved_history = result_saved_history;
save(savefn, 'saved', 'saved_history', 'saved_autoset');

end