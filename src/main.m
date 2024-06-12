% Phase field fracture
%% Start up
clc
clear all
format long
sympref('HeavisideAtOrigin','default');
warning('off','MATLAB:nearlySingularMatrix');
%% Pipeline
pipeline={};
if isempty(pipeline)
    solve('script');
else
    for i=1:length(pipeline)
%         try
            solve(cell2mat(pipeline(i)));
%         catch ME
%             warning(cell2mat(strcat('An error occurs in project ',{32},cell2mat(pipeline(i)))));
%             disp(ME.message);
%         end
    end
end