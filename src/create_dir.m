function [jobpath,figpath,datpath]=create_dir(filename,time)
% Create directories based on time and filename of .inp.
jobpath = strcat('output/',filename,'-',time);
figpath = strcat(jobpath,'/','figs');
datpath = strcat(jobpath,'/','dats');
mkdir(jobpath);
mkdir(figpath);
mkdir(datpath);
%Redirect command line to .txt file
diary(strcat(jobpath,'/','log.txt'));
diary on; 