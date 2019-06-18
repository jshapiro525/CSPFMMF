close all
clear all
clc

load('BPmfmodes.mat')
pdif1lam=abs((Z-Z_true)./Z_true)*100;

graph_element_distro(pdif1lam,'','Number of Entries in Z','Percent Error','Distribution of Percent Error')
print('betpicFMerror','-depsc')