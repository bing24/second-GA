clear 
close all


clear all;
close all;


% Load inputs
map = Map('map_image.bmp','resolution',200,'hieght',2000)

map.show('border')
represent(map)

% encoding(map)
% GA config
gaConfig = Configuration();
gaConfig.MaximumIterations = 50;
gaConfig.PopulationSize = 50;
gaConfig.PopulationType = 'random';
gaConfig.CrossoverRate = 0.8;
gaConfig.MutationRate = 0.03;
gaConfig.TournamentSize=10;
gaConfig.mutationProbability=0.01;
gaConfig.numberOfReplications = 2;
% Generatue initial populations

population = InitializePopulation(map, gaConfig) ;

Evaluating(population,map)

population.cover

obj=population;

temp=zeros(obj.gene_length,gaConfig.PopulationSize);
for i= 1: gaConfig.PopulationSize
count=1;
while temp(obj.gene_length,i)==0
temp(count,i)=map.mission_num(ceil(rand*size(map.mission_num,1)));
count=count+1;
end
end