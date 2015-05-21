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

% Evaluate 
Evaluating(population,map)
% Selection and Crossover
Selecting(population,gaConfig,0.5);
 % for i = 1:gaConfig.TournamentSize:gaConfig.PopulationSize
 %        %% TOURNAMENT SELECTION
 %        i1 = TournamentSelect(fitness,tournamentSelectionParameter,tournamentSize);
 %        i2 = TournamentSelect(fitness,tournamentSelectionParameter,tournamentSize);
 %        chromosome1 = population(i1,:);
 %        chromosome2 = population(i2,:);

 %        %% CROSS-OVER
 %        r = rand;
 %        if ( r < crossoverProbability)
 %            newChromosomePair = Cross(chromosome1, chromosome2);
 %            newPopulation(i,:) = newChromosomePair(1,:);
 %            newPopulation(i+1,:) = newChromosomePair(2,:);
 %        else
 %            newPopulation(i,:) = chromosome1;
 %            newPopulation(i+1,:) = chromosome2;
 %        end
 %    end


% Mutate
Mutating(population,map,gaConfig)
generation=200;
for i=1:generation
	Evaluating(population,map)
	Selecting(population,gaConfig,0.5)
	Mutating(population,map,gaConfig)
	Ploting(population,map)
 	% pause(0.5)
end


plot(map.mission_location(:,1),map.mission_location(:,2),'.')
