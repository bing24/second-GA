classdef InitializePopulation < handle

properties
cover
chromo
chromo_locationx
chromo_locationy
chromo_number
fitness
minimumFitness
bestIndividualIndex
newPopulation
gene_length=50;
start_point=14;
end_point=124;
end

methods

function population = InitializePopulation(map, gaConfig) 
population.chromo_number=population.gene_length;

population.chromo=zeros(population.gene_length,gaConfig.PopulationSize);
population.chromo(1,:)=population.start_point;
population.chromo(end,:)=population.end_point;
for i= 1: gaConfig.PopulationSize

count=2;

while population.chromo(population.gene_length-1,i)==0
population.chromo(count,i)=map.mission_num(ceil(rand*size(map.mission_num,1)));
count=count+1;
end
%Sort the chromosome
for j=1:size(population.chromo,1)

population.chromo_locationx(j,i)=map.location_matrix(population.chromo(j,i),1);
population.chromo_locationy(j,i)=map.location_matrix(population.chromo(j,i),2);
end
end
% Plot the current chromosome
% figure
% plot(obj.chromo_location(:,1),obj.chromo_location(:,2),'.')
% hold on 
% contour(obj.matrix)


end % function
function  Evaluating (obj,map)
% Calculate Euclidean between each nodes
for j= 1: size(obj.chromo,2)
for i= 1: obj.chromo_number-1
nodes_dis(i)=norm([obj.chromo_locationx(i,j),obj.chromo_locationy(i,j)]-[obj.chromo_locationx(i+1,j),obj.chromo_locationy(i+1,j)]);
% Insert charging period in each travel
% charging_number_node(i)=randi([floor(nodes_dis(i)/dis_per_battery_max),floor(nodes_dis(i)/dis_per_battery_min)]);
% if i==1
%     chromo_charging=obj.chromo(1);
% else
% chromo_charging=cat(1,chromo_charging,[obj.chromo(i);zeros(charging_number_node(i),1)]);
% end
end
% chromo_charging=[chromo_charging;obj.chromo(end)];
% Cost function (distance)

cost_dis(1,j)=sum(nodes_dis);
% Calculate path coverage
count=0;
for ii=1:size(map.mission_num,1)
	cc=map.mission_num(ii)~=obj.chromo(:,j);
	
	if sum(cc)~=size(obj.chromo,1)
		count=count+1;
	end
end
obj.cover(1,j)=count/size(map.mission_num,1);



obj.fitness(1,j)=cost_dis(1,j)/obj.cover(1,j);

% obj.fitness(1,j)=cost_dis(1,j);
[obj.minimumFitness, obj.bestIndividualIndex] = min(obj.fitness); 

end
fprintf('Minimum Fitness: %d\n',obj.minimumFitness);
fprintf('Minimum Fitness index: %d\n',obj.bestIndividualIndex);
end
function Selecting (obj,gaConfig,tournamentSelectionParameter)
    obj.newPopulation=obj.chromo;
	for i = 1:gaConfig.TournamentSize:gaConfig.PopulationSize

% iSelected = TournamentSelect( fitnessValues, tournamentSelectionParameter, tournamentSize)
%select 'tournamentSize' candidates for tournament
candidates = 1 + fix(rand(1,gaConfig.TournamentSize)*gaConfig.PopulationSize);
candidateFitnesses = obj.fitness(candidates);
[~, sortedIndexes] = sort(candidateFitnesses,1,'descend');
selectionProbabilityMatrix = tournamentSelectionParameter*((1-tournamentSelectionParameter).^(0:gaConfig.TournamentSize-2)');
r = rand;
iSelected = candidates(sortedIndexes(r>selectionProbabilityMatrix));
if isempty(iSelected)
    iSelected = candidates(sortedIndexes(end));
else
    iSelected = iSelected(1);
end

%% TOURNAMENT SELECTION
        
        chromosome1 = obj.chromo(:,iSelected);
        candidates = 1 + fix(rand(1,gaConfig.TournamentSize)*gaConfig.PopulationSize);
candidateFitnesses = obj.fitness(candidates);
[~, sortedIndexes] = sort(candidateFitnesses,1,'descend');
selectionProbabilityMatrix = tournamentSelectionParameter*((1-tournamentSelectionParameter).^(0:gaConfig.TournamentSize-2)');
r = rand;
iSelected = candidates(sortedIndexes(r>selectionProbabilityMatrix));
if isempty(iSelected)
    iSelected = candidates(sortedIndexes(end));
else
    iSelected = iSelected(1);
end
        chromosome2 = obj.chromo(:,iSelected);

        %% CROSS-OVER
        r = rand;
        if ( r < gaConfig.CrossoverRate)

            % newChromosomePair = Cross(chromosome1, chromosome2);
            n=round(rand(size(obj.chromo,1),1));
            for ii=1:length(n)
if n(ii)==1
    newChromosomePair(ii,1)=chromosome1(ii);
    newChromosomePair(ii,2)=chromosome2(ii);
else
    newChromosomePair(ii,1)=chromosome2(ii);
    newChromosomePair(ii,2)=chromosome1(ii);
end
            end

            obj.newPopulation(:,i) = newChromosomePair(:,1);
            obj.newPopulation(:,i+1) = newChromosomePair(:,2);
        else
            obj.newPopulation(:,i) = chromosome1;
            obj.newPopulation(:,i+1) = chromosome2;
        end

end
end

function Ploting(obj,map)
    hold off
     plot(obj.chromo_locationx(:,obj.bestIndividualIndex),obj.chromo_locationy(:,obj.bestIndividualIndex),'r')
     hold on
 contour(map.matrix)
 title(obj.minimumFitness)

	end

	function Mutating(obj,map,gaConfig)
indexes = rand(size(obj.chromo))<gaConfig.mutationProbability  ;               % Index for Mutations

temp=zeros(obj.gene_length,gaConfig.PopulationSize);
for i= 1: gaConfig.PopulationSize
count=1;
while temp(obj.gene_length,i)==0
temp(count,i)=map.mission_num(ceil(rand*size(map.mission_num,1)));
count=count+1;
end
end
%fix the starting point finishing point and best chromosome
indexes(1,:)=0;
indexes(end,:)=0;
indexes(:,obj.bestIndividualIndex)=0;
obj.newPopulation(indexes) =temp(indexes);

% = tempPopulation(indexes)*-1+1;                     % Bit Flip Occurs
%% PRESERVATION OF PREVIOUS BEST SOLUTION
    bestChromosome =obj.chromo(:,obj.bestIndividualIndex);

randIndexes = ceil(rand(1,gaConfig.numberOfReplications).*size(obj.chromo,2));
obj.newPopulation(:,randIndexes) = repmat(bestChromosome,1,gaConfig.numberOfReplications);
obj.chromo=obj.newPopulation;

for i= 1: gaConfig.PopulationSize
for j=1:size(obj.chromo,1)

obj.chromo_locationx(j,i)=map.location_matrix(obj.chromo(j,i),1);
obj.chromo_locationy(j,i)=map.location_matrix(obj.chromo(j,i),2);
end
end
	end

end % method

end % classdef