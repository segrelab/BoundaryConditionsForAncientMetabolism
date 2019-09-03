function [nam] = reducedNameSet(nam)

for i =1:length(nam)
j = splitString(nam{i},';');
nam{i} = j{1};
end
end