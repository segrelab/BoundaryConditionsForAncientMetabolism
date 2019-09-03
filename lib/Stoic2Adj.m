

function [A,k,R] = Stoic2Adj(S,rxn_iter,thresh)

% remove metabolites from adjacency
k = sum(~~S')./size(S,2) > thresh;
S(k,:) = [];

s = num2cell(~~S,2);
R = zeros(length(s));
A = zeros(length(s));

    for i =1:length(s)
        for j=i+1:length(s)
            t = s{i} & s{j};
            A(i,j) = any(t);
            A(j,i) = A(i,j);
            if A(i,j)
                
               R(i,j) = min(rxn_iter(t));
               R(j,i) = R(i,j);

        
            end
            
        end
    end
    
end