function [input] = append_input_reversible_matched(input)

input.matched  = [];
for i =1:length(input.rxns)
        % find reverse rxn
          r = input.rxns{i};
          isBackward = ~isempty(regexp(r,'\_b','ONCE'));
          isForward = ~isempty(regexp(r,'\_f','ONCE'));
          if ~isBackward && ~isForward
        % r is irreversible and forward
            r_new = NaN;
           
        
          elseif isBackward
             x = strcat(cell2mat(regexp(r,'R\d*','match')),'_f');
             r_new = find(strcmp(input.rxns,x));
        
          elseif isForward
             x = strcat(cell2mat(regexp(r,'R\d*','match')),'_b');
             r_new = find(strcmp(input.rxns,x));

          else
            error('none_matched');
        
          end
        input.matched = [input.matched,r_new];
 end
end