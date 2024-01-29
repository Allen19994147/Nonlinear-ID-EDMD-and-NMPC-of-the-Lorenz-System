%%% A library of candidate functions for nonlinear system identification
%%% You may extend this library to fit your use.

%%% X should be a m by n matrix, with m time sample state data, n state
function Theta = Candidate_Library(X)
[m,n] = size(X);

% Theta = [X(:,1) X(:,1).*X(:,1) X(:,2) X(:,2).^3];

% Build library
Theta = X;

% 2nd order
for i=1:n
    for j=i:n
        Theta = [Theta X(:,i).*X(:,j)];
    end
end

% 3rd order
for i=1:n
    for j=i:n
        for k=j:n
            Theta = [Theta X(:,i).*X(:,j).*X(:,k)];
        end
    end
end
% 4rd order
for i=1:n
    for j=i:n
        for k=j:n
            for ll=k:n
                Theta = [Theta X(:,i).*X(:,j).*X(:,k).*X(:,ll)];
            end
        end
    end
end

% sin & cos, freq up to w
for w = 1:5
    for i=1:n
        Theta = [Theta sin(w.*X(:,i)) cos(w.*X(:,i))];
    end
end


