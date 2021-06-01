function [partition, sumdegree] = greedy_partition(incidence_list, N, R, K)
    C = ceil(R/K);
    state_mat = zeros(C, N); 
    degree_max = zeros(1, N);
    partition = zeros(R, 1);
    partsize = zeros(C,1);
    picklist = randperm(R);
    for i = 1:R,
        templist = incidence_list{picklist(i)};
        mindelta = inf;
        pos = 0;
        for j = 1:C,
            if partsize(j) < K,
                delta = 0;
                for k = 1:length(templist),
                    if state_mat(j,templist(k)) == degree_max(templist(k)),
                        delta = delta + 1;
                    end
                end
                if delta < mindelta,
                    mindelta = delta;
                    pos = j;
                else if delta == mindelta && partsize(j) < partsize(pos),
                        pos = j;
                    end
                end
            end
        end
        partsize(pos) = partsize(pos) + 1;
        partition(picklist(i)) = pos;
        for k = 1:length(templist),
            state_mat(pos,templist(k)) = state_mat(pos,templist(k)) + 1;
            degree_max(templist(k)) = max(degree_max(templist(k)),...
                state_mat(pos,templist(k)));
        end
    end
    sumdegree = sum(degree_max);
end