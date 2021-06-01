for trial = 1:10
    lnum = length(labels)
    R_stats = zeros(lnum,3)
    seed_sets = spzeros(n,lnum)
    for lab = 1:length(labels)
        label = labels[lab]
        T = findnz(LabelMatrix[:,label])[1]
        nT = length(T)
        println("$lab \t $nT \t $label")

        # Generate a new seed set
        seednum = round(Int64,nT/20)
        grownum = round(Int64,min(nT*2,n/2-seednum))
        p = randperm(nT)
        Rstart = T[p[1:seednum]]
        OneHop = get_immediate_neighbors(H,Ht,Rstart)
        Rmore = BestNeighbors(H,d,Rstart,OneHop,grownum)
        R = union(Rmore,Rstart)
        Rs = findall(x->in(x,Rstart),R)     # Force seed nodes to be in output set
        prr, rer, f1r = PRF(T,R)
        R_stats[lab,:] = [prr, rer, f1r]
        seed_sets[Rstart,lab] .= 1
    end
# matwrite("Seed_sets_stack45_$trial.mat",Dict("R_stats"=>R_stats,"seed_sets"=>seed_sets))
end
