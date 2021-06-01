using MAT
using Random
using StatsBase
using SparseArrays
include("../include/HyperLocal/HyperLocal.jl")
include("../include/HyperLocal/SparseCardHyperLocal.jl")
include("../src/hypergraph-clustering-utils.jl")

## Load in the data
include("load_stackoverflow.jl")

## Generate seed sets
# NOTE: original seed sets are saved in Seeds folder

# include("generate_stackoverflow_seeds.jl")

## Run Hyperlocal+SparseCard with alternative hyperedge cut penalties
types = ["clique", "x9", "sqrt"]
lnum = length(labels)

for trial = 1:10
    SS = matread("Seeds/Seed_sets_stack45_$trial.mat")
    seed_sets = SS["seed_sets"]

    # Run for the three different alternate penalties
    for tt = 1:3
        epsilon = 1.0
        sparsityeps1 = 1.0
        sparsityeps2 = 0.1
        sparsityeps3 = 0.01
        type = types[tt]
        output = "Output/SCHL_$(type)_45_$trial.mat"

        if tt == 1
            EdgesW = clique_expansion_EdgesW(order)
        elseif tt == 2
            fun = x->x.^(.9)
            EdgesW = generic_EdgesW(order,fun)
        else
            funsqrt = x->sqrt.(x)
            EdgesW = generic_EdgesW(order,funsqrt)
        end

        sparse1_stats = zeros(lnum,6)
        sparse2_stats = zeros(lnum,6)
        sparse3_stats = zeros(lnum,6)

        for lab = 1:length(LabelNames)
            T = findnz(LabelMatrix[:,lab])[1]
            nT = length(T)
            println("$lab \t $nT \t $(LabelNames[label])")

            # Generate a new seed set
            seednum = round(Int64,nT/20)
            grownum = round(Int64,min(nT*2,n/2-seednum))
            Rstart = findall(x->x>0,seed_sets[:,lab])
            OneHop = get_immediate_neighbors(H,Ht,Rstart)
            Rmore = BestNeighbors(H,d,Rstart,OneHop,grownum)
            R = union(Rmore,Rstart)
            Rs = findall(x->in(x,Rstart),R)     # Force seed nodes to be in output set

            # Run HyperLocal + SparseCard with 3 different sparsity epsilons
            sparsities = [1.0, 0.1, 0.01]
            for sss = 1:3
                sparsityeps = sparsities[sss]
                s = time()
                S, lcond = SparseCardHyperLocal(H,Ht,EdgesW,Edges,sparsityeps,order,d,R,epsilon,Rs)
                timer = time()-s
                pr, re, f1 = PRF(T,S)
                nS =  length(setdiff(S,R))
                if sss == 1
                    sparse1_stats[lab,:] = [pr, re, f1, length(S), nS, timer]
                elseif sss == 2
                    sparse2_stats[lab,:] = [pr, re, f1, length(S), nS, timer]
                else
                    sparse3_stats[lab,:] = [pr, re, f1, length(S), nS, timer]
                end
            end
            # matwrite(output,Dict("sparse1_stats"=>sparse1_stats,"sparse2_stats"=>sparse2_stats,"sparse3_stats"=>sparse3_stats))
        end
    end

    # Run Hyperlocal with standard delta linear threshold penalty
    delta = 5000.0
    epsilon = 1.0
    delta_stats = zeros(lnum,6)
    output = "Output/deltalinear_45_$trial.mat"
    for lab = 1:length(labels)
        label = labels[lab]
        T = findnz(LabelMatrix[:,label])[1]
        nT = length(T)
        println("$lab \t $nT \t $label")

        # Generate a new seed set
        seednum = round(Int64,nT/20)
        grownum = round(Int64,min(nT*2,n/2-seednum))
        Rstart = findall(x->x>0,seed_sets[:,lab])
        OneHop = get_immediate_neighbors(H,Ht,Rstart)
        Rmore = BestNeighbors(H,d,Rstart,OneHop,grownum)
        R = union(Rmore,Rstart)
        Rs = findall(x->in(x,Rstart),R)     # Force seed nodes to be in output set

        # Run HyperLocal with delta = 5000.0
        s = time()
        S, lcond = HyperLocal(H,Ht,order,d,R,epsilon,delta,Rs,true)
        timer = time()-s
        pr, re, f1 = PRF(T,S)
        nS =  length(setdiff(S,R))
        delta_stats[lab,:] = [pr, re, f1, length(S), nS, timer]
        # matwrite(output,Dict("delta_stats"=>delta_stats))
    end

end
