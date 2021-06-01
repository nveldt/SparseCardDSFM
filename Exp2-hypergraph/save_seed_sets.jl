using MAT
using Random
using StatsBase
using SparseArrays
include("../include/HyperLocal/HyperLocal.jl")
include("../include/HyperLocal/SparseCardHyperLocal.jl")
include("../src/hypergraph-cond-utils.jl")
## Read in the matrix, if it isn't read in already

s = time()
@time M = matread(homedir()*"/Google Drive/Data Temp/stackoverflow_answer_H.mat")
LabelMatrix = M["LabelMatrix"]
LabelNames = M["LabelNames"]
MainLabels = M["MainLabels"]
H = M["H"]
order = round.(Int64,vec(sum(H,dims=2)))
d = vec(sum(H,dims=1))
volA = sum(d)
m,n = size(H)
Ht = sparse(H')
toc = time()-s
Edges = incidence2elist(H)
println("Done loading things into memory in $toc seconds.")


## Clusters: follows approach of Veldt et al. [KDD 2020]
Tconds = Vector{Float64}()
Tsizes = Vector{Int64}()
labels = Vector{Int64}()
for lab = 1:length(MainLabels)
    label = MainLabels[lab]
    T = findnz(LabelMatrix[:,label])[1]
    nT = length(T)
    condT, volT, cutT = tl_cond(H,T,d,1.0,volA,order)
    if condT < .2
        println("$nT \t $condT \t"*LabelNames[label])
        push!(Tconds,condT)
        push!(labels,label)
        push!(Tsizes,nT)
    end
end

## Save some seed sets
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
matwrite("Seed_sets_stack45.mat",Dict("R_stats"=>R_stats,"seed_sets"=>seed_sets))
