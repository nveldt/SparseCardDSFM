## Reads in data for stackoverflow hypergraph from
# https://www.cs.cornell.edu/~arb/data/stackoverflow-answers/
# and stores it in .mat format with main clusters

include("../src/hypergraph-clustering-utils.jl")
## Read label names
pathname = joinpath(dirname(dirname(@__FILE__)),"data")
dataname = "stackoverflow-answers"
open("$pathname/$dataname/label-names-$dataname.txt") do f
    for line in eachline(f)
        push!(names, line)
    end
end
## Read labels, save in matrix format
using SparseArrays

Labels = Vector{Vector{Int64}}()
pathname = joinpath(dirname(dirname(@__FILE__)),"data")
open("$pathname/$dataname/node-labels-$dataname.txt") do f
    for line in eachline(f)
        labels = parse.(Int64,split(line,','))
        push!(Labels, labels)
    end
end
N = length(names)
LM = elist2incidence(Labels,N)


## Read in the hypergraph
Edges = Vector{Vector{Int64}}()
pathname = joinpath(dirname(dirname(@__FILE__)),"data")
maxN = 0
open("$pathname/$dataname/hyperedges-$dataname.txt") do f
    global maxN
    for line in eachline(f)
        edge = [parse(Int64, v) for v in split(line, ',')]
        sort!(edge)
        if 2 <= length(edge)
            push!(Edges,edge)
            n = maximum(edge)
            if n > maxN
                maxN = n
            end
        end
    end
end

Hmat = elist2incidence(Edges,maxN)
order = round.(Int64,vec(sum(Hmat,dims=2)))
d = vec(sum(Hmat,dims=1))
volA = sum(d)

## Get main clusters, following conditions outlined by Veldt et al. [KDD 2020]
Tconds = Vector{Float64}()
Tsizes = Vector{Int64}()
labels = Vector{Int64}()
num = 1
for lab = 1:length(names)
    global num
    T = findnz(LM[:,lab])[1]
    nT = length(T)
    if nT < 10000 && nT > 2000
        condT, volT, cutT = tl_cond(Hmat,T,d,1.0,volA,order)
        if condT < .2
            println("$num\t  $nT \t $condT \t"*names[lab])
            num+= 1
            push!(Tconds,condT)
            push!(labels,lab)
            push!(Tsizes,nT)
        end
    end
end

LabelMat = LM[:,labels]
LabelNames = names[labels]

## Save in .mat format
matwrite("Stackoverflow_data.mat", Dict("H"=>Hmat,"LabelMat"=>LabelMat,
                                "LabelNames"=>LabelNames,"ClusSizes"=>Tsizes,
                                "ClusConds"=>Tconds))
