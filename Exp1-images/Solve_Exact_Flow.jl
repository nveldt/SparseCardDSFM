using MAT
include("../src/hypergraphs.jl")
include("../src/SparseCard.jl")
include("../src/pwl_approx.jl")
include("image_graphs_utils.jl")
## Load

datasets = ["020_smallplant","001_oct"]


for numsuperpix = [500; 200]

    if numsuperpix == 200
        lambdas = [0.3; 0.1]
        lambda2 = 0.002

    elseif numsuperpix == 500
        lambdas = [0.08; 0.01]

        lambda2 = 0.008
    end

for jj = 1:2

dataset = datasets[jj]
lambda1 = lambdas[jj]

M = matread("../data/data_reflections/$(dataset)_p$numsuperpix.mat")

# The graph is a grid: each node is a pixel in a picture
# This gives superpixel regions
sup_img = M["sup_img"]

# Load
M2 = matread("../data/data_reflections/$(dataset)_wts.mat")
W1 = M2["W1"] # edge weights for vertical edges
W2 = M2["W2"] # edge weights for horizontal edges
unary = M2["unary"] # this is the s-t edge weights, already linearized

tic = time()
# number of columns
col = size(sup_img,2)

# number of rows
row = size(sup_img,1)

# go from (i,j) node ID to linear node ID
node = (i,j)->i + (j-1)*row

# number of superpixels
region = length(unique(sup_img))

## Save superpixels in alternate format
Edges = Vector{Vector{Int64}}()
for e = 1:region
    # each superpixel defines a hyperedge
    push!(Edges,Vector{Int64}())
end

for i = 1:row
    for j = 1:col
        v = node(i,j)    # linear index of node at gridpoint (i,j)

        # the superpixel, i.e., hyperedge, it's in
        e = round(Int64,sup_img[i,j])

        push!(Edges[e],v)
    end
end

## Build edge information
I = Vector{Int64}()
J = Vector{Int64}()
W = Vector{Float64}()

# Type 1 edges: terminal edges
tvec = abs.((unary .> 0).* unary)
svec = abs.((unary .< 0).* unary)

# Type 2 edges: vertical edges, organized into rows
for i = 1:row-1
    for j = 1:col
        v1 = node(i,j)      # first node in edge
        v2 = node(i+1,j)    # second node in edge
        push!(I,v1)
        push!(J,v2)

        # Use edge weights as defined by Jegelka et al.
        push!(W,lambda1*W1[i,j])
    end
end

# Type 3 edges: horizontal edges, organized into columns
for i = 1:col-1
    for j = 1:row
        v1 = node(j,i)      # first node in edge
        v2 = node(j,i+1)    # second node in edge
        push!(I,v1)
        push!(J,v2)

        # Use edge weights as defined by Jegelka et al.
        push!(W,lambda1*W2[j,i])
    end
end

## Type 4 edges: cliques from clique potentials
for e in Edges
    edge = sort(e)
    for i = 1:length(e)-1
        v1 = e[i]
        for j = i+1:length(e)
            v2 = e[j]
            push!(I,v1)
            push!(J,v2)
            push!(W,lambda2)
        end
    end
end

# turn to sparse mat
N = length(unary)
A = sparse(I,J,W,N,N)
A = A+sparse(A')
formgraph = time()-tic

# Solve the optimal flow problem on it
total_edges = length(A.nzval)
tic = time()
F = maxflow(A,vec(svec),vec(tvec),0.0)
opttime = time()-tic
opt = F.cutvalue
S = source_nodes_min(F,1e-10)[2:end].-1
cval = check_cutval(A,svec,tvec,S)
@show abs(opt-cval), opt, cval

matwrite("OptimalSolutions/Optimal_$(dataset)_$numsuperpix.mat",
    Dict("S"=>S,"opttime"=>opttime,"opt"=>opt,"total_edges"=>total_edges,"formgraph"=>formgraph))

end
end
