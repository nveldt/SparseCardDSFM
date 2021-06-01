using MAT
include("../src/hypergraphs.jl")
include("../src/SparseCard.jl")

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

# number of columns
col = size(sup_img,2)

# number of rows
row = size(sup_img,1)

# go from (i,j) node ID to linear node ID
node = (i,j)->i + (j-1)*row

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
N = length(unary)
A = sparse(I,J,W,N,N)
A = A+sparse(A')
matwrite("Graph_noclique_$(dataset)_$numsuperpix.mat",Dict("A"=>A,"svec"=>svec, "tvec"=>tvec))

end
end
