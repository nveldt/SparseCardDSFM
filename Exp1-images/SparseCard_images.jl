# Image segmentation experiments, but using sparse reduction + maxflow
using MAT
include("../src/hypergraphs.jl")
include("../src/SparseCard.jl")
include("image_graphs_utils.jl")

## Load optimal solution to the problem, for comparison

for dataset = ["020_smallplant","001_oct"]

for numsuperpix = [200 500]

if numsuperpix == 200
    lambda = 0.002
else
    lambda = 0.008
end

## Load edges from non-superpixel (clique) regions
G = matread("GraphPotentials/Graph_noclique_$(dataset)_$(numsuperpix).mat")
An = G["A"]
tvec = G["tvec"]
svec = G["svec"]
n = size(An,1)

In,Jn,Wn = findnz(An)

## Load hyperedges from superpixel (clique) regions

Edges, EdgesW = load_superpixel_cliques(dataset,lambda,numsuperpix)
K = maximum(length(e) for e in Edges)

## Sparse Reduction Runs

#epsilons = [1.0;0.5;0.1;0.05;0.01;0.005;0.001];0.0005;0.0001]
small = 4
epsilons = 10 .^ range( 0, -small, length = 20)
runs = length(epsilons)
objvals = zeros(runs)
approxcuts = zeros(runs)
runtimes = zeros(runs)
numedges = zeros(runs)
preptime = zeros(runs)
println("")
for i = 1:runs
    epsi = epsilons[i]
    epsilon = epsi*ones(K)
    tic = time()
    Ac,Ic,Jc,Wc = SymmetricCard_reduction(Edges,EdgesW,n,epsilon)
    preptime[i] = time()-tic
    N = size(Ac,1)
    I = [In;Ic]
    J = [Jn;Jc]
    W = [Wn;Wc]
    A = sparse(I,J,W,N,N)

    ## Run max flow
    sparse_edges = n + length(I)

    Svec = vec([svec;zeros(N-n)])
    Tvec = vec([tvec;zeros(N-n)])
    tic = time()
    F = maxflow(A,Svec,Tvec,0.0)
    flowtime = time()-tic
    S = source_nodes_min(F,1e-10)[2:end].-1
    ct = F.cutvalue
    cval = check_cutval(A,Svec,Tvec,S)
    @assert(abs(ct-cval)<1e-8)
    Strue = intersect(S,1:n)
    original_obj = eval_fS_image(Edges,EdgesW,svec,tvec,An,Strue)

    # timeratio = flowtime/opttime
    # approx = original_obj/opt
    # sparsity = sparse_edges/total_edges
    objvals[i] = original_obj
    approxcuts[i] = ct     # cut value in the approximate reduced graph
    numedges[i] = sparse_edges
    runtimes[i] = flowtime

    println("$epsi \t $original_obj \t $sparse_edges \t $flowtime")

end

## Save the output
matwrite("Output/SparseCard_$(dataset)_$(numsuperpix)_$small.mat",Dict("objvals"=>objvals,
"numedges" => numedges, "runtimes"=>runtimes,"preptime"=>preptime,
"epsilons"=>epsilons,"approxcuts"=>approxcuts))

end

end
