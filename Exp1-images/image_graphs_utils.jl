
function check_cutval(A,svec,tvec,S)
    """
    Check the cut value associated with the set of nodes in graph A.
    Source and sink nodes are not explicitly included, but terminal
    edges weights are given in svec and tvec.
    """
    cut = 0
    n = size(A,1)
    Sbar = setdiff(1:n,S)   # complement set

    # terminal edge weight penalties
    cut += sum(tvec[S])
    cut += sum(svec[Sbar])
    AS = A[S,Sbar]
    cut += sum(AS)

    return cut
end


function clique_expansion_penalties(k)
    """
    Symmetric clique splitting penalty for hyperedge of size k.
    Only gives first r = floor(k/2) penalties
    """
    @assert(k > 1)
    r = floor(Int64,(k)/2)
    I = collect(0:k)
    w = I .* (k .- I)
    w = w[1:r+1]
    return w
end


function load_superpixel_cliques(dataset,lambda,numsuperpix=500)
    """
    Load superpixel regions for a certain image dataset from experiements of
    Jegalka et al.

    Save the superpixel regions in hyperedge list format, with clique
    expansion splitting penalties, weighted by lambda
    (use lambda values as given previously by Jegelka et al.)
    """

    # Load superpixel region information
    M = matread("../data/data_reflections/$(dataset)_p$numsuperpix.mat")
    sup_img = M["sup_img"]

    col = size(sup_img,2)                   # number of columns
    row = size(sup_img,1)                   # number of rows
    node = (i,j)->i + (j-1)*row             # go from (i,j) node ID to linear node ID
    region = length(unique(sup_img))        # number of superpixels

    # Save superpixels in alternate format
    Edges = Vector{Vector{Int64}}()
    for e = 1:region
        # each superpixel defines a hyperedge
        push!(Edges,Vector{Int64}())
    end

    EdgeLengths = zeros(Int64,region)
    for i = 1:row
        for j = 1:col
            v = node(i,j)    # linear index of node at gridpoint (i,j)

            # the superpixel, i.e., hyperedge, it's in
            e = round(Int64,sup_img[i,j])
            push!(Edges[e],v)
            EdgeLengths[e] += 1
        end
    end

    # Discard cliques with less than 2 nodes, they do nothing
    keep = findall(x->x>1,EdgeLengths)
    Edges = Edges[keep]
    EdgesW = Vector{Vector{Float64}}()
    for ei = 1:length(Edges)
        # We'll use clique expansion splitting penalties,
        # scaled with previously determined weights of Jegelka et al.
        e = Edges[ei]
        push!(EdgesW,lambda*clique_expansion_penalties(length(e)))
    end

    return Edges, EdgesW
end


function eval_fS_image(Edges,EdgesW,svec,tvec,A,S)
    """
    For a set S, evaluate the objective from a decomposable
    submodular function involving terminal edges as defined
    in svec and tvec, pairwise non-terminal edges
    as defined in A, and higher-order terms from Edges
    and EdgesW.
    """
    n = size(A,1)
    eS = zeros(n)
    eS[S] = ones(length(S))
    Sbar = setdiff(1:n,S)   # complement set

    # terminal edge weight penalties
    cutval = sum(tvec[S])
    cutval += sum(svec[Sbar])

    # interior edge penalties
    AS = A[S,Sbar]
    cutval += sum(AS)

    # higher-order cut penalties
    for i = 1:length(Edges)
        e = Edges[i]
        w = EdgesW[i]
        t = round(Int64,sum(eS[e]))

        # do this because it's a symmetric splitting function,
        # so w is only half the size
        t = min(t,length(e)-t)
        cutval += w[t+1]
    end

    return cutval

end

function sweep_cut_f(Edges,EdgesW,svec,tvec,A,y)
    if length(y) > 10
        println("This code is inefficient. Not recommended for use for vectors this big.")
        return
    end
     Order = sortperm(y,rev = true)
     n = length(y)
     @assert(n == size(A,1))
     Sbest = [Order[1]]
     fbest = eval_fS_image(Edges,EdgesW,svec,tvec,A,Sbest)
     for i = 2:n
         S = Order[1:i]
         f = eval_fS_image(Edges,EdgesW,svec,tvec,A,S)
         if f < fbest
             fbest = f
             Sbest = S
         end
     end

     return Sbest, fbest
 end



function Objective_from_Primals(Xmat,sweepi,Edges,EdgesW,svec,tvec,An)
    """
    Given a set of primal variables from a continuous DSFM method,
    stored as columns of Xmat, find the objective score resulting
    from the best level set S_lam = {i : x_i >= lam}.

    The value of sweepi tells us the best level set as
    determined by the Matlab/C++ code that Xmat comes from.
    """
    sourceside = sum(svec)
    t = size(Xmat,2)
    @assert(t == length(sweepi))
    Fvals = zeros(t)
    for i = 1:t
        x = Xmat[:,i]
        Order = sortperm(x,rev=true)
        wherei = findfirst(x->x==sweepi[i]+1,Order)
        S = Order[1:wherei]
        Fvals[i] = eval_fS_image(Edges,EdgesW,svec,tvec,An,S)
    end
    return Fvals
end
