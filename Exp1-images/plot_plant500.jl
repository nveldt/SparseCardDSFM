# Image segmentation experiments, but using sparse reduction + maxflow
using MAT
using LinearAlgebra
using StatsBase
include("../src/SparseCard.jl")
include("image_graphs_utils.jl")

## Load optimal solution to the problem, for comparison
numsuperpix = 500

if numsuperpix == 200
    lambda = 0.002
    K = 13
else
    lambda = 0.008
    K = 16
end

dataset = "020_smallplant"
short = "plant"

G = matread("GraphPotentials/Graph_noclique_$(dataset)_$numsuperpix.mat")
An = G["A"]
tvec = G["tvec"]
svec = G["svec"]
n = size(An,1)
In,Jn,Wn = findnz(An)
sourceside = sum(svec)

Edges, EdgesW = load_superpixel_cliques(dataset,lambda,numsuperpix)

Opt = matread("OptimalSolutions/Optimal_$(dataset)_$numsuperpix.mat")
S_opt = Opt["S"]
opttime = Opt["opttime"]
total_edges = Opt["total_edges"]
opt = eval_fS_image(Edges,EdgesW,svec,tvec,An,S_opt)


## Load

pl = plot()

c1 = 10
c2 = 25
T = 600

Col = [1 0 0; 0 .75 0; 0 0 1; 1 .5 0]
timelimit = 62
methods = ["APcompact_T_$T", "RCDM_$K", "ACDM_$(K)_$c1","ACDM_$(K)_$c2","AP_T_$T", "RCDM_greedy_$K"]
methodname = ["IAP","RCDM", "ACDM","ACDM-tuned", "AP","RCDM-G"]

lastone = 4
for j = 1:lastone
# j = 2
col = RGBA(Col[j,1],Col[j,2],Col[j,3],1)
ms = 3.5
method = methods[j]
if j == 1
    jump = 10
else
    jump = 20
end

method = methods[j]
C = matread("ContOutput/ContOutput_$(short)$(numsuperpix)/$(method)_1.mat")
levelcost = vec(C["levelcost"])
L = length(levelcost)
Runs = zeros(5,L)
Objs = zeros(5,L)

for jj = 1:5
    C = matread("ContOutput/ContOutput_$(short)$(numsuperpix)/$(method)_$jj.mat")
    levelcost = vec(C["levelcost"])
    timevec = vec(C["timevec"])/1000
    Runs[jj,:] = timevec
    Objs[jj,:] = levelcost .+ sourceside
end

timevec = vec(mean(Runs,dims = 1))
Fvals = vec(mean(Objs,dims = 1))
Upper = vec(maximum(Objs,dims = 1))./opt .-1
Lower = vec(minimum(Objs,dims = 1))./opt .-1
Approx = Fvals./opt .- 1

xrange = 1:jump:length(Approx)

# Truncate if over time limit
# if maximum(timevec)>timelimit
#     lst = findfirst(x->x>timelimit,timevec)
#     st = rem(lst,jump)
#     timevec = timevec[1:lst]
#     Approx = Approx[1:lst]
#     Lower = Lower[1:lst]
#     Upper = Upper[1:lst]
#     xrange = 1:jump:lst
# end
#
# # Truncate and adjust if exact convergence is reached
# if minimum(Approx) <= 0
#     converged = findfirst(x->x<= 0,Approx)
#     # avoids log(0) errors without changing plot
#     Approx[converged] = 1e-12
#     Lower[converged] = 1e-12
#     Upper[converged] = 1e-12
#     st = rem(converged,jump)
#     xrange = st:jump:converged
# end

# Truncate if over time limit
if maximum(timevec)>timelimit
    lst = findfirst(x->x>timelimit,timevec)
    st = rem(lst,jump)
    timevec = timevec[1:lst]
    Approx = Approx[1:lst]
    Lower = Lower[1:lst]
    Upper = Upper[1:lst]
    xrange = 1:jump:lst
end

# Truncate and adjust if exact convergence is reached
if minimum(Approx) <= 0
    # converged = findfirst(x->x<= 0,Approx)
    converged = findfirst(x->x<= 0,Lower)
    # avoids log(0) errors without changing plot
    Approx[converged] = 1e-12
    Lower[converged] = 1e-12
    Upper[converged] = 1e-12
    st = rem(converged,jump)
    xrange = st:jump:converged
end

label = methodname[j]
alp = 1.0
plot!(pl,timevec[xrange],Approx[xrange],
    label = label,alpha = alp,
    markershape = :circle,markersize = ms,markerstrokewidth = 0.0,
    markercolor = col, markerstrokecolor = col, color = col,yscale = :log10)

plot!(pl,timevec[xrange],Lower[xrange],alpha = 0,
    fillrange=Upper[xrange], fillalpha=0.3,
    markerstrokecolor = col, color = col,linewidth = 0.0,label = "")

end
pl

## Sparsecard
ms = 3.5
M = matread("Output/SparseCard_$(dataset)_$(numsuperpix)_4.mat")
runtimes = M["runtimes"]
objvals = M["objvals"]
epsilons = M["epsilons"]
flow_approx = (objvals./opt) .- 1

p = sortperm(runtimes)
plot!(pl, runtimes[p[1:end-2]], flow_approx[p[1:end-2]], label = "SparseCard",
    xlabel = "Runtime (seconds)",ylabel = "Approx Ratio - 1",yscale = :log10,
    markershape = :circle,markersize = ms,markerstrokewidth = 0.0,
    markercolor = :black, color = :black)
# @show epsilons[p[end-2]]

# annotate!(pl, runtimes[p[end-2]]-5,flow_approx[p[end-2]]+1e-5," ε = 0.0003",9 )

## Different size
s1 = 300
s2 = 250
gfs = 12
tfs = 12
plot!(pl,xlimit = [0,timelimit+1],ylimit = [1e-6,1],yscale = :log10,
        grid = false, size = (s1,s2),fontfamily = "helvetica",
        title = "\n plant, $numsuperpix superpixels",
        legend = false,tickfontsize=tfs,guidefontsize = gfs)

savefig("Figures/RunApprox_$(dataset)_$(numsuperpix).pdf")


## Table results

runtimes = round.(M["runtimes"],digits = 1)
approxes = round.(M["objvals"]./opt .- 1,digits = 9)
epsilons = round.(M["epsilons"],digits = 4)
numedges = M["numedges"]
sparsity = round.(numedges./total_edges,digits = 3)


## Table

for i = 1:4
    if i == 1
        v = epsilons
        print("\$\\varepsilon\$ &")
    elseif i == 2
        v = approxes
        print("Approx. &")
    elseif i == 3
        v = sparsity
        print("Sparsity &")
    else
        v = runtimes
        print("Runtime &")
    end
    for j = 1:3:20
        if i == 2
            print("$(v[j]) &")
        else
            print("$(v[j]) &")
        end
    end
    if i == 2
        println("$(v[end]) \\\\")
    else
        println("$(v[end]) \\\\")
    end
end
