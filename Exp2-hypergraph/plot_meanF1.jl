using MAT
using StatsBase
std = StatsBase.std
trials = 10
include("../include/HyperLocal/HyperLocal.jl")
include("../include/HyperLocal/SparseCardHyperLocal.jl")
include("../src/hypergraph-clustering-utils.jl")
include("load_stackoverflow.jl")

## Load Delta-linear results
Druns = zeros(45,trials)
Df1s = zeros(45,trials)

for trial = 1:trials
    M = matread("Output/deltalinear_45_$trial.mat")
    stats = M["delta_stats"]
    drun = stats[:,6]
    df1 = stats[:,3]
    Druns[:,trial] = drun
    Df1s[:,trial] = df1
end

d_run_mean = mean(Druns,dims = 2)
d_run_std = std(Druns,dims = 2)
d_f1_mean = mean(Df1s,dims = 2)
d_f1_std = std(Df1s,dims = 2)

## Other methods
run_mean = zeros(45,9)
run_std = zeros(45,9)

f1_mean = zeros(45,9)
f1_std = zeros(45,9)

types = ["clique","x9","sqrt"]

epsis = 3
for tt = 1:3
    type = types[tt]

    for jj = 1:epsis
        Druns = zeros(45,trials)
        Df1s = zeros(45,trials)
        for trial = 1:trials
            M = matread("Output/SCHL_$(type)_45_$trial.mat")
            stats = M["sparse$(jj)_stats"]
            drun = stats[:,6]
            df1 = stats[:,3]
            Druns[:,trial] = drun
            Df1s[:,trial] = df1
        end
        # println("$type \t $jj")
        index = 3*(tt-1)+jj
        run_mean[:,index] = mean(Druns,dims = 2)
        run_std[:,index] = std(Druns,dims = 2)
        f1_mean[:,index] = mean(Df1s,dims = 2)
        f1_std[:,index] = std(Df1s,dims = 2)
    end
end



## Mean results

F1s = [d_f1_mean'; f1_mean[:,1]'; f1_mean[:,4]'; f1_mean[:,7]']
BestF1 = zeros(Int64,size(F1s,1))
for i = 1:45
    a1 = F1s[1,i]
    a2 = F1s[2,i]
    a3 = F1s[3,i]
    a4 = F1s[4,i]
    if a1 == a2 || a1 == a3 || a1 == a4 || a2 == a3 || a2 == a4 || a3 == a4
        println("Same results for two methods")
    end
    mx, mxi = findmax(F1s[:,i])
    BestF1[mxi] += 1
end

## Plot mean f1 scores in each case
s1 = 700
s2 = 400
ms = 6
p = sortperm(vec(d_f1_mean))
pl = scatter(d_f1_mean[p],markersize = ms, color = :blue,markerstrokecolor = :blue,legend = :topleft,label="δ-linear" )
scatter!(pl,f1_mean[p,1],markersize = ms, color = :green,markerstrokecolor = :green,label ="clique")
scatter!(pl,f1_mean[p,4],markersize = ms, color = :orange,markerstrokecolor = :orange, label="x^{0.9}")
scatter!(pl,f1_mean[p,7],markersize = ms, color = :purple, markerstrokecolor = :purple,label="sqrt",size = (s1,s2), grid = false)

ln = "  ".*names[p]
stepy = 1
plot!(pl,xticks = (1:stepy:length(ln), ln[1:stepy:length(ln)]),xrotation = 40,
xtickfont=font(8),
ytickfont=font(11),
guidefont=font(12),
titlefont=font(10),
legendfont=font(9))

cl = f1_mean[:,1]
sum(cl .> d_f1_mean)
mean(cl - d_f1_mean)
pl
# savefig("Figures/MeanPlots.pdf")
