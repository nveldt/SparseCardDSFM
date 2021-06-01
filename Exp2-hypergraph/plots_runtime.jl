using MAT
using StatsBase
std = StatsBase.std
trials = 10
include("../include/HyperLocal/HyperLocal.jl")
include("../include/HyperLocal/SparseCardHyperLocal.jl")
include("../src/hypergraph-cond-utils.jl")
## Read in the matrix, if it isn't read in already
include("load_stackoverflow.jl")

## Mean runtime plots
for type = ["sqrt", "x9", "clique"]
typename = type
if type == "x9"
    typename = "x^{0.9}"
end
run_mean = zeros(45,3)
run_std = zeros(45,3)
f1_mean = zeros(45,3)
f1_std = zeros(45,3)
for jj = 1:3
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
    index = jj
    run_mean[:,index] = mean(Druns,dims = 2)
    run_std[:,index] = std(Druns,dims = 2)
    f1_mean[:,index] = mean(Df1s,dims = 2)
    f1_std[:,index] = std(Df1s,dims = 2)
end

function plotunif(pl, y, label, color,ms)
    scatter!(pl,y,markersize = ms, markerstrokecolor = color, color = color, markershape = :circle, label = label)
end
stepy = 2
s1 = 500
s2 = 350
p = sortperm(run_mean[:,1])
ln = names[p]
pl = plot(grid = false)
plotunif(pl,run_mean[p,3],"ε = 0.01",:black,5)
plotunif(pl,run_mean[p,2],"ε = 0.1",:grey,5)
plotunif(pl,run_mean[p,1],"ε = 1.0",:green,5)
plot!(pl,yscale = :log10,legend = :topleft,title = type*" function, mean runtime",
    ylabel = "Runtime (seconds)",size = (s1,s2))
mean(abs.(f1_mean[:,1]-f1_mean[:,2]))
plot!(pl,xticks = (1:stepy:length(ln), ln[1:stepy:length(ln)]),xrotation = 40,
xtickfont=font(8),
ytickfont=font(11),
guidefont=font(12),
titlefont=font(12),
legendfont=font(9))

savefig("Figures/$(type)_runtime_plot.pdf")
jump = mean(run_mean[:,2]./run_mean[:,1])
jump2 = mean(run_mean[:,3]./run_mean[:,2])
println("")
println("Runtime jumps = $jump, $jump2")
end


## Meanwhile, F1 doesn't change
for type = ["sqrt", "x9", "clique"]
run_mean = zeros(45,3)
run_std = zeros(45,3)
f1_mean = zeros(45,3)
f1_std = zeros(45,3)
for jj = 1:3
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
    index = jj
    run_mean[:,index] = mean(Druns,dims = 2)
    run_std[:,index] = std(Druns,dims = 2)
    f1_mean[:,index] = mean(Df1s,dims = 2)
    f1_std[:,index] = std(Df1s,dims = 2)
end

function plotunif(pl, y, label, color,ms)
    scatter!(pl,y,markersize = ms, markerstrokecolor = color, color = color, markershape = :circle, label = label)
end
stepy = 2
s1 = 500
s2 = 350
p = sortperm(run_mean[:,1])
ln = names[p]
pl = plot(grid = false)
plotunif(pl,f1_mean[p,3],"ε = 0.01",:black,5)
plotunif(pl,f1_mean[p,2],"ε = 0.1",:grey,5)
plotunif(pl,f1_mean[p,1],"ε = 1.0",:green,5)
typename = type
if type == "x9"
    typename = "x^{0.9}"
end
plot!(pl,yscale = :identity,legend = :topleft,title = typename*" function, mean F1",
    ylabel = "F1 Score",size = (s1,s2))
mean(abs.(f1_mean[:,1]-f1_mean[:,2]))
plot!(pl,xticks = (1:stepy:length(ln), ln[1:stepy:length(ln)]),xrotation = 40,
xtickfont=font(8),
ytickfont=font(11),
guidefont=font(12),
titlefont=font(12),
legendfont=font(9))

savefig("Figures/$(type)_f1_plot.pdf")
end
