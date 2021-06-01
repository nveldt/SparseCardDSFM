using MAT
using StatsBase
std = StatsBase.std
include("../include/HyperLocal/HyperLocal.jl")
include("../include/HyperLocal/SparseCardHyperLocal.jl")
include("../src/hypergraph-clustering-utils.jl")

## Read in the matrix, if it isn't read in already
include("load_stackoverflow.jl")

## Delta linear results
p = sortperm(Tsizes)
Rdel = zeros(45,trials)
Fdel = zeros(45,trials)

for trial = 1:trials
    M = matread("Output/deltalinear_45_$trial.mat")
    stats = M["delta_stats"]
    drun = stats[:,6]
    df1 = stats[:,3]
    Rdel[:,trial] = drun
    Fdel[:,trial] = df1
end
Rdel = Rdel[p,:]
Fdel = Fdel[p,:]
del_run_mean = mean(Rdel,dims = 2)
del_run_std = std(Rdel,dims = 2)
del_f1_mean = mean(Fdel,dims = 2)
del_f1_std = std(Fdel,dims = 2)

## Other results
Rx9 = zeros(45,trials)
Rsq = zeros(45,trials)
Rcl = zeros(45,trials)
Fx9 = zeros(45,trials)
Fsq = zeros(45,trials)
Fcl = zeros(45,trials)
types = ["sqrt", "x9", "clique"]
for tt = 1:3
    type = types[tt]
    for trial = 1:trials
        M = matread("Output/SCHL_$(type)_45_$trial.mat")
        stats = M["sparse1_stats"]
        drun = stats[:,6]
        if tt ==1
            Rsq[:,trial] = drun
        elseif tt == 2
            Rx9[:,trial] = drun
        else
            Rcl[:,trial] = drun
        end
        df = stats[:,3]
        if tt ==1
            Fsq[:,trial] = df
        elseif tt == 2
            Fx9[:,trial] = df
        else
            Fcl[:,trial] = df
        end
    end
end

Rx9 = Rx9[p,:]
Rsq = Rsq[p,:]
Rcl = Rcl[p,:]
Fx9 = Fx9[p,:]
Fsq = Fsq[p,:]
Fcl = Fcl[p,:]
## Write to file
open("tableoutput.txt", "w") do io
for ntab = 1:4
    if ntab == 1
        numcount = 9
        irange = 1:9
        size1 = Tsizes[p[1]]
        size2 = Tsizes[p[9]]
    else
        st = 10+(ntab-2)*12
        finish = st+11
        irange = st:finish
        # @show irange
        size1 = Tsizes[p[st]]
        size2 = Tsizes[p[finish]]
    end

        write(io,"\\begin{table}[t]
    	\\caption{Results for Stackoverflow clusters from size $size1 to size $size2.}
    	\\label{tab:stackmore$ntab}
    	\\centering
    	\\begin{tabular}{l l l l l l}
    		\\toprule
    		Cluster & Size & Penalty & F1 & Runtime  &\\# Best \\\\\n")


for i = irange

x9_mean_f1 = mean(Fx9[i,:])
cl_mean_f1 = mean(Fcl[i,:])
sq_mean_f1 = mean(Fsq[i,:])
del_mean_f1 = mean(Fdel[i,:])

x9_std_f1 = std(Fx9[i,:])
cl_std_f1 = std(Fcl[i,:])
sq_std_f1 = std(Fsq[i,:])
del_std_f1 = std(Fdel[i,:])

x9_mean_run = mean(Rx9[i,:])
cl_mean_run = mean(Rcl[i,:])
sq_mean_run = mean(Rsq[i,:])
del_mean_run = mean(Rdel[i,:])

x9_std_run = std(Rx9[i,:])
cl_std_run = std(Rcl[i,:])
sq_std_run = std(Rsq[i,:])
del_std_run = std(Rdel[i,:])

runs_m = [del_mean_run; cl_mean_run; sq_mean_run; x9_mean_run]
runs_std = [del_std_run; cl_std_run; sq_std_run; x9_std_run]

f1_m = [del_mean_f1; cl_mean_f1; sq_mean_f1; x9_mean_f1]
f1_std = [del_std_f1; cl_std_f1; sq_std_f1; x9_std_f1]

F1s = [Fdel[i,:]'; Fcl[i,:]'; Fsq[i,:]'; Fx9[i,:]']
BestF1 = zeros(Int64,size(F1s,1))
for i = 1:10
    mx, mxi = findmax(F1s[:,i])
    BestF1[mxi] += 1
end
@assert(sum(BestF1) == 10)

# \\f{\$\\pm$(0.01)\$}

mx, bestf1 = findmax(f1_m)
mx, bestrun = findmin(runs_m)
mx, countf1 = findmax(BestF1)

f1_m = round.(f1_m, digits = 3)
runs_m = round.(runs_m, digits = 1)
f1_std = round.(f1_std, digits = 2)
runs_std = round.(runs_std, digits = 1)

clusname = names[p[i]]
println("")
methods = ["\$\\delta\$-linear" "clique" "sqrt" "\$x^{0.9}\$"]
for r = 1:4
    method = methods[r]
    if r == 1
        write(io,"\\midrule $clusname & $(Tsizes[p[i]]) & $method &")
    else
        write(io,"& & $method &")
    end
    if r == bestf1
        write(io,"\\textbf{$(f1_m[r])}  \\f{\$\\pm$(f1_std[r])\$} & ")
    else
        write(io,"{$(f1_m[r])}  \\f{\$\\pm$(f1_std[r])\$} & ")
    end
    if r == bestrun
        write(io,"\\textbf{$(runs_m[r])} \\f{\$\\pm$(runs_std[r])\$} & ")
    else
        write(io,"$(runs_m[r]) \\f{\$\\pm$(runs_std[r])\$} & ")
    end
    if r == countf1
        write(io,"\\textbf{$(BestF1[r])} \\\\\n")
    else
        write(io,"$(BestF1[r]) \\\\\n")
    end
end
end

write(io,"\\bottomrule \\end{tabular}
\\end{table} \n")


end

end
