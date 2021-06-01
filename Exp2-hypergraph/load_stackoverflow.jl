# You need to do to the data folder and generate Stackoverflow_data.mat
# before this will work
s = time()
M = matread("../data/Stackoverflow_data.mat")
LabelMatrix = M["LabelMat"]
LabelNames = M["LabelNames"]
H = M["H"]
Tconds = M["ClusConds"]
Tsizes = M["ClusSizes"]
order = round.(Int64,vec(sum(H,dims=2)))
d = vec(sum(H,dims=1))
volA = sum(d)
m,n = size(H)
Ht = sparse(H')
toc = time()-s
println("Done loading things into memory in $toc seconds.")
