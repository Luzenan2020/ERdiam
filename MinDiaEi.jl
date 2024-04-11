include("graph.jl")
include("core.jl")
include("acc.jl")

using LinearAlgebra
using Laplacians
using SparseArrays
using Arpack
using Random
using LightGraphs


fname = open("filename.txt", "r")
str   = readline(fname);
nn    = parse(Int, str);


for ttttt=1:nn
    str = readline(fname);
    str = split(str);
    G   = get_graph(str[1]);
    on=G.n;om=G.m;
    Gc=findconnect(G)
    G=Gc;
    n=G.n;m=G.m;
    L=lapsp(G);
    A=adjsp(G);
    B=getB(G);
    pairx=0;
    pairy=0;

    t1=time();
    println(accR(A,n))
    kk=50;
    for i=1:kk
        A=adjsp(G);
        f = approxchol_lap(A);
        op = Laplacians.SqLinOp(true,1.0,size(A,1),f);
        e = eigs(op, which=:LM, nev=1);

        ei=e[2];
        #aei=eigs(L, nev = 2, which=:SM);

        pairx=argmax(ei)[1];
        pairy=argmin(ei)[1];
        #A[pairx,pairy]+=1;A[pairy,pairx]+=1;
        #disfi1=accR(A,n);
        # println(disfi1);
        #println(sum(A))

        Gv=G.v;
        push!(Gv,pairx);
        Gu=G.u;
        push!(Gu,pairy);
        G=Graph(G.n, G.m+1, Gv, Gu,G.nbr);
        n=G.n;m=G.m;
        A=adjsp(G);
        disfi2=accR(A,n);
        println(disfi2);
    end


    t2=time();

    




end

