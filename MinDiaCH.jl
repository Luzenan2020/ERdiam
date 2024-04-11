include("graph.jl")
include("core.jl")
include("bb.jl")
include("meth.jl")


using LinearAlgebra
using Laplacians
using SparseArrays
using PyPlot


fname = open("filename.txt", "r")
str   = readline(fname);
nn     = parse(Int, str);


for nnnn=1:nn 
    

    str = readline(fname);
    str = split(str);
    G   = get_graph(str[1]);
    on=G.n;om=G.m;
    Gc=findconnect(G)
    G=Gc;
    n=G.n;m=G.m;
    
    L=lapsp(G);
    A=adjsp(G);

    disfi=accR(A,n);
    println(disfi);
    timeacc=0;
    t1=time();
    kk=50;

    for i=1:kk
        
        eps=0.3;
        t=Int(24*round(log(n)/eps^2)+1);

        t1=time();
        Q=zeros(t,m); #t*m
        for i in 1:t
            for j in 1:m
                Q[i, j] = rand(0:1) == 1 ? 1 : -1
            end
        end
        QB=Q*B; #t*n
        Z=zeros(t,n);
        f=approxchol_lap(A);
        for i=1:t
            Z[i,:]=f(QB[i,:]);
        end
        new_lab=Z';

        K=100;
        index_this=bb(new_lab,K);
        ll=length(index_this);


        new_lab1=[];
        for i in index_this
            push!(new_lab1,new_lab[i,:]);
        end
        dis3=0;
        xx=0;
        yy=0;
        for i=1:ll
            for j=1:ll
                tmp_dis=sum((new_lab1[i] .- new_lab1[j]).^2)/t;
                if dis3 < tmp_dis
                    dis3=tmp_dis;
                    xx=index_this[i];
                    yy=index_this[j];
                end
            end
        end
        # println(xx,"  ",yy);
        #A[xx,yy]+=1;A[yy,xx]+=1;
        Gv=G.v;
        push!(Gv,xx);
        Gu=G.u;
        push!(Gu,yy);
        G=Graph(G.n, G.m+1, Gv, Gu, G.nbr);
        Gc=findconnect(G)
        G=Gc;
        n=G.n;m=G.m;
        A=adjsp(G);
        B=getB(G);
        if i%10==0
            t3=time();
            dis2=accR(A,n);
            println(dis2);
            t4=time();
            timeacc+=t4-t3;
        end
        
    end

    
    
    t2=time();
    println(t2-t1-timeacc);
    
    dis2=accR(A,n);
    println(dis2);



end
close(fname)
