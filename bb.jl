include("cc.jl")

using LightGraphs
using SimpleWeightedGraphs
using LinearAlgebra
using Random

function bb(mat_A,K)
    # t1=time();
    K=K;
    MN=size(mat_A);
    NN=MN[1];
    MM=MN[2];
    epsilon=0.8;
    eudis=sum(mat_A.*mat_A,dims=2);
    max_index=argmax(eudis)[1];
    index_set=zeros(Int,1,NN);
    index_set[max_index]=1;
    #println(eudis[max_index]);
    candidate_set=ones(Int,1,NN);
    candidate_set[max_index]=0;
    index_series2=union(1:NN);
    setdiff!(index_series2,max_index);
    alpha_zero=ones(sum(index_set),1);
    pre_index2=[];
    index_set2=[];
    candidate_set2=[];
    push!(index_set2,max_index);
    push!(candidate_set2,max_index);
    push!(pre_index2,max_index);
    # mustindex=[];
    # for i=1:MM
    #     amax=argmax(mat_A[:,i]);
    #     amin=argmin(mat_A[:,i]);
    #     push!(mustindex,amax);
    #     push!(mustindex,amin);
    # end
    # mustindex=union(mustindex)
    # println(mustindex);
    tmp_set=Int[];
    all_set=union(1:NN)
    selc=zeros(Int,NN);
    selc[max_index]=1;
    tmp_set=union(1:NN);
    setdiff!(tmp_set,max_index)
    tmp_add=Int[];
    tot=0;
    #t2=time();
    # tt1=t2-t1;
    # print("tt1");
    # println(tt1);
    # while sum(index_set)<K
    #     candidate_set=ones(Int,1,NN);
    #    for i in 1:size(index_set)[2]
    #        if index_set[1,i]==1
    #            candidate_set[1,i]=0;
    #            push!(candidate_set2,i);
    #        end
    #    end
        while sum(candidate_set)>0
            # println(sum(candidate_set));
            # println(length(index_set2));

            #t3=time();
            saz=size(alpha_zero);
            azN=saz[1];
            azM=saz[2];

            if length(index_set2)!=azN && length(pre_index2)!=0
                pre_series=pre_index2;
                now_index_set=index_set2;
                tmp_idx=sortperm(now_index_set);
                lennow=length(now_index_set);
                alpha_zero2=zeros(lennow,azM)
                for i in 1:lennow
                    if tmp_idx[i]!=lennow
                        alpha_zero2[i,:]=alpha_zero[tmp_idx[i],:];
                    else
                        alpha_zero2[i,:]=alpha_zero2[i,:];
                    end
                end
                alpha_zero=alpha_zero2;
                pre_index2=ones(Int,1,length(index_set2));
                for i in 1:length(index_set2)
                    pre_index2[i]=index_set2[i];
                end
            end
            # len=length(index_series2);

            # rnd_ind=rand(1:len,1,1);
            # this_index=index_series2[rnd_ind[1]];
            if tot%1000==0
                tmp_set=setdiff(tmp_set,tmp_add)
                tmp_add=Int[];
                #println(tot)
            end
            this_index=rand(tmp_set)
            while selc[this_index]==1
                this_index=rand(tmp_set)
            end


            this_data=mat_A[index_set2,:];
            p=mat_A[this_index,:];
            mn=size(mat_A);
            pp=zeros(1,mn[2]);#（1，100）
            for i in 1:mn[2]
                pp[1,i]=p[i];
            end
            alpha_zeronew=alpha_zero[:,1];
            #t4=time();
            # tt2=t3-t4;
            # print("tt2");
            # print(tt2);
            
            inorout,p_prime,alpha_coe,dist=cc(this_data,pp,epsilon,alpha_zeronew);
            #println(sqrt((pp'-p_prime)'*(pp'-p_prime))[1,1]);
            #t5=time();
            # tt3=t5-t4;
            # print("tt3");
            # println(tt3);
            if inorout==1
                #ta1=time();
                candidate_set[1,this_index]=0;
                push!(candidate_set2,this_index);
                push!(tmp_add,this_index);
                tot+=1;
                selc[this_index]=1;
                # setdiff!(index_series2,this_index);
                #ta2=time();
                #tta1=ta2-ta1;
                # print("tta1");
                # println(tta1);
            else
                direction=p_prime.-pp;
                index_series2=setdiff(index_series2,candidate_set2)
                S_data=mat_A[index_series2,:];
                # println(size(S_data));
                projected_val=S_data*direction;
                index_2=argmin(projected_val);
                index_candidate_2=index_series2[index_2[1]];
                if index_set[index_candidate_2]==0
                    index_set[index_candidate_2]=1;
                    push!(index_set2,index_candidate_2);
                    if length(index_set2)>=K
                       break;
                   end
                    candidate_set[index_candidate_2]=0;
                    push!(candidate_set2,index_candidate_2);
                    push!(tmp_add,index_candidate_2);
                    selc[index_candidate_2]=1;
                    tot+=1;
                    # setdiff!(index_series2,index_candidate_2);
                end
            end
            # t6=time();
            # tt4=t6-t5;
            # print("tt4");
            # println(tt4);
        end
    # end
    return index_set2;
end
