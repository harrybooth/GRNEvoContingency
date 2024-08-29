function brier_score(prob,label)
    total = 0
    for (n,p) in enumerate(prob)
        if n == label
            total += (p-1)^2
        else
            total += p^2
        end
    end

    return (1/length(prob))*total
end

function positive_label_prop(y,bin_id,bins,label)
    d = y[findall(x->x==bin_id,bins)]
    count(x->x==label,d) /length(d)
end

function positive_label_prop_binomial(y,bin_id,bins,label)
    d = y[findall(x->x==bin_id,bins)]

    n = length(d)
    x = count(x->x==label,d)

    return x/n, confint(BinomialTest(x, n), level = 0.95)
end

function produce_prob_calibration_curve!(fig,y_train_gtl,gtl_prob_train,y_test_gtl,gtl_prob_test,predict_colors)

    ax_train = Axis(fig[1,1],title = "Train (uncalibrated)")
    ax_test = Axis(fig[2,1],title = "Test (uncalibrated)")

    hist_edges = 0:0.1:1
    hist_edges_middle = 0.05:0.1:1.

    n_class = top_n + 1

    for label in 1:n_class

        all_cal = []
        all_mean = []

        h_fitness = fit(Histogram, gtl_prob_train[:,label], hist_edges; closed = :left) 

        pred_prob_bins = map(f->StatsBase.binindex(h_fitness, f),gtl_prob_train[:,label])

        for i in 1:10
            push!(all_cal,positive_label_prop(y_train_gtl,i,pred_prob_bins,label))
            push!(all_mean,mean(gtl_prob_train[pred_prob_bins .== i,label]))
        end

        CairoMakie.scatter!(ax_train,Float64.(all_mean),Float64.(all_cal), label = "Label " * string(label),color = predict_colors[label], markersize =15.,marker = 'x')
        CairoMakie.lines!(ax_train,Float64.(all_mean),Float64.(all_cal), label = "Label " * string(label),color = predict_colors[label], linewidth = 2.)

        all_cal = []
        all_mean = []

        h_fitness = fit(Histogram, gtl_prob_test[:,label], hist_edges; closed = :left) 

        pred_prob_bins = map(f->StatsBase.binindex(h_fitness, f),gtl_prob_test[:,label])

        for i in 1:10
            push!(all_cal,positive_label_prop(y_test_gtl,i,pred_prob_bins,label))
            push!(all_mean,mean(gtl_prob_test[pred_prob_bins .== i,label]))
        end
        
        CairoMakie.scatter!(ax_test,Float64.(all_mean),Float64.(all_cal), label = "Label " * string(label),color = predict_colors[label] , markersize =15.,marker = 'x')
        CairoMakie.lines!(ax_test,Float64.(all_mean),Float64.(all_cal), label = "Label " * string(label),color = predict_colors[label], linewidth = 2.)

    end

    CairoMakie.lines!(ax_train,hist_edges_middle |> collect,map(x->x,hist_edges_middle), label = "Perfect calibration", linestyle = :dash, color = :black)
    CairoMakie.lines!(ax_test,hist_edges_middle |> collect,map(x->x,hist_edges_middle), label = "Perfect calibration", linestyle = :dash, color = :black)

    ax_train.xticks = hist_edges_middle
    ax_test.xticks = hist_edges_middle
end

function produce_prob_calibration_curve_binomial!(fig,y_train_gtl,gtl_prob_train,y_test_gtl,gtl_prob_test,predict_colors)

    ax_train = Axis(fig[1,1],xlabel = L"\text{Predicted probabilities}",ylabel = L"\text{Realised frequencies}",title = L"\text{Train}")
    ax_test = Axis(fig[2,1],xlabel = L"\text{Predicted probabilities}",ylabel = L"\text{Realised frequencies}",title = L"\text{Test}")

    hist_edges = 0:0.1:1
    hist_edges_middle = 0.05:0.1:1.

    n_class = top_n + 1

    for label in 1:n_class

        all_cal = []
        all_cal_int = []
        all_mean = []

        h_fitness = fit(Histogram, gtl_prob_train[:,label], hist_edges; closed = :left) 

        pred_prob_bins = map(f->StatsBase.binindex(h_fitness, f),gtl_prob_train[:,label])

        for i in 1:10
            cal_val,cal_int = positive_label_prop_binomial(y_train_gtl,i,pred_prob_bins,label)
            push!(all_cal,cal_val)
            push!(all_cal_int,cal_int)
            push!(all_mean,mean(gtl_prob_train[pred_prob_bins .== i,label]))
        end

        CairoMakie.scatter!(ax_train,Float64.(all_mean),Float64.(all_cal), label = "Label " * string(label),color = predict_colors[label], markersize =15.,marker = 'x')
        CairoMakie.lines!(ax_train,Float64.(all_mean),Float64.(all_cal), label = "Label " * string(label),color = predict_colors[label], linewidth = 2.)
        CairoMakie.rangebars!(ax_train,Float64.(all_mean), first.(all_cal_int),  last.(all_cal_int), color = predict_colors[label],whiskerwidth = 10, direction = :y)

        all_cal = []
        all_cal_int = []
        all_mean = []

        h_fitness = fit(Histogram, gtl_prob_test[:,label], hist_edges; closed = :left) 

        pred_prob_bins = map(f->StatsBase.binindex(h_fitness, f),gtl_prob_test[:,label])

        for i in 1:10
            cal_val,cal_int = positive_label_prop_binomial(y_test_gtl,i,pred_prob_bins,label)
            push!(all_cal,cal_val)
            push!(all_cal_int,cal_int)
            push!(all_mean,mean(gtl_prob_test[pred_prob_bins .== i,label]))
        end
        
        CairoMakie.scatter!(ax_test,Float64.(all_mean),Float64.(all_cal), label = "Label " * string(label),color = predict_colors[label] , markersize =15.,marker = 'x')
        CairoMakie.lines!(ax_test,Float64.(all_mean),Float64.(all_cal), label = "Label " * string(label),color = predict_colors[label], linewidth = 2.)
        CairoMakie.rangebars!(ax_test,Float64.(all_mean), first.(all_cal_int),  last.(all_cal_int), color = predict_colors[label],whiskerwidth = 10, direction = :y)

    end

    CairoMakie.lines!(ax_train,hist_edges_middle |> collect,map(x->x,hist_edges_middle), label = "Perfect calibration", linestyle = :dash, color = :black)
    CairoMakie.lines!(ax_test,hist_edges_middle |> collect,map(x->x,hist_edges_middle), label = "Perfect calibration", linestyle = :dash, color = :black)

    ax_train.xticks = hist_edges_middle
    ax_test.xticks = hist_edges_middle
end

function produce_prob_calibration_curve_binomial_test!(fig,y_train_gtl,gtl_prob_train,y_test_gtl,gtl_prob_test,predict_colors)

    ax_test = Axis(fig[1,1],xlabel = L"\text{Predicted probabilities}",ylabel = L"\text{Realised frequencies}")

    hist_edges = 0:0.1:1
    hist_edges_middle = 0.05:0.1:1.

    n_class = top_n + 1

    for label in 1:n_class

        all_cal = []
        all_cal_int = []
        all_mean = []

        h_fitness = fit(Histogram, gtl_prob_train[:,label], hist_edges; closed = :left) 

        pred_prob_bins = map(f->StatsBase.binindex(h_fitness, f),gtl_prob_train[:,label])

        for i in 1:10
            cal_val,cal_int = positive_label_prop_binomial(y_train_gtl,i,pred_prob_bins,label)
            push!(all_cal,cal_val)
            push!(all_cal_int,cal_int)
            push!(all_mean,mean(gtl_prob_train[pred_prob_bins .== i,label]))
        end

        all_cal = []
        all_cal_int = []
        all_mean = []

        h_fitness = fit(Histogram, gtl_prob_test[:,label], hist_edges; closed = :left) 

        pred_prob_bins = map(f->StatsBase.binindex(h_fitness, f),gtl_prob_test[:,label])

        for i in 1:10
            cal_val,cal_int = positive_label_prop_binomial(y_test_gtl,i,pred_prob_bins,label)
            push!(all_cal,cal_val)
            push!(all_cal_int,cal_int)
            push!(all_mean,mean(gtl_prob_test[pred_prob_bins .== i,label]))
        end
        
        CairoMakie.scatter!(ax_test,Float64.(all_mean),Float64.(all_cal), label = "Label " * string(label),color = predict_colors[label] , markersize =15.,marker = 'x')
        CairoMakie.lines!(ax_test,Float64.(all_mean),Float64.(all_cal), label = "Label " * string(label),color = predict_colors[label], linewidth = 2.)
        CairoMakie.rangebars!(ax_test,Float64.(all_mean), first.(all_cal_int),  last.(all_cal_int), color = predict_colors[label],whiskerwidth = 10, direction = :y)

    end

    CairoMakie.lines!(ax_test,hist_edges_middle |> collect,map(x->x,hist_edges_middle), label = "Perfect calibration", linestyle = :dash, color = :black)

    ax_test.xticks = hist_edges_middle
end

function find_streak_gt(tr,m_type)

    pred_error = tr.gt_prediction_error[:,1]
    S0 = tr.H0

    if m_type == :fitness
        m = tr.fitness_traj
    else
        m = tr.weight_edits
    end

    l_id = filter(x->x!=1,sort(findall(p->p==1,pred_error[1:S0-1])))

    if length(l_id) > 0
        contender = []
        for l in l_id
            streak_length = 0 
            for p in pred_error[l:S0-1]
                if p == 1
                    streak_length+=1
                end
            end
            if streak_length == length(l:S0-1)
                push!(contender,m[l]/m[S0])
            end
        end

        if length(contender) > 0
            return contender[1]
        else
            return 1
        end
    else
        return 1
    end
end

function find_streak_tt(tr,m_type)

    pred_error = tr.tt_prediction_error[:,1]
    S0 = tr.H0

    if m_type == :fitness
        m = tr.fitness_traj
    else
        m = tr.weight_edits
    end

    l_id = filter(x->x!=1,sort(findall(p->p==1,pred_error[1:S0-1])))

    if length(l_id) > 0
        contender = []
        for l in l_id
            streak_length = 0 
            for p in pred_error[l:S0-1]
                if p == 1
                    streak_length+=1
                end
            end
            if streak_length == length(l:S0-1)
                push!(contender,m[l]/m[S0])
            end
        end

        if length(contender) > 0
            return contender[1]
        else
            return 1
        end
    else
        return 1
    end
end

function find_streak_hamming(tr)
    pred_error = tr.gt_prediction_error[:,1]
    S0 = tr.H0

    l_id = filter(x->x!=1,sort(findall(p->p==1,pred_error[1:S0-1])))

    if length(l_id) > 0
        contender = []
        for l in l_id
            streak_length = 0 
            for p in pred_error[l:S0-1]
                if p == 1
                    streak_length+=1
                end
            end
            if streak_length == length(l:S0-1)
                push!(contender,masked_hamming_distance(tr.geno_traj[l],tr.minimal_stripe_subgraphs[S0]))
            end
        end

        if length(contender) > 0
            return contender[1]
        else
            return -1
        end
    else
        return -1
    end
end

function calculate_pred_proportion(preds,top_n)
    
    total = length(preds)

    return [count(x->x==p,preds)/total for p in 1:top_n]
end