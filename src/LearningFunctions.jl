function create_train_test_split(X,y,train_test_split = 0.8)

    y_sample_id = [shuffle(findall(x->x==label,y)) for label in unique(y)]

    n_train_id = [Int(floor(train_test_split*length(y_sid))) for y_sid in y_sample_id]

    y_train_id = [y_sid[1:n_train] for (y_sid,n_train) in zip(y_sample_id,n_train_id)]

    y_test_id = [y_sid[n_train+1:end] for (y_sid,n_train) in zip(y_sample_id,n_train_id)]

    train_id = reduce(vcat,y_train_id )
    test_id = reduce(vcat,y_test_id )

    X_train = X[train_id,:]
    X_test = X[test_id,:]

    y_train = y[train_id]
    y_test = y[test_id]

    return X_train,X_test,y_train,y_test, train_id, test_id

end

function create_train_test_split_multilabel(X,y,train_test_split = 0.8)

    y_sample_id = [shuffle(findall(mapslices(row->row==label,y,dims = 2)[:,1])) for label in unique(eachrow(y))];

    n_train_id = [Int(floor(train_test_split*length(y_sid))) for y_sid in y_sample_id]

    y_train_id = [y_sid[1:n_train] for (y_sid,n_train) in zip(y_sample_id,n_train_id)]

    y_test_id = [y_sid[n_train+1:end] for (y_sid,n_train) in zip(y_sample_id,n_train_id)]

    train_id = reduce(vcat,y_train_id )
    test_id = reduce(vcat,y_test_id )

    X_train = X[train_id,:]
    X_test = X[test_id,:]

    y_train = y[train_id,:]
    y_test = y[test_id,:]

    return X_train,X_test,y_train,y_test, train_id, test_id

end

function create_train_test_id_split(y,train_test_split = 0.8)

    y_sample_id = [shuffle(findall(x->x==label,y)) for label in unique(y)]

    n_train_id = [Int(floor(train_test_split*length(y_sid))) for y_sid in y_sample_id]

    y_train_id = [y_sid[1:n_train] for (y_sid,n_train) in zip(y_sample_id,n_train_id)]

    y_test_id = [y_sid[n_train+1:end] for (y_sid,n_train) in zip(y_sample_id,n_train_id)]

    train_id = reduce(vcat,y_train_id )
    test_id = reduce(vcat,y_test_id )

    return train_id, test_id

end

function create_train_test_id_split_multilabel(y,train_test_split = 0.8)

    y_sample_id = [shuffle(findall(mapslices(row->row==label,y,dims = 2)[:,1])) for label in unique(eachrow(y))];

    n_train_id = [Int(floor(train_test_split*length(y_sid))) for y_sid in y_sample_id]

    y_train_id = [y_sid[1:n_train] for (y_sid,n_train) in zip(y_sample_id,n_train_id)]

    y_test_id = [y_sid[n_train+1:end] for (y_sid,n_train) in zip(y_sample_id,n_train_id)]

    train_id = reduce(vcat,y_train_id )
    test_id = reduce(vcat,y_test_id )

    return train_id, test_id

end
