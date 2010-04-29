function X = transitive_closures( X )

for d=1:4
    X.reach{d} = transitive_closure( X.contact{d} , find(X.taken_ids) ) ;
end

end