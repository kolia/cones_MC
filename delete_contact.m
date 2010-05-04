function X = delete_contact( X , id , d )

there                               = X.contact{d}(id,:) ;
X.contact{d}(id,there)              = false ;
X.contact{1+mod(d+1,4)}(there,id)   = false ;
X.contact{d}(id,id)                 = false ;

end