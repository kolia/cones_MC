function bestX = load_bestX( dir , pattern )

cwd = pwd ;
cd(dir)

lsd = ls(pattern) ;
filenames = {} ;

if ~strcmp(lsd,'ls: ')  % files exist
    filenames = regexp(lsd,'\s','split') ;
end

bestX  = {} ;
for i=1:length(filenames)
    if ~isempty(filenames{i})
        s = load(filenames{i}) ;
        bestX = [bestX ; s.bestX(:)] ;
    end
end

cd(cwd)

end