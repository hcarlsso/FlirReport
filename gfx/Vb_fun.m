function V=Vb_fun(t,ti)
    if t<ti
        V=3*ones(1,length(t));
    else
        V=0*ones(1,length(t));
    end
end