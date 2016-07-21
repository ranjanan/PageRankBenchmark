
function writeuv(fname1, fname2, u, v)
    
    n = length(u)

    f1 = open(fname1, "w")
    f2 = open(fname2, "w")
    @threads for i in [true, false]
        for j = 1:n 
            i ? write(f1, u[j]) : write(f2, v[j])
        end
    end
end
