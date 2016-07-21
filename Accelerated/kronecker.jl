using Base.Threads

function kronecker(scale::Int, edge_factor::Real)
    N = 1 << scale                
    M = Int(edge_factor * N)

    A, B, C = 0.57, 0.19, 0.19    

    # vertex arrays (edge list)
    v1 = ones(Int, M)
    v2 = ones(Int, M)

    ab = A + B
    c_norm = C/(1-(A+B))
    a_norm = A/(A+B)

    @threads for j = 1:M
        for ib = 1:scale          
            k = 1 << (ib-1)
            if rand() > ab
                v1[j] += k
                if rand() > c_norm
                    v2[j] += k
                end
            elseif rand() > a_norm
                v2[j] += k
            end
        end
    end

    v1, v2
end

