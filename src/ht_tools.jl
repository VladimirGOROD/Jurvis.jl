function fullsignal(vec::AbstractVector{<:Number})
    return dropdims(Hilbert.hilbert(reshape(vec, (length(vec),1))), dims = 2);
end

function envelope(vec::AbstractVector{<:Real})
    return abs.(fullsignal(vec))
end

function envelope(vec::AbstractVector{Complex})
    return abs.(vec);
end

function unwrap(vec::AbstractVector{<:Number})
    unwrapped = copy(vec)
    for i in 2:length(v)
        while unwrapped[i] - unwrapped[i-1] >= pi
        unwrapped[i] -= 2pi
        end
        while unwrapped[i] - unwrapped[i-1] <= -pi
        unwrapped[i] += 2pi
        end
    end
    return unwrapped
end
function unwrap!(vec::AbstractVector{<:Number})
    
    for i in 2:length(vec)
        while vec[i] - vec[i-1] >= pi
            vec[i] -= 2pi
        end
        while vec[i] - vec[i-1] <= -pi
            vec[i] += 2pi
        end
    end
    return vec
end

function instphase(vec::AbstractVector{<:Real})
    return unwrap!(angle.(fullsignal(vec)))    
end

function instphase(vec::AbstractVector{Complex})
    return unwrap!(angle.(vec))    
end



function instdamping(ω, A, dt)
    Ȧ=diff14(A,dt);
    dω=diff14(ω,dt);
    β = similar(ω);
    for i in eachindex(ω);
        β[i]=-Ȧ[i]/A[i]-(dω[i])/(2*ω[i]);
    end
    return β
end

function instfrequency(data::MeasuredData, channel)
    
end