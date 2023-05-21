function read_2d!(data::Array{T,2}, filepath::String) where {T<:Number}
    if prod(size(data))*sizeof(T) != stat(filepath).size
        throw(DimensionMismatch("Dimensions or datatype do not match file on disk"))
    end
    open(filepath,"r") do file
        read!(file,data)
    end
end

function read_bruker_p2d(filepath, f1dim, f2dim, f1length; stride_len = 1)
    rawdata = Array{Int32,2}(undef, f2dim, f1dim)
    read_2d!(rawdata,filepath)
    max = maximum(rawdata)/10.0
    s = copy(rawdata'[1:f1length,:]./max)

    stride_len = clamp(stride_len,1,f2dim)
    if stride_len == 1
        return s
    end

    s_avg = zeros(size(s))
    if iseven(stride_len)
        stride_len = stride_len + 1
    end
    for i in (stride_len÷2)+1:f2dim-(stride_len÷2)
        s_avg[:,i] .= sum(s[:,i-(stride_len÷2):i+(stride_len÷2)],dims=2)./stride_len
    end
    return s_avg
end