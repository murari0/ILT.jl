function read_bruker_p2d(filename, f1dim, f2dim, f1length; stride_len = 1)
    rawdata = Vector{Int32}(undef, f1dim*f2dim)
    open(filename,"r") do file
        read!(file,rawdata)
    end
    max = maximum(rawdata)/10.0

    s = copy(reshape(rawdata[1:f1length*f2dim]./max,f1length,f2dim)')
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