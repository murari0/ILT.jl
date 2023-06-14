"""
    read_2d!(data, filepath)

Reads data stored in the binary file at `filepath` into a 2-dimentional matrix `data`.
"""
function read_2d!(data::Array{T,2}, filepath::String) where {T<:Number}
    if prod(size(data))*sizeof(T) != stat(filepath).size
        throw(DimensionMismatch("Dimensions or datatype do not match file on disk"))
    end
    open(filepath,"r") do file
        read!(file,data)
    end
end

"""
    read_bruker_p2d(filepath, f1dim, f2dim, f1length; stride_len = 1)

Reads a binary file written by Bruker's TopSpin software located at `filepath` and returns the data in a 2-dimensional matrix.
    
Utility function to read a Fourier transformed and phased pseudo 2-dimensional spectrum saved by Bruker's TopSpin software to the file at `filepath` (generally a 2rr file). The length in the f2 dimension (frequency domain) is given by `f2dim`, the length (processed) in the f1 dimension is `f1dim` and the actual number of points acquired in the f1 dimension is `f1length`. Optionally, the keyword argument `stride_len` instructs the function to smooth the data in the f2 dimension using a moving average with a window length of `stride_len`.
"""
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
