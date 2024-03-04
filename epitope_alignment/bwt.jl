module BWTs

export BWT, findKmer

function fillBWTCO!(bwt::Array{UInt8},
                    c::Dict{UInt8,Int64},
                    o::Dict{UInt8,Array{Int64}},
                    seq::Array{UInt8},
                    sa::Array{Int64},
                    BLOCK_SIZE::Int64)::Nothing
    oSize = ((length(bwt)-1) ÷ BLOCK_SIZE) + 1
    @inbounds for i=1:length(bwt)
        if sa[i] > 1
            ichar = seq[ sa[i]-1 ]
        else
            ichar = b"!"[1]
        end
        bwt[i] = ichar
        if !(ichar in keys(c))
            c[ichar] = 0
        end
        c[ichar] += 1
        if !(ichar in keys(o))
            o[ichar] = zeros(Int64, oSize)
        end
        o[ichar][1 + ((i-1) ÷ BLOCK_SIZE)] += 1
    end
    cumul = 0
    cKeys = sort(collect(keys(c)))
    @inbounds for ichar in cKeys
        cumul += c[ichar]
        c[ichar] = cumul - c[ichar]
        ochar = o[ichar]
        for i=2:oSize
            ochar[i] += ochar[i-1]
        end
    end
    return nothing
end

mutable struct BWT
    seq::Array{UInt8}
    sa::Array{Int64}
    bwt::Array{UInt8}
    c::Dict{UInt8,Int64}
    o::Dict{UInt8,Array{Int64}}

    BLOCK_SIZE::Int64

    function BWT(seq::Array{UInt8}, sa::Array{Int64})
        this = new()
        this.seq = seq
        this.sa = sa
        this.bwt = Array{UInt8}(undef, length(this.seq))
        this.c = Dict{UInt8,Int64}()
        this.o = Dict{UInt8,Int64}()
        this.BLOCK_SIZE = 100
        fillBWTCO!(this.bwt, this.c, this.o, this.seq, this.sa, this.BLOCK_SIZE)
        return this
    end
end

function handleBlocking(bwt::Array{UInt8},
                        c::Dict{UInt8,Int64},
                        o::Dict{UInt8,Array{Int64}},
                        i::Int64,
                        ichar::UInt8,
                        blockSize::Int64)::Int64
    iRounded = ((i-1) ÷ blockSize)
    iNew = c[ichar] + (iRounded > 0 ? o[ichar][iRounded] : 0) + 1
    for j=(1+(blockSize*iRounded)):(i-1)
        if bwt[j] == ichar
            iNew += 1
        end
    end
    return iNew
end

function findKmer(this::BWT, kmer::Array{UInt8})::Tuple{Int64,Int64}
    k = length(kmer)
    start = 1
    end_ = length(this.bwt) + 1
    for i=k:-1:1
        ichar = kmer[i]
        start = handleBlocking(
            this.bwt, this.c, this.o, start, ichar, this.BLOCK_SIZE
        )
        end_ = handleBlocking(
            this.bwt, this.c, this.o, end_, ichar, this.BLOCK_SIZE
        )
    end
    return (start, end_)
end

function calculateD(rev::BWT, kmer::Array{UInt8})::Array{Int64}
    k = length(kmer)
    d = zeros(k)
    start = 1
    end_ = length(rev.bwt) + 1
    z = 0
    for i=1:k
        ichar = kmer[i]
        start = handleBlocking(rev.bwt, rev.c, rev.o, start, ichar, rev.BLOCK_SIZE)
        end_ = handleBlocking(rev.bwt, rev.c, rev.o, end_, ichar, rev.BLOCK_SIZE)
        if end_ <= start
            start = 1
            end_ = length(rev.bwt) + 1
            z += 1
        end
        d[i] = z
    end
    return d
end

function findKmer(this::BWT, kmer::Array{UInt8},
                  nMismatch::Int64, d::Array{Int64},
                  i::Int64, start::Int64, end_::Int64;
                  indels::Bool=false)::Array{Tuple{Int64,Int64}}
    if i < 1
        if nMismatch >= 0
            return [(start, end_)]
        else
            return Array{Tuple{Int64,Int64}}(undef, 0)
        end
    end    
    if nMismatch < d[i]
        return Array{Tuple{Int64,Int64}}(undef, 0)
    end
    I = Array{Tuple{Int64,Int64}}(undef, 0)
    if indels
        append!(I, findKmer(
            this, kmer, nMismatch-1, d, i-1, start, end_, indels=true
        ))
    end
    ichar = kmer[i]
    for bchar in sort(collect(keys(this.c)))
        bstart = handleBlocking(this.bwt, this.c, this.o, start, bchar, this.BLOCK_SIZE)
        bend = handleBlocking(this.bwt, this.c, this.o, end_, bchar, this.BLOCK_SIZE)
        if bend > bstart
            if indels
                append!(I, findKmer(
                    this, kmer, nMismatch-1, d, i, bstart, bend, indels=true
                ))
            end
            if bchar == ichar
                append!(I, findKmer(
                    this, kmer, nMismatch, d, i-1, bstart, bend, indels=indels
                ))
            else
                append!(I, findKmer(
                    this, kmer, nMismatch-1, d, i-1, bstart, bend, indels=indels
                ))
            end
        end
    end
    return I
end

function findKmer(this::BWT, rev::BWT, kmer::Array{UInt8},
                  nMismatch::Int64; indels::Bool=false)::Array{Tuple{Int64,Int64}}
    d = calculateD(rev, kmer)
    return findKmer(this, kmer, nMismatch, d,
                    length(kmer), 1, length(this.bwt)+1, indels=indels)
end

end # module
