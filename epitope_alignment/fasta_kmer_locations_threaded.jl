#!/usr/bin/env julia

include("bwt.jl")
include("skewsuffixarray.jl")

using DataFrames
using DataStructures
using GZip

using Main.SkewSuffixArray
using Main.BWTs

function readFasta(filename::String)::OrderedDict{String,String}
    out = OrderedDict{String,String}()
    curName = ""
    sb = IOBuffer()
    fgz = (endswith(lowercase(filename), ".gz") ||
           endswith(lowercase(filename), ".gzip"))
    fh = (fgz ? GZip.open(filename, "r") : open(filename, "r"));
    while !eof(fh)
        curLine = rstrip(readline(fh))
        if startswith(curLine, ">")
            if curName != ""
                out[curName] = String(take!(sb))
            end
            curName = curLine[2:length(curLine)]
        else
            write(sb, curLine)
        end
    end
    if curName != "" && !(curName in keys(out))
        out[curName] = String(take!(sb))
    end
    close(fh)
    close(sb)
    return out
end

function constructSuffixArray(str::Array{UInt8})
    chars = convert.(Int64, [str[i] for i in 1:length(str)])
    append!(chars, [0, 0, 0])
    SA = zeros(Int64, length(chars))
    suffixArray!(chars, SA, length(str), maximum(chars))
    return SA[1:length(str)] .+ 1
end

mutable struct KmerLocator
    refSeqs::OrderedDict{String,String}
    seqNames::Array{String}
    catSeq::Array{UInt8}
    bounds::Array{Int64}
    sa::Array{Int64}
    bwt::BWT
    revBwt::BWT

    function KmerLocator(fasta::String)
        this = new()
        this.refSeqs = readFasta(fasta)
        this.seqNames = collect(keys(this.refSeqs))
        sb = IOBuffer()
        this.bounds = Array{Int64}(undef, length(this.seqNames)+1)
        this.bounds[1] = 1
        i = 2
        for seq in values(this.refSeqs)
            write(sb, seq)
            write(sb, "\$")
            this.bounds[i] = this.bounds[i-1] + length(seq) + 1
            i += 1
        end
        this.catSeq = take!(sb)
        close(sb)
        this.sa = constructSuffixArray(this.catSeq)
        this.bwt = BWT(this.catSeq, this.sa)
        qesTac = reverse(this.catSeq)
        qesTac = qesTac[2:length(qesTac)]
        push!(qesTac, '\$')
        revSA = constructSuffixArray(qesTac)
        this.revBwt = BWT(qesTac, revSA)
        return this
    end
end

function sourceBlock(this::KmerLocator, s::Int64,
                   mindex::Int64=1, maxdex::Int64=-1)::Int64
    if maxdex == -1
        maxdex = length(this.bounds) + 1
    elseif maxdex == (mindex + 1)
        return mindex
    end
    middex = (mindex + maxdex) ÷ 2
    if s < this.bounds[middex]
        return sourceBlock(this, s, mindex, middex)
    else
        return sourceBlock(this, s, middex, maxdex)
    end
end
function sourceBlock(this::KmerLocator, s0::Array{Int64})::Array{Int64}
    maxdex = length(this.bounds) + 1
    return [sourceBlock(this, s, 1, maxdex) for s in s0]
end

function findKmerInner(this::KmerLocator, kmer::Array{UInt8},
                       nMismatch::Int64=0; indels::Bool=false)::DataFrame
    saIntervals = Main.BWTs.findKmer(this.bwt, this.revBwt, kmer,
                                     nMismatch, indels=indels)
    # saIntervals = [Main.BWTs.findKmer(this.bwt, kmer)]
    out = Array{DataFrame}(undef, 0)
    push!(out, DataFrame(query=Array{String}(undef, 0),
                         seqname=Array{String}(undef, 0),
                         start=Array{Int64}(undef, 0)))
    kms = String(kmer)
    for saInterval in saIntervals
        seqname = Array{String}(undef, saInterval[2]-saInterval[1])
        start = Array{Int64}(undef, saInterval[2]-saInterval[1])
        for i = saInterval[1]:(saInterval[2]-1)
            s = this.sa[i]
            block = sourceBlock(this, s)
            seqname[i-saInterval[1]+1] = this.seqNames[block]
            start[i-saInterval[1]+1] = s - this.bounds[block] + 1
        end
        if length(seqname) > 0
            push!(out, DataFrame(query=kms, seqname=seqname, start=start))
        else
            push!(out, DataFrame(query=Array{String}(undef, 0),
                                 seqname=seqname, start=start))
        end
    end
    catted = vcat(out...)
    catted.target = Array{String}(undef, size(catted, 1))
    overlaps = falses(size(catted, 1))
    for i = 1:size(catted, 1)
        ist = catted[i, :start]
        ien = ist + length(catted[i, :query]) - 1
        if ien <= length(this.refSeqs[catted[i, :seqname]])
            catted[i, :target] = this.refSeqs[catted[i, :seqname]][ist:ien]
        else
            overlaps[i] = true
        end
    end
    return catted[.!overlaps, :]
end

function findKmer(this::KmerLocator, kmers::Array{String},
                  nMismatch::Int64=0; indels::Bool=false)::DataFrame
    results = Array{DataFrame}(undef, length(kmers))
    Threads.@threads for i in 1:length(kmers)
        results[i] = findKmerInner(this, Array{UInt8}(kmers[i]),
                                   nMismatch, indels=indels)
    end
    return vcat(results...)
end


## =============================================================================
# ARGS_ = ["UP000005640_9606.fasta.gz", "UP000000354_694009.fasta.gz", "10", "1"]

kmLocator = KmerLocator(ARGS[1]);

queryFasta = readFasta(ARGS[2]);

k = parse(Int64, ARGS[3]);
kmers = Set{String}();
for seq in values(queryFasta)
    for i = 1:(length(seq)-k+1)
        push!(kmers, seq[i:(i+k-1)])
    end
end

nMismatch = 0
if length(ARGS) > 3
    nMismatch = parse(Int64, ARGS[4])
end

locs = findKmer(kmLocator, collect(kmers), nMismatch);
println("query\tseqname\tstart\ttarget")
for i = 1:size(locs, 1)
    println(locs[i, 1] * "\t" * locs[i, 2] * "\t" * string(locs[i, 3]) * "\t" * locs[i, 4])
end

# julia fasta_kmer_locations.jl UP000005640_9606.fasta.gz UP000000354_694009.fasta.gz 10 1
