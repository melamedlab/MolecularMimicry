module SkewSuffixArray

export suffixArray!

function radixPass!(a::Array{Int64}, b::Array{Int64}, r::Array{Int64},
                    rShift::Int64, n::Int64, K::Int64)::Nothing
    Kp1 = K+1
    ## count occurrences
    c = zeros(Int64, Kp1)
    @inbounds for i = 1:n
        c[r[rShift + a[i] + 1] + 1] += 1
    end
    ## exclusive prefix sums
    sum = 0
    for i = 1:Kp1
        t = c[i]
        c[i] = sum
        sum += t
    end
    ## sort
    @inbounds for i = 1:n
        b[c[r[rShift + a[i] + 1] + 1] + 1] = a[i]
        c[r[rShift + a[i] + 1] + 1] += 1        
    end
    return nothing
end

function suffixArray!(s::Array{Int64}, SA::Array{Int64}, n::Int64, K::Int64)::Nothing
    n0=(n+2)÷3; n1=(n+1)÷3; n2=n÷3; n02=n0+n2
    s12 = zeros(Int64, n02+3)
    SA12 = zeros(Int64, n02+3)
    s0 = Array{Int64}(undef, n0)
    SA0 = Array{Int64}(undef, n0)
    ##
    ## generate positions of mod 1 and mod 2 suffixes
    ## the "+(n0-n1)" adds a dummy mod 1 suffix if n%3 == 1
    let j = 1
        @inbounds for i = 0:(n+(n0-n1)-1)
            if i%3 != 0
                s12[j] = i
                j += 1            
            end
        end
    end
    ##
    ## lsb radix sort the mod 1 and mod 2 triples
    radixPass!(s12, SA12, s, 2, n02, K)
    radixPass!(SA12, s12, s, 1, n02, K)
    radixPass!(s12, SA12, s, 0, n02, K)
    ##
    ## find lexicographic names of triples
    let name = 0, c0 = -1, c1 = -1, c2 = -1
        @inbounds for i = 1:n02
            if s[SA12[i]+1] != c0 || s[SA12[i]+2] != c1 || s[SA12[i]+3] != c2
                name += 1
                c0 = s[SA12[i]+1]
                c1 = s[SA12[i]+2]
                c2 = s[SA12[i]+3]
            end
            if (SA12[i] % 3) == 1
                ## left half
                s12[SA12[i]÷3 + 1] = name
            else
                ## right half
                s12[SA12[i]÷3 + n0+1] = name
            end
        end
        ##
        ## recurse if names are not yet unique
        if name < n02
            suffixArray!(s12, SA12, n02, name)
            ## store unique names in s12 using the suffix array
            @inbounds for i = 1:n02
                s12[SA12[i]+1] = i
            end
        else
            ## generate the suffix array of s12 directly
            @inbounds for i = 1:n02
                SA12[s12[i]] = i-1
            end
        end
    end
    ##
    ## stably sort the mod 0 suffixes from SA12 by their first character
    let j = 1
        @inbounds for i = 1:n02
            if SA12[i] < n0
                s0[j] = 3 * SA12[i]
                j += 1
            end
        end
    end
    radixPass!(s0, SA0, s, 0, n0, K)
    ##
    ## merge sorted SA0 suffixes and sorted SA12 suffixes
    let p=0; t=n0-n1; k=0
        while k < n
            ## pos of current offset 12 suffix
            i = (SA12[t+1] < n0 ? SA12[t+1] * 3 + 1 : (SA12[t+1] - n0) * 3 + 2)
            ## pos of current offset 0 suffix
            j = SA0[p+1]
            ##
            if (SA12[t+1] < n0 ?
                (s[i+1], s12[SA12[t+1]+n0+1]) < (s[j+1], s12[j÷3+1]) :
                (s[i+1], s[i+2], s12[SA12[t+1]-n0+2]) < (s[j+1], s[j+2], s12[j÷3+n0+1]))
                ## suffix from SA12 is smaller
                SA[k+1] = i
                t += 1
                if t == n02
                    ## done --- only SA0 suffixes left
                    k += 1
                    while p < n0
                        SA[k+1] = SA0[p+1]
                        p += 1
                        k += 1
                    end
                end
            else
                SA[k+1] = j
                p += 1                
                if p == n0
                    ## done --- only SA12 suffixes left
                    k += 1
                    while t < n02
                        SA[k+1] = (SA12[t+1] < n0 ?
                                   SA12[t+1] * 3 + 1 :
                                   (SA12[t+1] - n0) * 3 + 2)
                        t += 1
                        k += 1
                    end
                end
            end
            k += 1
        end
    end
    return nothing
end

end # module
