"""
    Sponge{R,T<:NTuple,Fxfm}

Cryptographic sponge holding a fixed-length permutation function, state of that length,
and an offset for byte-wise input and output.

A `Sponge` allows absorbing and squeezing arbitrary-sized data chunks. Once per every rate
bytes absorbed or squeezed, the permutation function is invoked. Note that a `Sponge` is
immutable; any operations updating it return an updated `Sponge` instead of mutating the
given one.

# Type parameters
* `T`: type of the `NTuple` holding the sponge contents
* `R`: the rate in units of the state elements (`eltype(T)`)
* `Fxfm`: type of the permutation function

# Fields
* `transform`: the permutation function invoked for every rate bytes absorbed or squeezed.
* `state`: the sponge contents
* `k`: the write (while absorbing) or read (while squeezing) byte offset inside the state,
  between `0` (inclusive) and `R*sizeof(eltype(T))` (exclusive)

# SIMD support
If `T<:SIMD.Vec` (i.e. the state tuple holds `SIMD.Vec`s), absorbing and squeezing treats
the individual data paths separately (and the rate in bytes is given by
`R*sizeof(eltype(eltype(T)))`) in this case.)
"""
struct Sponge{R,T<:NTuple,Fxfm}
    transform::Fxfm
    state::T
    k::Int
end
Sponge{R}(transform::Fxfm, state::T, k) where {R,T,Fxfm} =
    Sponge{R,T,Fxfm}(transform, state, k)

"""
    update(sponge::Sponge, state, k)

Return a new `Sponge` of the same type as `sponge` with the same `transform` function, but
with `state` and `k` replaced with the given values.
"""
update(sponge::S, state::T, k) where {R,T,S<:Sponge{R,T}} = S(sponge.transform, state, k)

"""
    rate(sponge::Sponge)

Returns the rate (in bytes) of the given `sponge`.
"""
rate(::Sponge{R,NTuple{K,T}}) where {R,K,ELT<:Unsigned,T<:Union{ELT,Vec{<:Any,ELT}}} =
    R * sizeof(ELT)

# Return a tuple `out` of length `blocklen`, such that
# `out[destoffset+1:blocklen] == data[1+srcoffset:blocklen-destoffset+srcoffset]`
# for indices within `data` and `0x00` outside
@inline function copy_subchunk(
    data::Union{AbstractVector{UInt8},NTuple{<:Any,UInt8}},
    srcoffset, destoffset,
    ::Val{blocklen}
) where {blocklen}
    ntuple(Val(blocklen)) do i
        @inline
        if !checkindex(Bool, eachindex(data), i-1-destoffset+srcoffset+firstindex(data))
            return 0x00
        else
            return data[i-1-destoffset+srcoffset+firstindex(data)]
        end
    end
end

# return data[srcoffset+1:srcoffset+sizeof(T)*R]` as an `R`-tuple of `T`s
@inline function copy_chunk(data::AbstractVector{UInt8}, srcoffset, ::Type{T}, ::Val{R}) where {T,R}
    srcoffset += firstindex(data)
    return reinterpret(NTuple{R,T}, view(data, (srcoffset:srcoffset+sizeof(T)*R-1)))[1]
end
@inline function copy_chunk(data::NTuple{<:Any,UInt8}, srcoffset, ::Type{T}, ::Val{R}) where {T,R}
    return reinterpret(NTuple{R,T}, ntuple(l -> data[l+srcoffset], Val(sizeof(T)*R)))
end

"""
    absorb(sponge::Sponge, data::Union{AbstractVector{UInt8},NTuple{<:Any,UInt8}})

Absorbs the provided `data` into the `sponge` and returns the updated sponge.

The provided `data` does not have to have a length that is a multiple of the sponge rate.
However, before the first `squeeze`, appropriate padding should be performed.

If the sponge holds `SIMD.Vec`s, the same `data` is used for every data path.
"""
function absorb(
    sponge::Sponge{R,NTuple{K,T}},
    data::Union{AbstractVector{UInt8},NTuple{<:Any,UInt8}}
) where {R,K,ELT<:Unsigned,T<:Union{ELT,Vec{<:Any,ELT}}}
    st = sponge.state
    k = sponge.k
    n = firstindex(data)
    if k != 0
        block = reinterpret(NTuple{R,ELT}, copy_subchunk(data, 0, k, Val(rate(sponge))))
        st = let st=st, block=block
            ntuple(l -> l <= R ? st[l] ⊻ block[l] : st[l], Val(K))
        end
        Δn = min(length(data), rate(sponge)-k)
        n += Δn
        k += Δn
        if k == rate(sponge)
            st = sponge.transform(st)
            k = 0
        end
    end
    while n + rate(sponge) - 1 <= lastindex(data)
        block = copy_chunk(data, n-1, ELT, Val(R))
        st = let st=st, block=block
            ntuple(l -> l <= R ? st[l] ⊻ block[l] : st[l], Val(K))
        end
        Δn = rate(sponge)
        n += Δn
        st = sponge.transform(st)
    end
    if n <= lastindex(data)
        block = reinterpret(NTuple{R,ELT}, copy_subchunk(data, n-1, 0, Val(rate(sponge))))
        st = let st=st, block=block
            ntuple(l -> l <= R ? st[l] ⊻ block[l] : st[l], Val(K))
        end
        Δn = lastindex(data)-n+1
        k += Δn
    end
    return update(sponge, st, k)
end

"""
    absorb(sponge::Sponge, data1, ..., dataN)

Absorbs the provided data into the SIMD `sponge` and returns the updated sponge.

The sponge has to hold `SIMD.Vec{N}`s, matching the number of provided `data` parameters,
and their lengths must be equal, i.e. `length(data1) == ... == length(dataN)` has to hold.

The provided data does not have to have a length that is a multiple of the sponge rate.
However, before the first `squeeze`, appropriate padding should be performed.

"""
function absorb(
    sponge::Sponge{R,NTuple{K,T}},
    data::Vararg{Union{AbstractVector{UInt8},NTuple{<:Any,UInt8}},N}
) where {R,K,N,ELT<:Unsigned,T<:Vec{N,ELT}}
    if !allequal(map(length, data))
        throw(DimensionMismatch("data arguments provided to `absorb` have different length"))
    end
    st = sponge.state
    k = sponge.k
    n = firstindex(data[1])
    if k != 0
        blocks = let k=k
            reinterpret(
                NTuple{N,NTuple{R,ELT}},
                map(@inline(d -> copy_subchunk(d, 0, k, Val(rate(sponge)))), data),
            )
        end
        st = let st=st, blocks=blocks
            ntuple(l -> l <= R ? st[l] ⊻ Vec{N,ELT}(map(b->b[l], blocks)) : st[l], Val(K))
        end
        Δn = min(length(data[1]), rate(sponge)-k)
        n += Δn
        k += Δn
        if k == rate(sponge)
            st = sponge.transform(st)
            k = 0
        end
    end
    while n + rate(sponge) - 1 <= lastindex(data[1])
        blocks = let n=n; map(@inline(d -> copy_chunk(d, n-1, ELT, Val(R))), data); end
        st = let st=st, blocks=blocks
            ntuple(l -> l <= R ? st[l] ⊻ Vec{N,ELT}(map(b->b[l], blocks)) : st[l], Val(K))
        end
        Δn = rate(sponge)
        n += Δn
        st = sponge.transform(st)
    end
    if n <= lastindex(data[1])
        blocks = let n=n
            reinterpret(
                NTuple{N,NTuple{R,ELT}},
                map(@inline(d -> copy_subchunk(d, n-1, 0, Val(rate(sponge)))), data)
            )
        end
        st = let st=st, blocks=blocks
            ntuple(l -> l <= R ? st[l] ⊻ Vec{N,ELT}(map(b->b[l], blocks)) : st[l], Val(K))
        end
        Δn = lastindex(data[1])-n+1
        k += Δn
    end
    return update(sponge, st, k)
end

# prevent ambiguity and solve be reinterpreting `Vec{1,T}` as `T`
function absorb(
    sponge::Sponge{R,NTuple{K,T}},
    data::Union{AbstractVector{UInt8},NTuple{<:Any,UInt8}}
) where {R,K,ELT<:Unsigned,T<:Vec{1,ELT}}
    sponge′ = @inline absorb(
        Sponge{R}(
            sponge.transform,
            reinterpret(NTuple{K,ELT}, sponge.state),
            sponge.k,
        ),
        data,
    )
    return Sponge{R}(
        sponge′.transform,
        reinterpret(NTuple{K,T}, sponge′.state),
        sponge′.k,
    )
end

"""
    sponge′, data = squeeze(sponge, len)

Squeezes `len` bytes from the `sponge` and returns the updated sponge and the obtained data.

If `len` is a `Val`, the `data` is returned as a tuple, otherwise as a `Memory{UInt8}`.

If the `sponge` holds `SIMD.Vec{N}`s, the returned `data` is an `N`-tuple of the data
squeezed form the respective data paths.
"""
function squeeze end

function squeeze(sponge::Sponge{R,NTuple{K,T}}, len) where {R,K,T<:Unsigned}
    J = sizeof(T)
    data = Memory{UInt8}(undef, len)
    k = sponge.k
    st = sponge.state
    n = 1
    while n <= len
        curlen = min(len-n+1, rate(sponge)-k)
        copyto!(data, n, reinterpret(NTuple{J*K,UInt8}, st), k+1, curlen)
        k += curlen
        n += curlen
        if k == rate(sponge)
            st = sponge.transform(st)
            k = 0
        end
    end
    return update(sponge, st, k), data
end

function squeeze(
    sponge::Sponge{R,NTuple{K,T}},
    len,
) where {R,K,N,ELT<:Unsigned,T<:Vec{N,ELT}}
    data = ntuple(_ -> Memory{UInt8}(undef, len), Val(N))
    k = sponge.k
    st = sponge.state
    n = 1
    while n <= len
        curlen = min(len-n+1, rate(sponge)-k)
        for i in 1:N
            sti = let i=i, st=st
                ntuple(j -> st[j][i], Val(R))
            end
            copyto!(data[i], n, reinterpret(NTuple{rate(sponge),UInt8}, sti), k+1, curlen)
        end
        k += curlen
        n += curlen
        if k == rate(sponge)
            st = sponge.transform(st)
            k = 0
        end
    end
    return update(sponge, st, k), data
end

function squeeze(sponge::Sponge{R,NTuple{K,T}}, ::Val{len}) where {R,K,T<:Unsigned,len}
    J = sizeof(T)
    k = sponge.k
    st = sponge.state
    if len+k <= rate(sponge)
        data = let k=k, st=reinterpret(NTuple{J*K,UInt8}, st)
            ntuple(i -> st[i+k], Val(len))
        end
        k += len
        if k == rate(sponge)
            st = sponge.transform(st)
            k = 0
        end
        return update(sponge, st, k), data
    elseif len+k <= 2*rate(sponge)
        st′ = sponge.transform(st)
        data = let k=k, st1=reinterpret(NTuple{J*K,UInt8}, st), r=rate(sponge)
            st2=reinterpret(NTuple{J*K,UInt8}, st′)
            ntuple(i -> i+k <= r ? st1[i+k] : st2[i+k-r], Val(len))
        end
        k += len - rate(sponge)
        if k == rate(sponge)
            st′ = sponge.transform(st′)
            k = 0
        end
        return update(sponge, st′, k), data
    else
        sponge, data′ = squeeze(sponge, len)
        return sponge, ntuple(i -> data′[i], Val(len))
    end
end

function squeeze(
    sponge::Sponge{R,NTuple{K,T}},
    ::Val{len}
) where {R,K,N,ELT<:Unsigned,T<:Vec{N,ELT},len}
    k = sponge.k
    st = sponge.state
    if len+k <= rate(sponge)
        states = let st=st
            ntuple(@inline(n -> ntuple(i -> st[i][n], Val(R))), Val(N))
        end
        data = let k=k, states=states, r=rate(sponge)
            ntuple(Val(N)) do n
                let st=reinterpret(NTuple{r,UInt8}, states[n])
                    ntuple(i -> st[i+k], Val(len))
                end
            end
        end
        k += len
        if k == rate(sponge)
            st = sponge.transform(st)
            k = 0
        end
        return update(sponge, st, k), data
    elseif len+k <= 2*rate(sponge)
        states = let st=st
            ntuple(@inline(n -> ntuple(i -> st[i][n], Val(R))), Val(N))
        end
        st′ = sponge.transform(st)
        states′ = let st=st′
            ntuple(@inline(n -> ntuple(i -> st[i][n], Val(R))), Val(N))
        end
        data = let k=k, states=states, r=rate(sponge)
            @inline ntuple(Val(N)) do n
                @inline
                let st1=reinterpret(NTuple{r,UInt8}, states[n]),
                    st2=reinterpret(NTuple{r,UInt8}, states′[n])
                    ntuple(i -> i+k <= r ? st1[i+k] : st2[i+k-r], Val(len))
                end
            end
        end
        k += len - rate(sponge)
        if k == rate(sponge)
            st′ = sponge.transform(st′)
            k = 0
        end
        return update(sponge, st′, k), data
    else
        sponge, data′ = squeeze(sponge, len)
        return sponge, ntuple(n -> ntuple(i -> data′[n][i], Val(len)), Val(N))
    end
end
