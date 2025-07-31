"""
    KeccakP{nrounds}

Convenience wrapper for `keccak_p` with `nrounds` fixed.

This allows replacing constructions like `state -> keccak_p(state, Val(nrounds))`
with `KeccakP{nrounds}()`:
```jldoctest; setup=:(using Keccak; state=Tuple(zeros(UInt64,25))), filter=r"^#\\d*"
julia> const f1 = KeccakP{23}()
KeccakP{23}()

julia> const f2 = state -> keccak_p(state, Val(23))
#2 (generic function with 1 method)

julia> f1(state) == f2(state) == keccak_p(state, Val(23)) # for suitable state
true
```
"""
struct KeccakP{nrounds} end
(::KeccakP{nrounds})(state::NTuple{25}) where {nrounds} = keccak_p(state, Val(nrounds))

"""
    KeccakPad

A padding function performing the Kᴇᴄᴄᴀᴋ pad10*1 padding with optional suffix bits inserted
first for domain separation. These suffix bits are specified by providing a byte that
includes the first 1 bit of the padding.

E.g. for SHA3, the domain separation suffix is `01`, so the padding appends
`011000...0001`. Noting the LSB-first convention, the first byte of the padding thus has to
be `0b00000110 == 0x06`, and the corresponding padding function would be instantiated with
`KeccakPad(0x06)`.

!!! note
    Only up to six domain separation suffix bits are supported, so including the first
    padding 1 bit, the first byte has be between 0b00000001 (no suffix bits) and 0b01111111
    (suffix 111111).
"""
struct KeccakPad
    firstbyte::UInt8
    function KeccakPad(firstbyte=0x01)
        if !(0x01 <= firstbyte <= 0x7f)
            throw(ArgumentError("invalid first byte of padding"))
        end
        return new(firstbyte)
    end
end
function (pad::KeccakPad)(sponge::Sponge{R, NTuple{K, T}} where {T}) where {R, K}
    ELT = lanetype(sponge)
    J = sizeof(ELT)
    i, j = divrem(sponge.k, J)
    st = let st=sponge.state
        ntuple(Val(K)) do l
            s = st[l]
            if l == i+1
                s ⊻= ELT(pad.firstbyte) << (8j)
            end
            if l == R
                s ⊻= ELT(0x80) << (8*(J-1))
            end
            return s
        end
    end
    st = sponge.transform(st)
    return Keccak.update(sponge, st, 0)
end

"""
    KeccakSponge{R,T,nrounds}

A `Sponge` specialization for Kᴇᴄᴄᴀᴋ:
- The permutation function is `KeccakP{nrounds}`.
- The sponge contents are a 25-tuple of `T`, which must be an unsigned integer or a
  `SIMD.Vec` thereof.
- The rate in bits is given by `8*sizeof(T)*R` (or `8*sizeof(eltype(T))*R)` if
  `T<:SIMD.Vec`).
"""
const KeccakSponge{R,T<:Union{Unsigned,<:Vec{<:Any,<:Unsigned}},nrounds} =
    Sponge{R,NTuple{25,T},KeccakP{nrounds},KeccakPad}

"""
    KeccakSponge{R,T,nrounds}(pad)

Return a zero-initialized `KeccakSponge{R,T,nrounds}` with the provided padding
function `pad`.
"""
function KeccakSponge{R,T,nrounds}(pad::KeccakPad) where {R,T,nrounds}
    Sponge{R}(
        KeccakP{nrounds}(),
        pad,
        ntuple(_ -> zero(T), Val(25)),
        0,
    )
end

"""
    KeccakSponge{R,T}(pad)

Return a zero-initialized `KeccakSponge{R,T,nrounds}` with the default `nrounds=12+2ℓ(T)`
and the provided padding function `pad`.
"""
KeccakSponge{R,T}(pad::KeccakPad) where {R,T} = KeccakSponge{R,T,12+2ℓ(T)}(pad)
