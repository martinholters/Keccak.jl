"""
    ℓ(::Type{T<:Unsigned})

Returns the base-2 logarithm of the bitwidth of values of type `T`.

# Examples
```jldoctest; setup=:(using Keccak: ℓ)
julia> ℓ(UInt32)
5

julia> 2^5
32
```
"""
function ℓ end

ℓ(::Type{UInt8}) = 3
ℓ(::Type{UInt16}) = 4
ℓ(::Type{UInt32}) = 5
ℓ(::Type{UInt64}) = 6

"""
    ℓ(::Type{SIMD.Vec{N,T<:Unsigned}})

Returns the base-2 logarithm of the bitwidth of values of type `T`, i.e. the element type of
the `Vec`.

# Examples
```jldoctest; setup=:(using Keccak: ℓ; import SIMD)
julia> ℓ(SIMD.Vec{4,UInt32})
5

julia> 2^5
32
```
"""
ℓ(::Type{Vec{N, T}}) where {N, T<:Unsigned} = ℓ(T)

"""
    θ(state)

Perform the step mapping θ of the Keccak algorithm and return the updated state.

The state is a 25-tuple of unsigned integers, related to state array A in FIPS-202 by

    A[x,y,z] == state[5y+x+1] >> z & 1

for x = 0,...,4, y = 0,...,4 and z = 0,...,w-1, where w is the bitwidth of the integers in
`state`.

Alternatively, the state can contain `SIMD.Vec`s of unsigned integers to perform a parallel
transformation of multiple state instances.
"""
@inline function θ(state::NTuple{25})
    C = ntuple(Val(5)) do i
        return state[i] ⊻ state[i + 5] ⊻ state[i + 10] ⊻ state[i + 15] ⊻ state[i + 20]
    end
    D = ntuple(@inline(i -> C[mod1(i - 1, 5)] ⊻ bitrotate(C[mod1(i + 1, 5)], 1)), Val(5))
    return ntuple(k -> state[k] ⊻ D[mod1(k, 5)], Val(25))
end

"""
    ρ(state)

Perform the step mapping ρ (bit rotation in each lane) of the Keccak algorithm and return
the updated state.

The state is a 25-tuple of unsigned integers, related to state array A in FIPS-202 by

    A[x,y,z] == state[5y+x+1] >> z & 1

for x = 0,...,4, y = 0,...,4 and z = 0,...,w-1, where w is the bitwidth of the integers in
`state`.

Alternatively, the state can contain `SIMD.Vec`s of unsigned integers to perform a parallel
transformation of multiple state instances.
"""
@inline function ρ(state::NTuple{25})
    SHA3_ROTC = ( # table 2 in FIPS-202, reordered to 0...4 indexing
          0,   1, 190,  28,  91,
         36, 300,   6,  55, 276,
          3,  10, 171, 153, 231,
        105,  45,  15,  21, 136,
        210,  66, 253, 120,  78,
    )
    return ntuple(k -> bitrotate(state[k], SHA3_ROTC[k]), Val(25))
end

"""
    π(state)

Perform the step mapping π (permutation of the lanes) of the Keccak algorithm and return
the updated state.

The state is a 25-tuple of unsigned integers, related to state array A in FIPS-202 by

    A[x,y,z] == state[5y+x+1] >> z & 1

for x = 0,...,4, y = 0,...,4 and z = 0,...,w-1, where w is the bitwidth of the integers in
`state`.

Alternatively, the state can contain `SIMD.Vec`s of unsigned integers to perform a parallel
transformation of multiple state instances.
"""
@inline π(state::NTuple{25}) =
    ntuple(Val(25)) do k
        @inline
        y, x = divrem(k-1, 5)
        state[mod(x+3y, 5) + 5*x + 1]
    end

"""
    χ(state)

Perform the nonlinear step mapping χ of the Keccak algorithm and return the updated state.

The state is a 25-tuple of unsigned integers, related to state array A in FIPS-202 by

    A[x,y,z] == state[5y+x+1] >> z & 1

for x = 0,...,4, y = 0,...,4 and z = 0,...,w-1, where w is the bitwidth of the integers in
`state`.

Alternatively, the state can contain `SIMD.Vec`s of unsigned integers to perform a parallel
transformation of multiple state instances.
"""
@inline χ(state::NTuple{25}) =
    ntuple(Val(25)) do k
        @inline
        y, x = divrem(k-1, 5)
        state[k] ⊻ (~state[rem(x+1, 5) + 5y + 1] & state[rem(x+2, 5) + 5y + 1])
    end

"""
    rc(t)

Evaluate `t` rounds of a linear feedback shift register and return the leading bit after
the final round as specified in Algorithm 5 of FIPS-202.

This function is used for computing the round constants by `round_consts`.
"""
Base.@assume_effects :terminates_locally function rc(t)
    if t % 255 == 0
        return true
    end
    R = 0x80
    for _ in 1:mod(t, 255)
        R8 = (R & 0x01)
        R >>= 1
        R ⊻= (R8 << 7) | (R8 << 3) | (R8 << 2) | (R8 << 1)
    end
    return R & 0x80 != 0
end

"""
    round_consts(::Type{T}, ::Val{nrounds}=Val(12+2ℓ(T)))

Compute the Keccak round constants for lane length w corresponding to the bitwidth of `T`
and the given number fo rounds `nrounds`.

The return value is a tuple of `T`s containing the round constants for round indices
`12+2ℓ(T)-nrounds` to `12+2ℓ(T)-1` in that order. The default for `nrounds` results in the first round index being 0.

!!! note
    It is valid to have `nrounds>12+2ℓ(T)`, so that round index starts negative.

# Examples
```jldoctest; setup=:(using Keccak: round_consts)
julia> round_consts(UInt16)
(0x0001, 0x8082, 0x808a, 0x8000, 0x808b, 0x0001, 0x8081, 0x8009, 0x008a, 0x0088, 0x8009, 0x000a, 0x808b, 0x008b, 0x8089, 0x8003, 0x8002, 0x0080, 0x800a, 0x000a)

julia> round_consts(UInt8)
(0x01, 0x82, 0x8a, 0x00, 0x8b, 0x01, 0x81, 0x09, 0x8a, 0x88, 0x09, 0x0a, 0x8b, 0x8b, 0x89, 0x03, 0x02, 0x80)

julia> round_consts(UInt8, Val(10))
(0x8a, 0x88, 0x09, 0x0a, 0x8b, 0x8b, 0x89, 0x03, 0x02, 0x80)

julia> round_consts(UInt8, Val(20))
(0x02, 0x8a, 0x01, 0x82, 0x8a, 0x00, 0x8b, 0x01, 0x81, 0x09, 0x8a, 0x88, 0x09, 0x0a, 0x8b, 0x8b, 0x89, 0x03, 0x02, 0x80)
```
"""
@inline function round_consts(::Type{T}, ::Val{nrounds}=Val(12+2ℓ(T))) where {T,nrounds}
    lastround = 12+2ℓ(T)-1
    return ntuple(Val(nrounds)) do round
        @inline
        ir = round - nrounds + lastround
        reduce(|, ntuple(@inline(j -> convert(T, rc((j-1)+7*ir))<<(2^(j-1)-1)), Val(ℓ(T)+1)))
    end
end
@inline round_consts(::Type{Vec{N,T}}, nr::Val) where {N,T} = round_consts(T, nr)

"""
    ι(state, RC)

Perform the step mapping ι (round-dependent modification of lane[0,0]) of the Keccak
algorithm and return the updated state.

Instead of the round index as in FIPS-202, the corresponding round constant `RC` is expected
as the second argument.

The state is a 25-tuple of unsigned integers, related to state array A in FIPS-202 by

    A[x,y,z] == state[5y+x+1] >> z & 1

for x = 0,...,4, y = 0,...,4 and z = 0,...,w-1, where w is the bitwidth of the integers in
`state`.

Alternatively, the state can contain `SIMD.Vec`s of unsigned integers to perform a parallel
transformation of multiple state instances.
"""
@inline ι(state::NTuple{25,Union{T,Vec{<:Any,T}}}, RC::T) where {T} =
    (state[1] ⊻ RC, state[2:end]...)

"""
    keccak_p(state, ::Val{nrounds})

The Kᴇᴄᴄᴀᴋ-p permutation function as defined in FIPS-202, with the number of rounds given
by `nrounds`. The input state as well as the returned updated state are 25-tuples of
unsigned integers of type `T`. Their LSB-first serialization corresponds to the bitstrings S
and S′ in FIPS-202; its length b is determined by the type `T`. E.g. for `T==UInt64`, one
obtains the commonly used b = 25 ⋅ 64 = 1600.

For example, for b=800 (hence `T==UInt32`), S=1000...000 corresponds to
`state == (0x00000001, (0x00000000 for _ in 2:25)...)`.

Alternatively, the state can contain `SIMD.Vec`s of unsigned integers to perform a parallel
transformation of multiple state instances. E.g, for `state::NTuple{25,SIMD.Vec{4,UInt64}}`,
it holds that
```jldoctest; setup = :(import SIMD; using Keccak; state = ((SIMD.Vec(rand(UInt64,4)...) for _ in 1:25)...,))
julia> SIMD.Vec.((keccak_p(map(s -> s[n], state)) for n in 1:4)...) == keccak_p(state)
true
```
"""
function keccak_p(state::NTuple{25,T}, ::Val{nrounds}=Val(12+2ℓ(T))) where {T,nrounds}
    RCs = round_consts(T, Val(nrounds))
    for ir in (12+2ℓ(T)-nrounds):(12+2ℓ(T)-1)
        state = θ(state)
        state = ρ(state)
        state = π(state)
        state = χ(state)
        RC = RCs[ir+nrounds-12-2ℓ(T)+1]
        state = ι(state, RC)
    end
    return state
end

"""
    KeccakP{nrounds}

Convenience wrapper for `keccak_p` with `nrounds` fixed.

This allows replacing constructions like `state -> keccak_p(state, Val(nrounds))`
with `KeccakP{nrounds}()`:
```jldoctest; setup=:(using Keccak; state=Tuple(zeros(UInt64,25)))
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
    KeccakSponge{R,T,nrounds}

A `Sponge` specialization for Kᴇᴄᴄᴀᴋ:
- The permutation function is `KeccakP{nrounds}`.
- The sponge contents are a 25-tuple of `T`, which must be an unsigned integer or a
  `SIMD.Vec` thereof.
- The rate in bits is given by `8*sizeof(T)*R` (or `8*sizeof(eltype(T))*R)` if
  `T<:SIMD.Vec`).
"""
const KeccakSponge{R,T<:Union{Unsigned,<:Vec{<:Any,<:Unsigned}},nrounds} =
    Sponge{R,NTuple{25,T},KeccakP{nrounds}}

"""
    KeccakSponge{R,T,nrounds}(pad)

Return a zero-initialized `KeccakSponge{R,T,nrounds}` with the provided padding
function `pad`.
"""
function KeccakSponge{R,T,nrounds}(pad::F) where {R,T,nrounds,F}
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
KeccakSponge{R,T}(pad::F) where {R,T,F} = KeccakSponge{R,T,12+2ℓ(T)}(pad)
