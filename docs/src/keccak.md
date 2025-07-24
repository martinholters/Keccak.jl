# Keccak-based sponge construction

The Keccak algorithm family as defined in
[NIST FIPS-202](https://doi.org/10.6028/NIST.FIPS.202) is made up of a mapping function
(Kᴇᴄᴄᴀᴋ-p) and a padding rule wrapped in a cryptographic sponge.
All Keccak-based hashing algorithms are based on this construction.

The permutation function Kᴇᴄᴄᴀᴋ-p[b,nᵣ] is defined on a bit-string of length b=25⋅2ˡ for
l ∈ {0, …, 6}. This bit-string is treated as a sequence of 25 lanes of length 2ˡ each,
and the corresponding  Keccak.jl implementation [`keccak_p`](@ref) expects an `NTuple{25}`
of unsigned integers holding those lanes. Given the available unsigned integer types,
it follows that the implementation does not support l<3. As all relevant uses of Keccak
set l=6, this restriction does not pose a limitation in practice. As second argument to
[`keccak_p`](@ref), the desired number of rounds nᵣ can be specified (as a `Val`), with the
default of nᵣ=12+2l resulting in the permutation function Kᴇᴄᴄᴀᴋ-f[b].

For example, applying Kᴇᴄᴄᴀᴋ-f[200] to the all-zero state would look as follows:
```jldoctest; setup=:(using Keccak)
julia> keccak_p(ntuple(_ -> zero(UInt8), Val(25)))
(0x3c, 0x28, 0x26, 0x84, 0x1c, 0xb3, 0x5c, 0x17, 0x1e, 0xaa, 0xe9, 0xb8, 0x11, 0x13, 0x4c, 0xea, 0xa3, 0x85, 0x2c, 0x69, 0xd2, 0xc5, 0xab, 0xaf, 0xea)
```

To avoid repeated anonymous functions `s -> keccak_p(s, Val(nr))` fixing the number of
rounds, a wrapper [`KeccakP{nrounds}()`](@ref) is provided with the same effect.

If the input to [`keccak_p`](@ref) is an `NTuple{25}` of `SIMD.Vec`s, it is treated as
multiple states is parallel. That is
```jldoctest; setup=:(using Keccak; import SIMD)
state1 = Tuple(rand(UInt64, 25))
state2 = Tuple(rand(UInt64, 25))
simdstate = SIMD.Vec.(state1, state2) # now an NTuple{25, SIMD.Vec{2, UInt64}}
SIMD.Vec.(keccak_p(state1), keccak_p(state2)) == keccak_p(simdstate)

# output

true
```

Keccak employs the "pad10*1" padding rule, i.e. a 1-bit, an appropriate number of 0-bits,
and another 1-bit are appended to the input to obtain a total length that is a multiple of
the sponge rate.
Typical uses like SHA-3, SHAKE, or cSHAKE insert additional fixed bits immediately before
the padding to achieve domain separation.
This padding with optional domain separation is encapsulated in [`KeccakPad`](@ref), see
there for a discussion how to specify the domain separation bits.
Instances of `KeccakPad` are applied to sponges rather than input data and directly absorb
the padding for performance reasons.
They are used by [`pad`](@ref) internally:
```jldoctest; setup=:(using Keccak)
sha3pad = KeccakPad(0b110) # domain separation 01, first padding 1-bit, but LSB-first, hence reversed
sponge = sha3_256_sponge()
sponge = absorb(sponge, rand(UInt8, 23))
pad(sponge) == sha3pad(sponge)

# output

true
```

Sponges with a combination of [`KeccakP`](@ref) as transformation function and
[`KeccakPad`](@ref) as padding rule are defined by [`KeccakSponge`](@ref).
A `KeccakSponge` is instantiated by calling `KeccakSponge{R,T,nrounds}(pad)`
where `T` is the storage type of the lanes (i.e. unsigned integers of 2ˡ bits) or
`SIMD.Vec`s thereof, `R` is the sponge rate _in number of lanes_, `nrounds` is
the number of rounds to use in Kᴇᴄᴄᴀᴋ-p (defaulting to 12+2l if omitted) and `pad`
is an instance of `KeccakPad` to use for padding.
For example, SHA-3 uses l=6, so `T=UInt64`, and a rate of 1600-2d bits where d is
the digest length. For a digest length of 256, that gives 1088 bits corresponding to
1088/64=17 lanes.
Thus [`sha3_256_sponge()`](@ref) is equal to `KeccakSponge{17,UInt64}(KeccakPad(0b110))`:
```jldoctest; setup=:(using Keccak)
julia> data = "a message to compute the digest of";

julia> digest = squeeze(pad(absorb(KeccakSponge{17,UInt64}(KeccakPad(0b110)), data)), Val(32))[2]
(0xfa, 0x24, 0x49, 0x91, 0xa1, 0x94, 0x44, 0x1e, 0xd6, 0xa8, 0x34, 0x19, 0xd2, 0x00, 0x44, 0x12, 0x15, 0xa9, 0x7c, 0x49, 0xc6, 0xcf, 0x5e, 0x28, 0x7f, 0x75, 0x75, 0x45, 0x73, 0x6d, 0x25, 0x23)

julia> digest == sha3_256(data)
true
```
