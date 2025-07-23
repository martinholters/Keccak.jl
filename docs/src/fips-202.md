# SHA-3 (and SHAKE)

## SHA-3

The most commonly used functionality contained in Keccak.jl probably is computation of
digests according to the Secure Hash Algorithm-3 (SHA-3) as specified in
[NIST FIPS-202](https://doi.org/10.6028/NIST.FIPS.202).
SHA-3 is a family of hashing functions differing in the length of the generated digits,
offering 224 bits, 256 bits, 384 bits, or 512 bits of output length. These map to the
functions [`sha3_224`](@ref), [`sha3_256`](@ref), [`sha3_384`](@ref), and
[`sha3_512`](@ref), respectively. In the following, examples will be limited to the
256 bits version; the other lengths work similarly.

### Basic usage

We start with a very basic example:
```jldoctest
julia> using Keccak: sha3_256

julia> sha3_256("a message to compute the digest of")
(0xfa, 0x24, 0x49, 0x91, 0xa1, 0x94, 0x44, 0x1e, 0xd6, 0xa8, 0x34, 0x19, 0xd2, 0x00, 0x44, 0x12, 0x15, 0xa9, 0x7c, 0x49, 0xc6, 0xcf, 0x5e, 0x28, 0x7f, 0x75, 0x75, 0x45, 0x73, 0x6d, 0x25, 0x23)
```
As can be seen, the output is a tuple of `UInt8`s, namely 32 of them for a total of the
required 32⋅8=256 bits. The interface is very similar to that of the `SHA` standard library:
```jldoctest
julia> using SHA: sha3_256

julia> sha3_256("a message to compute the digest of")
32-element Vector{UInt8}:
 0xfa
 0x24
 0x49
 0x91
 0xa1
 0x94
 0x44
 0x1e
 0xd6
 0xa8
    ⋮
 0x28
 0x7f
 0x75
 0x75
 0x45
 0x73
 0x6d
 0x25
 0x23
```
The `SHA` standard library produces the same output, of course, but stores it in a `Vector`.
!!! note
    But `SHA` and `Keccak` export functions of the same name, so care must be taken when
    importing both of them to prevent name clashes. It is recommended to import only those
    symbols actually used, e.g. `using Keccak: sha3_256` instead of `using Keccak`.

Apart from the return type, there are two noteworthy differences between the `SHA` and the
`Keccak` implementations:
* `Keccak` only allows `AbstractVector{UInt8}` und `Tuple{Vararg{UInt8}}` input in addition
  to the `String` input demonstrated above, while `SHA` also supports `IO` to hash all data
  coming from an `IO` object (e.g. a file).
* Performance:
  ```jldoctest; filter=r"\d*(\.\d* .?s)?"
  using Chairmarks: @b
  import Keccak, SHA
  @b rand(UInt8, 1_000_000) Keccak.sha3_256, SHA.sha3_256

  # output

  (1.803 ms, 6.684 ms (9 allocs: 720 bytes))
  ```
  The `Keccak` implementation is faster and avoids allocations.
  (Comparison was done using Julia v1.12.0-rc1).

### Chunked input (and output)

When e.g. needing to hash a large file, it may be inappropriate to read it into memory as a
whole. Rather, one would like to process the data in reasonably-sized chunks. This is
possible using the [sponge-based interface](@ref "Sponge operations") of `Keccak`.
(The term "sponge" stems from the
algorithm family underlying SHA-3 -- Keccak -- using a so-called cryptographic sponge
construction.)

Using the sponge-based interface entails the following steps:
1. Obtain a suitable sponge with e.g. [`sha3_256_sponge`](@ref).
2. Process the input data with zero or more calls to [`absorb`](@ref).
   (Zero calls correspond to an empty input.)
3. Call [`pad`](@ref) _exactly_ once.
4. Produce the output by one or more calls to [`squeeze`](@ref).
   (Zero calls are technically permitted, too, but pointless.)

!!! warning
    At present, adhering to above sequence is not enforced, but the behavior for any other
    sequence of [`absorb`](@ref), [`pad`](@ref) and [`squeeze`](@ref) invocations is
    unspecified.

An important feature of the sponges used by `Keccak` is that they are immutable. Therefore,
none of the operations above mutate the given sponge in-place; rather, they return an
updated sponge.

The following example shows both chunked input (multiple calls to `absorb`) as well as
chunked output (multiple calls to `squeeze`), although the latter is certainly less useful
is this context:
```jldoctest; setup=:(using Keccak)
julia> sponge = sha3_256_sponge();

julia> sponge = absorb(sponge, 0x00:0x04); # absorb first data chunk

julia> sponge = absorb(sponge, 0x05:0x09); # absorb another data chunk

julia> sponge = pad(sponge); # absorb appropriate padding

julia> sponge, out1 = squeeze(sponge, Val(16)); # first part of output

julia> sponge, out2 = squeeze(sponge, Val(16)); # second part of output

julia> (out1..., out2...) == sha3_256(0x00:0x09) # same result
true
```
Note that the desired output length of [`squeeze`](@ref) has to be given in bytes.
If the number is passed directly instead of wrapped in a `Val`, the output will be a
`Vector{UInt8}` instead of a tuple.
!!! note
    Squeezing fewer or more bytes from the sponge than the standard demands is perfectly
    valid technically, but obviously not standard-compliant. And as the output length
    matches the security strength, squeezing more bytes will not produce a more secure
    digest. An application scenario where squeezing more bytes makes sense is to use the
    hashing function as a pseudo-random function generator. If this is your aim, consider
    [SHAKE](@ref), [cSHAKE](@ref), or the extensible output variants of [KMAC](@ref),
    [TupleHash](@ref), or [ParallelHash](@ref).

## SHAKE

The [NIST FIPS-202 standard](https://doi.org/10.6028/NIST.FIPS.202) specifies an
extensible-output variant of SHA-3, called SHAKE, for security strengths 128 bits and
256 bits. These are available in `Keccak` as [`shake_128`](@ref) and
[`shake_256`](@ref), respectively.

Similar to the [SHA-3](@ref) functions, one can directly compute a SHAKE-digest from
input data, but has to pass the desired output length:
```jldoctest; setup=:(using Keccak)
julia> msg = "a message to compute the digest of";

julia> shake_128(msg, 10) # returns a vector
10-element Vector{UInt8}:
 0xd4
 0xa7
 0x25
 0x77
 0xde
 0x29
 0x05
 0x20
 0x9b
 0x35

julia> shake_128(msg, Val(10)) # returns a tuple
(0xd4, 0xa7, 0x25, 0x77, 0xde, 0x29, 0x05, 0x20, 0x9b, 0x35)

julia> shake_128(msg, Val(15)) # longer output, first 10 bytes equal
(0xd4, 0xa7, 0x25, 0x77, 0xde, 0x29, 0x05, 0x20, 0x9b, 0x35, 0x64, 0x68, 0x9a, 0x96, 0xad)
```

Chunked input and output is possible in the same way as for SHA-3 (see
[above](@ref "Chunked input (and output)")), replacing `sha3_256_sponge` with
[`shake_128_sponge`](@ref) (or [`shake_256_sponge`](@ref)). For convenience, one can
also obtain a sponge ready for squeezing by calling [`shake_128(data)`](@ref) without
specifying the output length:
```jldoctest; setup=:(using Keccak)
julia> msg = "a message to compute the digest of";

julia> sponge = shake_128(msg);

julia> sponge, out1 = squeeze(sponge, Val(10));

julia> out1 # as above
(0xd4, 0xa7, 0x25, 0x77, 0xde, 0x29, 0x05, 0x20, 0x9b, 0x35)

julia> sponge, out2 = squeeze(sponge, Val(5));

julia> out2 # last five bytes of the length-15 example above
(0x64, 0x68, 0x9a, 0x96, 0xad)
```
