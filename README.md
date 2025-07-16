# Keccak.jl - Keccak-based hashing (SHA3, SHAKE, ...)

[![CI](https://github.com/martinholters/Keccak.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/martinholters/Keccak.jl/actions/workflows/ci.yml)

See the
[SHA-3 Standard: Permutation-Based Hash and Extendable-Output Functions][NIST-FIPS-202]
for details on the algorithms implemented in this package.

## Basic usage

### SHA3

SHA3 hashes can be computed with `sha3_224`, `sha3_256`, `sha3_384`, and `sha3_512`, which
behave similar to their cousins in `SHA.jl`, with two important differences:
* The input is restricted to being an `AbstractVector{UInt8}`, a `Tuple{Vararg{UInt8}}` or
  a `String`.
* The output is an `NTuple{N,UInt8}` (with `N` depending on hash length, e.g.
  `N==256÷8==32` for SHA3-256).

For example:
```julia-repl
julia> using Keccak

julia> import SHA

julia> msg = rand(UInt8, 1_000_000);

julia> sha3_256(msg)
(0xe0, 0x66, 0xc7, 0xa3, 0xc8, 0xd9, 0xe0, 0xbb, 0xac, 0x21, 0x08, 0xd1, 0x2b, 0x8b, 0xb8, 0x77, 0xf4, 0x15, 0xd2, 0x02, 0x13, 0x1f, 0x8a, 0x29, 0xf6, 0x6e, 0xaa, 0xce, 0x40, 0x86, 0xe6, 0x24)

julia> sha3_256(msg) == Tuple(SHA.sha3_256(msg))
true
```

A motivation for using `Keccak.jl` instead of `SHA.jl` is that it does not allocate and is
more efficient:
```julia-repl
julia> using Chairmarks

julia> @b msg sha3_256, SHA.sha3_256
(1.774 ms, 6.411 ms (9 allocs: 704 bytes))

julia> VERSION
v"1.11.5"
```

### SHAKE

SHAKE hashes can be computed similarly with `shake_128` and `shake_256`, but the output
length (in bytes) has to be given explicitly, like so:
```julia-repl
julia> shake_128(msg, 10) # returns a vector
10-element Vector{UInt8}:
 0x5c
 0x6b
 0xe2
 0x0c
 0xe4
 0xe1
 0x7c
 0x14
 0xbb
 0xf8

julia> shake_128(msg, Val(10)) # returns a tuple
(0x5c, 0x6b, 0xe2, 0x0c, 0xe4, 0xe1, 0x7c, 0x14, 0xbb, 0xf8)

julia> shake_128(msg, Val(15)) # longer output, starting equal
(0x5c, 0x6b, 0xe2, 0x0c, 0xe4, 0xe1, 0x7c, 0x14, 0xbb, 0xf8, 0x40, 0x8d, 0xf9, 0x30, 0x89)
```

Note that SHAKE-128 and SHAKE-256 differ in cryptographic strength (and the produced
output, of course), but do not prescribe a certain output length.

## Advanced Usage

### Piece-wise input/output

The Keccak hashes are based on a "sponge" construction, where data is first "absorbed" and
the hash is then "squeezed" from the sponge. However, the input data has to be padded
appropriately. `Keccak.jl` makes these operations available as `absorb`, `squeeze`, and
`pad`, where the latter is applied to the sponge rather than the input data for efficiency.
The `Keccak.jl` sponges are immutable, so instead of in-place modification, an updated
sponge is returned by these functions.

For example:
```julia-repl
julia> sponge = sha3_256_sponge();

julia> sponge = absorb(sponge, 0x00:0x04); # absorb first data chunk

julia> sponge = absorb(sponge, 0x05:0x09); # absorb another data chunk

julia> sponge = pad(sponge); # absorb appropriate padding

julia> sponge, out1 = squeeze(sponge, Val(16)); # first part of output

julia> sponge, out2 = squeeze(sponge, Val(16)); # second part of output

julia> (out1..., out2...) == sha3_256(0x00:0x09) # same result
true
```

Note that squeezing more (or less) data from a SHA-3 sponge than demanded by the
[standard][NIST-FIPS-202] is perfectly valid, it is just not a standard-conforming hash.
For SHAKE, producing piecewise output is arguably more useful and follows the same
approach, using `shake_128_sponge` (or `shake_256_sponge`) to set up the sponge.

### SIMD

`Keccak.jl` allows exploiting SIMD capabilities to simultaneously hash multiple messages,
which however have to be of the same size:
```julia-repl
julia> msg1 = rand(UInt8, 1_000_000); msg2 = rand(UInt8, 1_000_000);

julia> sha3_256(msg1, msg2) == (sha3_256(msg1), sha3_256(msg2))
true
```
The advantage is increased performance, in particular if the number of messages fits the
hardware, i.e. usually two or four:
```julia-repl
julia> @b sha3_256($msg1, $msg2), (sha3_256($msg1), sha3_256($msg2))
(2.709 ms, 3.612 ms)

julia> @b sha3_256($msg1, $msg2, $msg2), (sha3_256($msg1), sha3_256($msg2), sha3_256($msg2))
(3.678 ms, 4.880 ms)

julia> @b sha3_256($msg1, $msg2, $msg2, $msg2), (sha3_256($msg1), sha3_256($msg2), sha3_256($msg2), sha3_256($msg2))
(3.177 ms, 7.260 ms)
```

To create a sponge suitable for SIMD operation, pass a `Val` with the desired
SIMD-multiplicity when constructing the sponge. Then, `absorb` accepts multiple data
arguments (or just one reused for every data path), and `squeeze` returns a tuple of
(partial) hashes in addition to the sponge. An example inspired by ML-KEM could look like
this:
```julia
sponge = shake_128_sponge(Val(4)); # 4-fold SIMD
sponge = absorb(sponge, ρ) # ρ holds 32 pseudo-random bytes
sponge = absorb(sponge, (0x00, 0x00), (0x00, 0x01), (0x01, 0x00), (0x01, 0x01))
sponge = pad(sponge)
while some_condition
    sponge, (C00, C01, C10, C11) = squeeze(sponge, Val(3))
    # ... do something with the Cs, each holding three bytes ...
end
```

## Further Keccak-based hashes

The hashing functions described in
[SHA-3 Derived Functions: cSHAKE, KMAC, TupleHash and ParallelHash][NIST-SP-800-185]
are available as
* `cshake_128`, `cshake_128_sponge`, `cshake_256`, `cshake_256_sponge`
* `kmac_128`, `kmac_xof_128`, `kmac_128_sponge`, `kmac_256`, `kmac_xof_256`, `kmac_256_sponge`
* `tuplehash_128`, `tuplehash_xof_128`, `tuplehash_256`, `tuplehash_xof_256`
* `parallelhash_128`, `parallelhash_xof_128`, `parallelhash_256`, `parallelhash_xof_256`

See the respective docstrings for further information.

[NIST-SP-800-185]: https://doi.org/10.6028/NIST.SP.800-185
[NIST-FIPS-202]: https://doi.org/10.6028/NIST.FIPS.202
