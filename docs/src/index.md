# Keccak.jl Documentation

## Why?

The two most important, i.e. most commonly used, algorithms this package provides, namely
SHA-3 and SHAKE, are available is the `SHA` stdlib as well. So why use `Keccak`?

* Flexibility. Not only does `Keccak` offer
  [more readily available algorithms](@ref "SHA-3 Derived Functions"), it also provides the
  [building blocks](@ref "Keccak-based sponge construction") to implement further
  Keccak-based algorithms outside of the package. (Also, it allows things like turning
  SHA-3 into a SHAKE-like PRF by producing more output data from it than intended by the
  standard. One might consider this a bit too much flexibility, though.)
* Performance. The `Keccak` implementations are (at present) faster and do not allocate.
  While the run-time of the `SHA` implementations might be improved in the future, the fact
  that they allocate memory is to some extent baked into the API and hence, not likely to
  change. Furthermore, in admittedly limited application scenarios, `Keccak` allows SIMD
  processing to further speed up simultaneous computation of multiple hashes.

## For the impatient

Great, more performance is always a good thing. I'd like to replace `SHA` with `Keccak`
for computing SHA-3 hashes. How do I do it?

Assuming you have added `Keccak` to your project using `Pkg`, your ready to compute your
first SHA-3 digest:
```jldoctest
julia> using Keccak: sha3_256

julia> sha3_256("a message to compute the digest of")
(0xfa, 0x24, 0x49, 0x91, 0xa1, 0x94, 0x44, 0x1e, 0xd6, 0xa8, 0x34, 0x19, 0xd2, 0x00, 0x44, 0x12, 0x15, 0xa9, 0x7c, 0x49, 0xc6, 0xcf, 0x5e, 0x28, 0x7f, 0x75, 0x75, 0x45, 0x73, 0x6d, 0x25, 0x23)
```
In fact, you can use [`sha3_224`](@ref), [`sha3_256`](@ref), [`sha3_384`](@ref), and
[`sha3_512`](@ref) almost as drop-in replacements for their counterparts from `SHA`,
paying attention to the following:
* `SHA` and `Keccak` export functions with the same names, so watch out for name clashes
  when using both. (Explicit `using` states as above are recommended anyway.)
* The digest is returned as a tuple rather than a vector. Ideally, you can process the
  tuple like you processed the vector from `SHA`. If needed, you can `collect` the tuple
  into a vector (at the cost of memory allocation, of course).
* Directly hashing from `IO` objects is not supported (yet).

## Going further

To get to know everything offered by `Keccak` it is recommended to start reading the
[user guide from the beginning](@ref "SHA-3 (and SHAKE)"), as its sections not only build
on each other, but are also (more or less) ordered by importance.