for d in [224, 256, 384, 512]
    R = (1600-2*d) ÷ 64
    spongefunc = Symbol("sha3_$(d)_sponge")
    shafunc = Symbol("sha3_$(d)")
    @eval begin
        """
            $($spongefunc)()

        Returns a sponge suitable for SHA3-$($d) computation.
        """
        $spongefunc() = KeccakSponge{$R,UInt64}(KeccakPad(0x06))

        """
            $($spongefunc)(::Val{N})

        Returns a sponge suitable for SHA3-$($d) computation using `N`-fold SIMD, i.e. a
        sponge with `SIMD.Vec{N}` contents.
        """
        $spongefunc(::Val{N}) where {N} = KeccakSponge{$R,Vec{N,UInt64}}(KeccakPad(0x06))

        """
            $($shafunc)(data)

        Computes the SHA3-$($d) hash of `data`.

        The provided `data` must be an `AbstractVector` or a `Tuple` of `UInt8`s, or a
        `String` (which is converted to bytes by `codeunits`).
        The hash is returned as a $($(d÷8))-tuple of `UInt8`s.

        # Examples

        ```julia-repl; setup = :(using Keccak)
        julia> $($shafunc)(fill(0xa3, 200))
        $($shafunc(fill(0xa3, 200)))
        ```
        """
        function $shafunc(data::AbsorbableData)
            sponge = $spongefunc()
            sponge = absorb(sponge, data)
            sponge = pad(sponge)
            return squeeze(sponge, Val($(d÷8)))[2]
        end

        """
            $($shafunc)(data1, ..., dataN)

        Computes the SHA3-$($d) hashes of `data1`, ..., `dataN` (for `N >= 2`).

        The provided data must be `AbstractVector`s or `Tuple`s of `UInt8`s or `Strings`,
        and they all must have the same length in bytes (using `codeunits` for `String`s)
        as SIMD is used to speed up the hash computation.
        The hashes are returned as an `N`-tuple of $($(d÷8))-tuples of `UInt8`s.
        """
        function $shafunc(data::Vararg{AbsorbableData,N}) where {N}
            if !allequal(map(absorblength, data))
                throw(DimensionMismatch("data arguments provided to `$($shafunc)` have different length"))
            end
            sponge = $spongefunc(Val(N))
            sponge = absorb(sponge, data...)
            sponge = pad(sponge)
            return squeeze(sponge, Val($(d÷8)))[2]
        end
    end
end
