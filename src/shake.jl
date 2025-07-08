for d in [128, 256]
    R = (1600-2*d) ÷ 64
    spongefunc = Symbol("shake_$(d)_sponge")
    shakefunc = Symbol("shake_$(d)")
    @eval begin
        """
            $($spongefunc)()

        Returns a sponge suitable for SHAKE-$($d) computation.
        """
        $spongefunc() = KeccakSponge{$R,UInt64}(KeccakPad(0x1f))

        """
            $($spongefunc)(::Val{N})

        Returns a sponge suitable for SHAKE-$($d) computation using `N`-fold SIMD, i.e. a
        sponge with `SIMD.Vec{N}` contents.
        """
        $spongefunc(::Val{N}) where {N} = KeccakSponge{$R,Vec{N,UInt64}}(KeccakPad(0x1f))

        """
            $($shakefunc)(data)

        Returns a SHAKE-$($d) sponge with `data` absorbed, ready to `squeeze` from.

        The provided `data` must be an `AbstractVector` or a `Tuple` of `UInt8`s, or a
        `String` (which is converted to bytes by `codeunits`).
        Appropriate padding is applied by `$($shakefunc)`.
        """
        function $shakefunc(data::AbsorbableData)
            sponge = $spongefunc()
            sponge = absorb(sponge, data)
            return pad(sponge)
        end

        """
            $($shakefunc)(data, len)

        Computes the SHAKE-$($d) hash of `data`.

        The provided `data` must be an `AbstractVector` or a `Tuple` of `UInt8`s, or a
        `String` (which is converted to bytes by `codeunits`).
        The hash is returned as a tuple of `UInt8`s if `len` is a `Val`, as a
        `Vector{UInt8}` otherwise.

        # Examples

        ```julia-repl; setup = :(using Keccak)
        julia> $($shakefunc)(fill(0xa3, 200), Val(20))
        $($shakefunc(fill(0xa3, 200), Val(20)))
        ```
        """
        function $shakefunc(data::AbsorbableData, len::Union{Val,Integer})
            sponge = $spongefunc()
            sponge = absorb(sponge, data)
            sponge = pad(sponge)
            return squeeze(sponge, len)[2]
        end

        """
            $($shakefunc)(data1, ..., dataN)

        Absorbs `data1`, ..., `dataN` (for `N >= 2`) into a new `N`-fold-SIMD SHAKE-$($d)
        sponge and `pad`s it.

        The provided data must be `AbstractVector`s or `Tuple`s of `UInt8`s or `Strings`,
        and they all must have the same length in bytes (using `codeunits` for `String`s)
        as SIMD is used to speed up the hash computation.
        The data is absorbed into a fresh (SIMD) SHAKE-$($d) sponge and padding is applied.

        The returned sponge is ready for being `squeeze`d.
        """
        function $shakefunc(data::Vararg{AbsorbableData,N}) where {N}
            if !allequal(map(absorblength, data))
                throw(DimensionMismatch("data arguments provided to `$($shakefunc)` have different length"))
            end
            sponge = $spongefunc(Val(N))
            sponge = absorb(sponge, data...)
            return pad(sponge)
        end

        """
            $($shakefunc)(data1, ..., dataN, len)

        Computes the SHAKE-$($d) hash of `data1`, ..., `dataN` (for `N >= 2`).

        The provided data must be `AbstractVector`s or `Tuple`s of `UInt8`s or `Strings`,
        and they all must have the same length in bytes (using `codeunits` for `String`s)
        as SIMD is used to speed up the hash computation.

        An `N`-tuple of the individual hashes of the desired length is returned. The hashes
        are tuples of `UInt8`s if `len` is a `Val`, `Vector{UInt8}`s otherwise.
        """
        function $shakefunc(args::Vararg{Any,N}) where N
            data = args[1:end-1]::NTuple{N-1,AbsorbableData}
            len = args[end]
            sponge = $shakefunc(data...)
            return squeeze(sponge, len)[2]
        end
    end
end
