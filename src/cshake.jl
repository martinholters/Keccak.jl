for d in [128, 256]
    R = (1600-2*d) ÷ 64
    cshake_spongefunc = Symbol("cshake_$(d)_sponge")
    cshakefunc = Symbol("cshake_$(d)")
    @eval begin
        """
            $($cshake_spongefunc)(N=(), S=())

        Returns a sponge suitable for cSHAKE-$($d) computation for function-name `N` and
        customization `S`. Both `N` and `S` can be `Tuple`s or `AbstractVector`s of
        `UInt8`, or `String` (of which their `codeunits` are used) and default to `()`,
        i.e. the empty string.
        """
        function $(cshake_spongefunc)(N::AbsorbableData=(), S::AbsorbableData=())
            if isempty(N) && isempty(S)
                return $(Symbol("shake_$(d)_sponge"))()
            end
            sponge = KeccakSponge{$R,UInt64}(KeccakPad(0x04))
            sponge = absorb_ratepadded(sponge) do sponge
                sponge = absorb_encoded_string(sponge, N)
                return absorb_encoded_string(sponge, S)
            end
            return sponge
        end

        """
            $($cshakefunc)(data, N=(), S=())

        Returns a cSHAKE-$($d) sponge for function-name `N` and
        customization `S` with `data` absorbed, ready to `squeeze` from.

        The provided `data`, `N`, and `S` must be `AbstractVector`s or a `Tuple`s of
        `UInt8`s, or a `String` (which is converted to bytes by `codeunits`).
        Appropriate padding is applied by `$($cshakefunc)`.
        """
        function $cshakefunc(data::AbsorbableData, N::AbsorbableData=(), S::AbsorbableData=())
            sponge = $cshake_spongefunc(N, S)
            sponge = absorb(sponge, data)
            return pad(sponge)
        end

        """
            $($cshakefunc)(data, len, N=(), S=())

        Computes the cSHAKE-$($d) hash of `data` for function-name `N` and
        customization `S`.

        The provided `data`, `N`, and `S` must be `AbstractVector`s or a `Tuple`s of
        `UInt8`s, or a `String` (which is converted to bytes by `codeunits`).
        The hash is returned as a tuple of `UInt8`s if `len` is a `Val`, as a
        `Vector{UInt8}` otherwise.
        """
        function $cshakefunc(
            data::AbsorbableData,
            len::Union{Val,Integer},
            N::AbsorbableData=(),
            S::AbsorbableData=(),
        )
            sponge = $cshake_spongefunc(N, S)
            sponge = absorb(sponge, data)
            sponge = pad(sponge)
            return squeeze(sponge, len)[2]
        end
    end
end
