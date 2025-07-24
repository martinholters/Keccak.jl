const KMAC_N = Tuple(codeunits("KMAC")) # cSHAKE function string for KMAC

"""
    KMACPad(basepad::KeccakPad, len::UInt64)

A padding function that first appends (i.e. absorbs) `8*len` after "right_encode"ing it and
then invokes `basepad`.
"""
struct KMACPad
    basepad::KeccakPad
    len::UInt64
end

function (kmacpad::KMACPad)(sponge)
    sponge = absorb_right_encoded(sponge, 8*kmacpad.len)
    return kmacpad.basepad(sponge)
end

for d in [128, 256]
    R = (1600-2*d) ÷ 64
    cshake_spongefunc = Symbol("cshake_$(d)_sponge")
    kmac_spongefunc = Symbol("kmac_$(d)_sponge")
    kmacfunc = Symbol("kmac_$(d)")
    kmacxoffunc = Symbol("kmac_xof_$(d)")
    @eval begin
        function $(kmacxoffunc) end

        """
            $($(kmac_spongefunc))(K, len=0, S=())

        Creates a sponge for computing the Kᴇᴄᴄᴀᴋ Message Authentication Code (KMAC$($(d)))
        using key `K` and (optional) customization string `S`. The length `len` (in bytes)
        is absorbed into the sponge, so the output depends on it, though it does not
        prescribe the amount of data to `squeeze`. However, to compute a
        standard-conforming KAMC$($(d)), the lengths must match, or `len` must be set to zero
        to indicate an arbitrary output length KMACXOF$($(d)).
        """
        function $(kmac_spongefunc)(K::AbsorbableData, len::Integer, S::AbsorbableData = ())
            sponge = $cshake_spongefunc(KMAC_N, S)
            sponge = absorb_ratepadded(sponge) do sponge
                return absorb_encoded_string(sponge, K)
            end
            return Sponge{$R}(
                sponge.transform,
                KMACPad(sponge.pad, len),
                sponge.state,
                sponge.k,
            )
        end
        $(kmac_spongefunc)(K::AbsorbableData, S::AbsorbableData = ()) =
            $(kmac_spongefunc)(K, 0, S)

        """
            $($(kmacfunc))(K, data, len, S=())

        Computes the Kᴇᴄᴄᴀᴋ Message Authentication Code (KMAC$($(d))) of `data` using
        key `K`, output length `len` (in bytes) and (optional) customization string `S`.

        If `len` is provided as a `Val`, the KMAC is returned as a `Tuple`, otherwise as a
        `Vector`.

        !!! note
            The output length is included in the hash, so changing it will not only produce
            more or fewer bytes of the same hash, but rather change the hash completely.
            See `$($(kmacxoffunc))` for an alternative.
        """
        function $(kmacfunc)(
            K::AbsorbableData,
            data::AbsorbableData,
            len::Union{Val{L},Integer},
            S::AbsorbableData = (),
        ) where {L}
            intlen = len isa Val ? L : len
            sponge = $(kmac_spongefunc)(K, intlen, S)
            sponge = absorb(sponge, data)
            sponge = pad(sponge)
            return squeeze(sponge, len)[2]
        end

        """
            $($(kmacxoffunc))(K, data, S=())

        Creates a sponge for the extensible output Kᴇᴄᴄᴀᴋ Message Authentication Code
        (KMACXOF$($(d))) of `data` using key `K` and (optional) customization string `S`.

        The sponge is ready to squeeze output of the desired length from.
        """
        function $(kmacxoffunc)(
            K::AbsorbableData,
            data::AbsorbableData,
            S::AbsorbableData = (),
        )
            sponge = $(kmac_spongefunc)(K, 0, S)
            sponge = absorb(sponge, data)
            return pad(sponge)
        end

        """
            $($(kmacxoffunc))(K, data, len, S=())

        Computes the extensible output Kᴇᴄᴄᴀᴋ Message Authentication Code (KMACXOF$($(d)))
        of `data` using key `K`, output length `len` (in bytes) and (optional)
        customization string `S`.

        If `len` is provided as a `Val`, the KMAC is returned as a `Tuple`, otherwise as a
        `Vector`. Contrary to `$($(kmacfunc))`, the output length is not included in the
        hash, so the result for a smaller `len` is the prefix of the result for a larger
        `len`.
        """
        function $(kmacxoffunc)(
            K::AbsorbableData,
            data::AbsorbableData,
            len::Union{Val,Integer},
            S::AbsorbableData = (),
        )
            return squeeze($(kmacxoffunc)(K, data, S), len)[2]
        end
    end
end
