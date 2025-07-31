const TupleHash_N = Tuple(codeunits("TupleHash"))

for d in [128, 256]
    R = (1600-2*d) ÷ 64
    cshake_spongefunc = Symbol("cshake_$(d)_sponge")
    tuplehashhelperfunc = Symbol("tuplehashhelper_$(d)")
    tuplehashfunc = Symbol("tuplehash_$(d)")
    tuplehashxoffunc = Symbol("tuplehash_xof_$(d)")
    @eval begin
        $(tuplehashhelperfunc)(data, ::Val{L}, S) where {L} =
            $(tuplehashhelperfunc)(data, L, S)
        function $(tuplehashhelperfunc)(data, len::Integer, S)
            sponge = $(cshake_spongefunc)(TupleHash_N, S)
            sponge = foldl(absorb_encoded_string, data, init = sponge)
            sponge = absorb_right_encoded(sponge, 8 * len)
            return pad(sponge)
        end

        """
            $($(tuplehashfunc))(data, len, S=())

        Computes the TupleHash$($(d)) of length `len` from `data` and the optional
        customization string `S`.

        The `data` can be any iterable object that produces `Tuple`s or `AbstractVector`s
        of `UInt8`s or `String`s. (I.e. contrary to what the name suggests, `data` being a
        `Tuple` is just one option).

        The hash is returned as a tuple of `UInt8`s if `len` is a `Val`, as a
        `Vector{UInt8}` otherwise.

        # Example
        ```jldoctest; setup = :(using Keccak)
        julia> $($(tuplehashfunc))(["abcd", "ef"], Val(16))
        $($(tuplehashfunc)(["abcd", "ef"], Val(16)))

        julia> $($(tuplehashfunc))(["abc", "def"], Val(16)) # different partitioning of same data -> different hash
        $($(tuplehashfunc)(["abc", "def"], Val(16)))

        julia> $($(tuplehashfunc))(["abc", "def"], Val(12)) # different output length -> different hash
        $($(tuplehashfunc)(["abc", "def"], Val(12)))
        ```
        """
        function $(tuplehashfunc)(data, len::Union{Val, Integer}, S::AbsorbableData = ())
            sponge = $(tuplehashhelperfunc)(data, len, S)
            return squeeze(sponge, len)[2]
        end

        """
            $($(tuplehashxoffunc))(data, S=())

        Creates a sponge for the arbitrary-length output TupleHash$($(d)) for `data`
        and the optional customization string `S`.

        The `data` can be any iterable object that produces `Tuple`s or `AbstractVector`s
        of `UInt8`s or `String`s. (I.e. contrary to what the name suggests, `data` being a
        `Tuple` is just one option).

        The returned sponge is ready to squeeze output of the desired length from.
        """
        $(tuplehashxoffunc)(data, S::AbsorbableData = ()) = $(tuplehashhelperfunc)(data, 0, S)

        """
            $($(tuplehashxoffunc))(data, len, S=())

        Computes the arbitrary-length output TupleHash$($(d)) of length `len` from `data`
        and the optional customization string `S`.

        The `data` can be any iterable object that produces `Tuple`s or `AbstractVector`s
        of `UInt8`s or `String`s. (I.e. contrary to what the name suggests, `data` being a
        `Tuple` is just one option).

        The hash is returned as a tuple of `UInt8`s if `len` is a `Val`, as a
        `Vector{UInt8}` otherwise.

        The difference to `$($(tuplehashfunc))` is that the `len` is not itself hashed.

        # Example
        ```jldoctest; setup = :(using Keccak)
        julia> $($(tuplehashxoffunc))(["abc", "def"], Val(16))
        $($(tuplehashxoffunc)(["abc", "def"], Val(16)))

        julia> $($(tuplehashxoffunc))(["abc", "def"], Val(12)) # prefix of the above
        $($(tuplehashxoffunc)(["abc", "def"], Val(16))[1:12])
        ```
        """
        function $(tuplehashxoffunc)(data, len::Union{Val,Integer}, S::AbsorbableData = ())
            return squeeze($(tuplehashxoffunc)(data, S), len)[2]
        end
    end
end
