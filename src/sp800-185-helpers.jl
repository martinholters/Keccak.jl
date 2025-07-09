"""
    absorb_right_encoded(sponge::Sponge, x::Integer)

Absorbs the non-negative integer `x` into the given `sponge` after applying the
"right_encode" operation from [NIST SP.800-185](https://doi.org/10.6028/NIST.SP.800-185)
to it.

!!! note
    The current implementation requires `x ≤ typemax(UInt64)`.
"""
absorb_right_encoded(sponge::Sponge, x::Integer) =
    absorb_right_encoded(sponge, convert(UInt64, x))
function absorb_right_encoded(sponge::Sponge, x::UInt64)
    x = reinterpret(NTuple{8, UInt8}, x)
    if x[8] != 0
        sponge = absorb(sponge, (reverse(x)..., 0x08))
    elseif x[7] != 0
        sponge = absorb(sponge, (reverse(x[1:7])..., 0x07))
    elseif x[6] != 0
        sponge = absorb(sponge, (reverse(x[1:6])..., 0x06))
    elseif x[5] != 0
        sponge = absorb(sponge, (reverse(x[1:5])..., 0x05))
    elseif x[4] != 0
        sponge = absorb(sponge, (reverse(x[1:4])..., 0x04))
    elseif x[3] != 0
        sponge = absorb(sponge, (reverse(x[1:3])..., 0x03))
    elseif x[2] != 0
        sponge = absorb(sponge, (reverse(x[1:2])..., 0x02))
    else
        sponge = absorb(sponge, (x[1], 0x01))
    end
    return sponge
end

"""
    absorb_left_encoded(sponge::Sponge, x::Integer)

Absorbs the non-negative integer `x` into the given `sponge` after applying the
"left_encode" operation from [NIST SP.800-185](https://doi.org/10.6028/NIST.SP.800-185)
to it.
Returns the updated sponge.

!!! note
    The current implementation requires `x ≤ typemax(UInt64)`.
"""
absorb_left_encoded(sponge::Sponge, x::Integer) =
    absorb_left_encoded(sponge, convert(UInt64, x))
function absorb_left_encoded(sponge::Sponge, x::UInt64)
    x = reinterpret(NTuple{8, UInt8}, x)
    if x[8] != 0
        sponge = absorb(sponge, (0x08, reverse(x)...))
    elseif x[7] != 0
        sponge = absorb(sponge, (0x07, reverse(x[1:7])...))
    elseif x[6] != 0
        sponge = absorb(sponge, (0x06, reverse(x[1:6])...))
    elseif x[5] != 0
        sponge = absorb(sponge, (0x05, reverse(x[1:5])...))
    elseif x[4] != 0
        sponge = absorb(sponge, (0x04, reverse(x[1:4])...))
    elseif x[3] != 0
        sponge = absorb(sponge, (0x03, reverse(x[1:3])...))
    elseif x[2] != 0
        sponge = absorb(sponge, (0x02, reverse(x[1:2])...))
    else
        sponge = absorb(sponge, (0x01, x[1]))
    end
    return sponge
end

"""
    absorb_encoded_string(sponge::Sponge, str::AbsorbableData)

Absorbs the string `str` (which can also be a `Tuple` or `AbstractVector` of `UInt8`s) into
the given `sponge` after applying the "encode_string" operation from
[NIST SP.800-185](https://doi.org/10.6028/NIST.SP.800-185) to it.
Returns the updated sponge.
"""
function absorb_encoded_string(sponge::Sponge, str::AbsorbableData)
    sponge = absorb_left_encoded(sponge, 8*absorblength(str))
    return absorb(sponge, str)
end

"""
    absorb_ratepadded(f, sponge::Sponge)

Absorb data into the given sponge and zero-pad to a multiple of the sponges rate.
The padding uses the "bytepad" operation from
[NIST SP.800-185](https://doi.org/10.6028/NIST.SP.800-185), but is limited to "w" being the
sponge's rate.

The data is provided implicitly by passing a function `f` that is called with a sponge and
has to return an updated sponge with the data absorbed.

# Example
```jldoctest; setup=:(using Keccak: absorb, absorb_encoded_string, absorb_ratepadded, rate, sha3_256_sponge; sponge = sha3_256_sponge())
absorb_ratepadded(sponge) do sponge
    return absorb_encoded_string(sponge, "Test")
end == absorb(sponge, [
    0x01, UInt8(rate(sponge)), # left_encode'd pad width w=rate(sponge)
    # start of data produced by `f` (i.e. absorb_encoded_string)
    0x01, UInt8(4*8),
    codeunits("Test")...,
    # end of data produced by `f`
    zeros(UInt8, rate(sponge)-8)... # absorbed 8 bytes so far, zero-pad up to rate
])

# output

true
```
"""
function absorb_ratepadded(f::F, sponge::Sponge) where F
    old_k = sponge.k # we pad to a multiple of the rate, so overall, k mustn't change
    sponge = absorb_left_encoded(sponge, rate(sponge))
    sponge = f(sponge)
    if sponge.k <= old_k
        # implicitly zero-pad the input by just progressing k
        return update(sponge, sponge.state, old_k)
    else
        # zero-padding wraps around, so need to transform once
        return update(sponge, sponge.transform(sponge.state), old_k)
    end
end
