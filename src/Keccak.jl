module Keccak

if reinterpret(UInt64, (0x08, 0x07, 0x06, 0x05, 0x04, 0x03, 0x02, 0x01)) != 0x0102030405060708
    # Julia itself at present only supports little-endian architectures, so this should never happen
    error("Keccak.jl requires a little-endian architecture")
end

using SIMD: Vec

export KeccakP, KeccakPad, KeccakSponge, absorb, keccak_p, pad, squeeze
export sha3_224, sha3_224_sponge, sha3_256, sha3_256_sponge, sha3_384, sha3_384_sponge,
    sha3_512, sha3_512_sponge
export shake_128, shake_128_sponge, shake_256, shake_256_sponge
export cshake_128, cshake_128_sponge, cshake_256, cshake_256_sponge

include("sponge.jl")
include("keccakp.jl")
include("keccaksponge.jl")
include("sha3.jl")
include("shake.jl")
include("sp800-185-helpers.jl")
include("cshake.jl")

end # module Keccak
