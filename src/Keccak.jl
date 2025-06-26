module Keccak

if reinterpret(UInt64, (0x08, 0x07, 0x06, 0x05, 0x04, 0x03, 0x02, 0x01)) != 0x0102030405060708
    # Julia itself at present only supports little-endian architectures, so this should never happen
    error("Keccak.jl requires a little-endian architecture")
end

using SIMD: Vec

export KeccakP, KeccakSponge, absorb, keccak_p, squeeze

include("sponge.jl")
include("keccakp.jl")

end # module Keccak
