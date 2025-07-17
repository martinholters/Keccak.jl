# Keccak.jl Documentation

```@docs
sha3_224
sha3_256
sha3_384
sha3_512
shake_128
shake_256
sha3_224_sponge
sha3_256_sponge
sha3_384_sponge
sha3_512_sponge
shake_128_sponge
shake_256_sponge
cshake_128
cshake_256
cshake_128_sponge
cshake_256_sponge
kmac_128
kmac_256
kmac_xof_128
kmac_xof_256
kmac_128_sponge
kmac_256_sponge
tuplehash_128
tuplehash_256
tuplehash_xof_128
tuplehash_xof_256
parallelhash_128
parallelhash_256
parallelhash_xof_128
parallelhash_xof_256
absorb
pad
squeeze
keccak_p
KeccakP
KeccakPad
KeccakSponge
```

## Internals
### Sponge
```@docs
Keccak.Sponge
Keccak.AbsorbableData
Keccak.update
Keccak.rate
```

### Kᴇᴄᴄᴀᴋ Helpers
```@docs
Keccak.ℓ
Keccak.rc
Keccak.round_consts
```
### Kᴇᴄᴄᴀᴋ step mappings
These functions follow the specification in section 3.2 of FIPS-202.
```@docs
Keccak.ι
Keccak.θ
Keccak.π
Keccak.χ
Keccak.ρ
```

### SP.800-185 helpers
```@docs
Keccak.absorb_right_encoded
Keccak.absorb_left_encoded
Keccak.absorb_encoded_string
Keccak.absorb_ratepadded
Keccak.KMACPad
```
