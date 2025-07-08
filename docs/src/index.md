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
