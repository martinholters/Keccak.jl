# Keccak.jl Documentation

```@docs
absorb
squeeze
keccak_p
KeccakP
KeccakSponge
```

## Internals
### Sponge
```@docs
Keccak.Sponge
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
