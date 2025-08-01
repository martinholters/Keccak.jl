# Internals

API described here is intended for internal use.
Changing it in incompatible ways will not be considered a breaking change w.r.t versioning.
If you feel the need to use any of this outside of Keccak.jl, feel free to open a Github
issue about it.

## Sponge
```@docs
Keccak.Sponge
Keccak.AbsorbableData
Keccak.update
Keccak.lanetype
Keccak.rate
```

## Kᴇᴄᴄᴀᴋ Helpers
```@docs
Keccak.ℓ
Keccak.rc
Keccak.round_consts
```
## Kᴇᴄᴄᴀᴋ step mappings
These functions follow the specification in section 3.2 of FIPS-202.
```@docs
Keccak.ι
Keccak.θ
Keccak.π
Keccak.χ
Keccak.ρ
```

## SP.800-185 helpers
```@docs
Keccak.absorb_right_encoded
Keccak.absorb_left_encoded
Keccak.absorb_encoded_string
Keccak.absorb_ratepadded
Keccak.KMACPad
```
