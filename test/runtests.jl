using Keccak
using OffsetArrays: Origin
using SIMD: Vec
using Test: @test, @test_throws, @testset

@testset "Sponge" begin
    @testset "identity sponge" for (maybe_tup, maybe_val) in ((identity, identity), (Tuple, Val))
        sponge = Keccak.Sponge{3}(identity, identity, Tuple(zeros(UInt16,7)), 0)
        sponge = Keccak.absorb(sponge, maybe_tup([0x01, 0x02, 0x03, 0x04, 0x05, 0x06]))
        @test sponge.state == (0x0201, 0x0403, 0x0605, 0x0000, 0x0000, 0x0000, 0x0000)
        sponge = Keccak.absorb(sponge, maybe_tup([0x07, 0x08, 0x09]))
        @test sponge.state == (
            0x0201 ⊻ 0x0807, 0x0403 ⊻ 0x0009, 0x0605,
            0x0000, 0x0000, 0x0000, 0x0000
        )
        sponge = Keccak.absorb(sponge, maybe_tup([0x0a, 0x0b, 0x0c, 0x0d]))
        @test sponge.state == (
            0x0201 ⊻ 0x0807 ⊻ 0x000d, 0x0403 ⊻ 0x0a09, 0x0605 ⊻ 0x0c0b,
            0x0000, 0x0000, 0x0000, 0x0000
        )
        sponge = Keccak.absorb(sponge, maybe_tup(zeros(UInt8, 5)))
        ref_output = [0x01 ⊻ 0x07 ⊻ 0x0d, 0x02 ⊻ 0x08, 0x03 ⊻ 0x09, 0x04 ⊻ 0x0a, 0x05 ⊻ 0x0b, 0x06 ⊻ 0x0c]
        i = 1
        for l in 0:20
            sponge, output = Keccak.squeeze(sponge, maybe_val(l))
            @test length(output) == l
            for j in 1:l
                @test output[j] == ref_output[i]
                i = mod1(i+1, 6)
            end
        end
    end
    @testset "rot sponge" for (maybe_tup, maybe_val) in ((identity, identity), (Tuple, Val))
        f(state) = (state[end], state[1:end-1]...)
        sponge = Keccak.Sponge{3}(f, identity, Tuple(zeros(UInt16,7)), 0)
        sponge = Keccak.absorb(sponge, maybe_tup([0x01, 0x02, 0x03, 0x04, 0x05, 0x06]))
        @test sponge.state == (0x0000, 0x0201, 0x0403, 0x0605, 0x0000, 0x0000, 0x0000)
        sponge = Keccak.absorb(sponge, maybe_tup([0x07, 0x08, 0x09]))
        @test sponge.state == (
            0x0807, 0x0201 ⊻ 0x0009, 0x0403,
            0x0605, 0x0000, 0x0000, 0x0000
        )
        sponge = Keccak.absorb(sponge, maybe_tup([0x0a, 0x0b, 0x0c, 0x0d]))
        @test sponge.state == (
            0x000d, 0x0807, 0x0201 ⊻ 0x0a09, 0x0403 ⊻ 0x0c0b,
            0x0605, 0x0000, 0x0000
        )
        sponge = Keccak.absorb(sponge, maybe_tup(zeros(UInt8, 5)))
        ref_output = [
            0x00, 0x00, 0x0d, 0x00, 0x07, 0x08, 0x01 ⊻ 0x09, 0x02 ⊻ 0x0a,
            0x03 ⊻ 0x0b, 0x04 ⊻ 0x0c, 0x05, 0x06, 0x00, 0x00
        ]
        i = 1
        for l in 0:20
            sponge, output = Keccak.squeeze(sponge, maybe_val(l))
            @test length(output) == l
            for j in 1:l
                @test output[j] == ref_output[i]
                i = mod1(i+1, 6)
                if i == 1
                    circshift!(ref_output, 2)
                end
            end
        end
    end

    @testset "$N-fold SIMD rot sponge" for N in (1,2,4)
        # compare an N-SIMD sponge to N scalar sponges on random data
        f(state) = (state[end], state[1:end-1]...)
        ref_sponges = [Keccak.Sponge{3}(f, identity, Tuple(zeros(UInt16,7)), 0) for _ in 1:N]
        simd_sponge = Keccak.Sponge{3}(f, identity, Tuple(zeros(Vec{N,UInt16},7)), 0)
        for _ in 1:100
            # single (same) message for all N instances
            len = rand(0:8)
            msg = rand(UInt8, len)
            if rand(Bool)
                msg = Tuple(msg)
            end
            ref_sponges = map(sponge -> Keccak.absorb(sponge, msg), ref_sponges)
            simd_sponge = Keccak.absorb(simd_sponge, msg)
            for n in 1:N
                @test map(v -> v[n], simd_sponge.state) == ref_sponges[n].state
            end
        end
        for _ in 1:100
            # different message per instance
            len = rand(0:8)
            msgs = [rand(UInt8, len) for _ in 1:N]
            if rand(Bool)
                msgs = Tuple.(msgs)
            end
            ref_sponges = map((sponge, msg) -> Keccak.absorb(sponge, msg), ref_sponges, msgs)
            simd_sponge = Keccak.absorb(simd_sponge, msgs...)
            for n in 1:N
                @test map(v -> v[n], simd_sponge.state) == ref_sponges[n].state
            end
        end
        if N > 1
            @test_throws DimensionMismatch Keccak.absorb(simd_sponge, [rand(UInt8, len) for len in 1:N]...)
        end
        # padding
        len = 6 - simd_sponge.k
        msgs = [zeros(UInt8, len) for _ in 1:N]
        if rand(Bool)
            msgs = Tuple.(msgs)
        end
        ref_sponges = map((sponge, msg) -> Keccak.absorb(sponge, msg), ref_sponges, msgs)
        simd_sponge = Keccak.absorb(simd_sponge, msgs...)
        for n in 1:N
            @test map(v -> v[n], simd_sponge.state) == ref_sponges[n].state
        end
        # verify squeezed output
        for _ in 1:100
            len = rand(0:20)
            simd_sponge, simd_data = Keccak.squeeze(simd_sponge, len)
            for n in 1:N
                ref_sponges[n], ref_data = Keccak.squeeze(ref_sponges[n], len)
                @test simd_data[n] == ref_data
            end
        end
        for _ in 1:100
            len = rand(0:20)
            simd_sponge, simd_data = Keccak.squeeze(simd_sponge, Val(len))
            for n in 1:N
                ref_sponges[n], ref_data = Keccak.squeeze(ref_sponges[n], Val(len))
                @test simd_data[n] == ref_data
            end
        end
    end
end

@testset "LFSR rc" begin
    @test Keccak.rc(0) == true
    R = Origin(0)(Bool[1, 0, 0, 0, 0, 0, 0, 0])
    for t in 1:254
        R = Origin(0)([false; R...])
        R[0] = R[0] ⊻ R[8]
        R[4] = R[4] ⊻ R[8]
        R[5] = R[5] ⊻ R[8]
        R[6] = R[6] ⊻ R[8]
        @test Keccak.rc(t) == R[8]
        resize!(R, 8)
    end
    for t in -100:-1
        @test Keccak.rc(t) == Keccak.rc(t+255)
    end
    for t in 255:17:1234
        @test Keccak.rc(t) == Keccak.rc(mod(t,255))
    end
end

@testset "Keccak-p[$(25w)]" for (T, w) in [
    (UInt8, 8),
    (UInt16, 16),
    (UInt32, 32),
    (UInt64, 64),
]
    l = round(Int, log2(w))
    state_to_array(state) = Origin(0)([Bool(state[5y+x+1] >> z & 1) for x in 0:4, y in 0:4, z in 0:(w-1)])

    @testset "theta" begin
        state = Tuple(rand(T, 25))
        A = state_to_array(state)
        C = similar(A, 0:4, 0:w-1)
        for x in 0:4, z in 0:w-1
            C[x,z] = A[x,0,z] .⊻ A[x,1,z] .⊻ A[x,2,z] .⊻ A[x,3,z] .⊻ A[x,4,z]
        end
        D = similar(A, 0:4, 0:w-1)
        for x in 0:4, z in 0:w-1
            D[x,z] = C[mod(x-1, 5), z] ⊻ C[mod(x+1, 5), mod(z-1, w)]
        end
        A′ = similar(A)
        for x in 0:4, y in 0:4, z in 0:w-1
            A′[x,y,z] = A[x,y,z] ⊻ D[x,z]
        end
        @test state_to_array(Keccak.θ(state)) == A′
    end

    @testset "rho" begin
        state = Tuple(rand(T, 25))
        A = state_to_array(state)
        A′ = similar(A)
        for z in 0:w-1
            A′[0,0,z] = A[0,0,z]
        end
        x, y = 1, 0
        for t in 0:23
            for z in 0:w-1
                A′[x,y,z] = A[x,y,mod(z-(t+1)*(t+2)÷2, w)]
            end
            x, y = y, mod(2x+3y, 5)
        end
        @test state_to_array(Keccak.ρ(state)) == A′
    end

    @testset "pi" begin
        state = Tuple(rand(T, 25))
        A = state_to_array(state)
        A′ = similar(A)
        for x in 0:4, y in 0:4, z in 0:w-1
            A′[x,y,z] = A[mod(x+3y, 5),x,z]
        end
        @test state_to_array(Keccak.π(state)) == A′
    end

    @testset "chi" begin
        state = Tuple(rand(T, 25))
        A = state_to_array(state)
        A′ = similar(A)
        for x in 0:4, y in 0:4, z in 0:w-1
            A′[x,y,z] = A[x,y,z] ⊻ ((A[mod(x+1, 5),y,z] ⊻ true) & A[mod(x+2, 5),y,z])
        end
        @test state_to_array(Keccak.χ(state)) == A′
    end

    @testset "round_consts" begin
        l = round(Int, log2(w))
        rounds = -30:12+2l-1
        RCs = Keccak.round_consts(T, Val(length(rounds)))
        @test RCs isa NTuple{length(rounds),T}
        for ir in rounds
            RC = Origin(0)(falses(w))
            for j in 0:l
                RC[2^j-1] = Keccak.rc(j + 7ir)
            end
            @test Origin(0)([Bool(RCs[ir-rounds[1]+1] >> z & 1) for z in 0:(w-1)]) == RC
        end
    end

    @testset "iota" begin
        state = Tuple(rand(T, 25))
        for ir in -30:12+2l-1
            A = state_to_array(state)
            A′ = copy(A)
            RC = Origin(0)(falses(w))
            for j in 0:l
                RC[2^j-1] = Keccak.rc(j + 7ir)
            end
            for z in 0:w-1
                A′[0,0,z] = A′[0,0,z] ⊻ RC[z]
            end
            @test state_to_array(Keccak.ι(state, Keccak.round_consts(T, Val(30+12+2l))[ir+30+1])) == A′
        end
    end

    @testset "kessak_p for $nrounds rounds" for nrounds in [12 + l, 12 + 2l, 12 + 3l]
        state = Tuple(rand(T, 25))
        state′ = state
        for ir = 12+2l-nrounds:12+2l-1
            RC = zero(T)
            for j in 0:l
                RC |= T(Keccak.rc(j + 7ir)) << (2^j-1)
            end
            state′ = Keccak.ι(Keccak.χ(Keccak.π(Keccak.ρ(Keccak.θ(state′)))), RC)
        end
        @test keccak_p(state, Val(nrounds)) == state′
    end

    @testset "$n-element SIMD kessak_p" for n in [1, 2, 4, 8], nrounds in 12 .+ [l, 2l, 3l]
        # verify that a SIMD-ed kessak_p gives the same results as multiple calls to scalar
        # kessak_p
        simd_state = Tuple(Vec(rand(UInt64,n)...) for _ in 1:25)
        simd_state′ = keccak_p(simd_state, Val(nrounds))
        for i in 1:n
            state = Tuple(s[i] for s in simd_state)
            state′ = Tuple(s[i] for s in simd_state′)
            @test state′ == keccak_p(state, Val(nrounds))
        end
    end
end

@testset "KeccakSponge" begin
    shakepad = KeccakPad(0x1f)
    # test data from
    # https://csrc.nist.gov/projects/cryptographic-standards-and-guidelines/example-values
    # for SHAKE-128, 0-bit input
    empty_ref = hex2bytes("7F9C2BA4E88F827D616045507605853ED73B8093F6EFBC88EB1A6EACFA66EF263CB1EEA988004B93103CFB0AEEFD2A686E01FA4A58E8A3639CA8A1E3F9AE57E235B8CC873C23DC62B8D260169AFA2F75AB916A58D974918835D25E6A435085B2BADFD6DFAAC359A5EFBB7BCC4B59D538DF9A04302E10C8BC1CBF1A0B3A5120EA17CDA7CFAD765F5623474D368CCCA8AF0007CD9F5E4C849F167A580B14AABDEFAEE7EEF47CB0FCA9767BE1FDA69419DFB927E9DF07348B196691ABAEB580B32DEF58538B8D23F87732EA63B02B4FA0F4873360E2841928CD60DD4CEE8CC0D4C922A96188D032675C8AC850933C7AFF1533B94C834ADBB69C6115BAD4692D8619F90B0CDF8A7B9C264029AC185B70B83F2801F2F4B3F70C593EA3AEEB613A7F1B1DE33FD75081F592305F2E4526EDC09631B10958F464D889F31BA010250FDA7F1368EC2967FC84EF2AE9AFF268E0B1700AFFC6820B523A3D917135F2DFF2EE06BFE72B3124721D4A26C04E53A75E30E73A7A9C4A95D91C55D495E9F51DD0B5E9D83C6D5E8CE803AA62B8D654DB53D09B8DCFF273CDFEB573FAD8BCD45578BEC2E770D01EFDE86E721A3F7C6CCE275DABE6E2143F1AF18DA7EFDDC4C7B70B5E345DB93CC936BEA323491CCB38A388F546A9FF00DD4E1300B9B2153D2041D205B443E41B45A653F2A5C4492C1ADD544512DDA2529833462B71A41A45BE97290B6F")
    sponge = Keccak.absorb(Keccak.KeccakSponge{21,UInt64}(shakepad), ())
    sponge = Keccak.pad(sponge)
    @test Keccak.squeeze(sponge, length(empty_ref))[2] == empty_ref

    # piece-wise output
    k = 1
    while k <= length(empty_ref)
        len = min(rand(0:200), length(empty_ref)-k+1)
        sponge′, data = Keccak.squeeze(sponge, Val(len))
        @test data == Tuple(empty_ref[k:k+len-1])
        sponge, data = Keccak.squeeze(sponge, len)
        @test data == empty_ref[k:k+len-1]
        @test sponge == sponge′
        k += len
    end

    # test data from
    # https://csrc.nist.gov/projects/cryptographic-standards-and-guidelines/example-values
    # for SHAKE-128, 1600-bit input
    len1600_ref = hex2bytes("131AB8D2B594946B9C81333F9BB6E0CE75C3B93104FA3469D3917457385DA037CF232EF7164A6D1EB448C8908186AD852D3F85A5CF28DA1AB6FE3438171978467F1C05D58C7EF38C284C41F6C2221A76F12AB1C04082660250802294FB87180213FDEF5B0ECB7DF50CA1F8555BE14D32E10F6EDCDE892C09424B29F597AFC270C904556BFCB47A7D40778D390923642B3CBD0579E60908D5A000C1D08B98EF933F806445BF87F8B009BA9E94F7266122ED7AC24E5E266C42A82FA1BBEFB7B8DB0066E16A85E0493F07DF4809AEC084A593748AC3DDE5A6D7AAE1E8B6E5352B2D71EFBB47D4CAEED5E6D633805D2D323E6FD81B4684B93A2677D45E7421C2C6AEA259B855A698FD7D13477A1FE53E5A4A6197DBEC5CE95F505B520BCD9570C4A8265A7E01F89C0C002C59BFEC6CD4A5C109258953EE5EE70CD577EE217AF21FA70178F0946C9BF6CA8751793479F6B537737E40B6ED28511D8A2D7E73EB75F8DAAC912FF906E0AB955B083BAC45A8E5E9B744C8506F37E9B4E749A184B30F43EB188D855F1B70D71FF3E50C537AC1B0F8974F0FE1A6AD295BA42F6AEC74D123A7ABEDDE6E2C0711CAB36BE5ACB1A5A11A4B1DB08BA6982EFCCD716929A7741CFC63AA4435E0B69A9063E880795C3DC5EF3272E11C497A91ACF699FEFEE206227A44C9FB359FD56AC0A9A75A743CFF6862F17D7259AB075216C0699511643B6439")
    sponge = Keccak.absorb(Keccak.KeccakSponge{21,UInt64}(shakepad), fill(0xa3, 200))
    sponge = Keccak.pad(sponge)
    @test Keccak.squeeze(sponge, length(len1600_ref))[2] == len1600_ref

    # piece-wise output
    k = 1
    while k <= length(len1600_ref)
        len = min(rand(0:200), length(len1600_ref)-k+1)
        sponge′, data = Keccak.squeeze(sponge, Val(len))
        @test data == Tuple(len1600_ref[k:k+len-1])
        sponge, data = Keccak.squeeze(sponge, len)
        @test data == len1600_ref[k:k+len-1]
        @test sponge == sponge′
        k += len
    end

    # byte-wise input
    sponge = Keccak.KeccakSponge{21,UInt64}(shakepad)
    for _ in 1:200
        sponge = Keccak.absorb(sponge, tuple(0xa3))
    end
    sponge = Keccak.pad(sponge)
    @test Keccak.squeeze(sponge, length(len1600_ref))[2] == len1600_ref

    # piece-wise input
    sponge = Keccak.KeccakSponge{21,UInt64}(shakepad)
    k = 0
    while k < 200
        len = min(rand(0:10), 200-k)
        sponge = Keccak.absorb(sponge, fill(0xa3, len))
        k += len
    end
    sponge = Keccak.pad(sponge)
    @test Keccak.squeeze(sponge, length(len1600_ref))[2] == len1600_ref
end

@testset "SHA3-$d" for (d, shafunc, spongefunc, empty_ref, len1600_ref) in [
    # reference data from
    # https://csrc.nist.gov/projects/cryptographic-standards-and-guidelines/example-values
    # for input message lengths 0 bit and 1600 bit
    (
        224,
        sha3_224,
        sha3_224_sponge,
        hex2bytes("6B4E03423667DBB73B6E15454F0EB1ABD4597F9A1B078E3F5B5A6BC7"),
        hex2bytes("9376816ABA503F72F96CE7EB65AC095DEEE3BE4BF9BBC2A1CB7E11E0"),
    ),
    (
        256,
        sha3_256,
        sha3_256_sponge,
        hex2bytes("A7FFC6F8BF1ED76651C14756A061D662F580FF4DE43B49FA82D80A4B80F8434A"),
        hex2bytes("79F38ADEC5C20307A98EF76E8324AFBFD46CFD81B22E3973C65FA1BD9DE31787"),
    ),
    (
        384,
        sha3_384,
        sha3_384_sponge,
        hex2bytes("0C63A75B845E4F7D01107D852E4C2485C51A50AAAA94FC61995E71BBEE983A2AC3713831264ADB47FB6BD1E058D5F004"),
        hex2bytes("1881DE2CA7E41EF95DC4732B8F5F002B189CC1E42B74168ED1732649CE1DBCDD76197A31FD55EE989F2D7050DD473E8F"),
    ),
    (
        512,
        sha3_512,
        sha3_512_sponge,
        hex2bytes("A69F73CCA23A9AC5C8B567DC185A756E97C982164FE25859E0D1DCC1475C80A615B2123AF1F5F94C11E3E9402C3AC558F500199D95B6D3E301758586281DCD26"),
        hex2bytes("E76DFAD22084A8B1467FCF2FFA58361BEC7628EDF5F3FDC0E4805DC48CAEECA81B7C13C30ADF52A3659584739A2DF46BE589C51CA1A4A8416DF6545A1CE8BA00"),
    ),
]
    @test shafunc(UInt8[]) == shafunc(()) == Tuple(empty_ref)

    @test shafunc(fill(0xa3, 200)) == Tuple(len1600_ref)

    sp = spongefunc()
    for _ in 1:200
        sp = absorb(sp, tuple(0xa3))
    end
    @test squeeze(pad(sp), Val(d÷8))[2] == Tuple(len1600_ref)

    sp = spongefunc()
    for _ in 1:200
        sp = absorb(sp, [0xa3])
    end
    @test squeeze(pad(sp), Val(d÷8))[2] == Tuple(len1600_ref)

    # now verify with random input - where mis-indexing would matter - that chunked input
    # behaves identical to en-bloc input
    input = rand(UInt8, 10_000)
    sp = spongefunc()
    k = 0
    while k < length(input)
        l = min(rand(0:300), length(input)-k)
        chunk = input[k+1:k+l]
        if rand(Bool) # randomly choose between vector and tuple input
            sp = absorb(sp, chunk)
        else
            sp = absorb(sp, Tuple(chunk))
        end
        k += l
    end
    @test squeeze(pad(sp), Val(d÷8))[2] == shafunc(input)

    # now similarly for SIMD
    for N in 2:4
        len = 10_000
        input = [rand(UInt8, len) for _ in 1:N]
        sp = spongefunc(Val(N))
        k = 0
        while k < len
            l = min(rand(0:300), len-k)
            chunks = [inp[k+1:k+l] for inp in input]
            if l > 16 || rand(Bool) # randomly choose between vector and tuple input for short chunks
                sp = absorb(sp, chunks...)
            else
                sp = absorb(sp, Tuple.(chunks)...)
            end
            k += l
        end
        # compare chunk-wise SIMD, en-bloc SIMD, and en-bloc non-SIMD cases
        @test squeeze(pad(sp), Val(d÷8))[2] == shafunc(input...) == Tuple(shafunc(inp) for inp in input)
    end
end

@testset "SHAKE-$d" for (d, shakefunc, spongefunc, empty_ref, len1600_ref) in [
    # reference data from
    # https://csrc.nist.gov/projects/cryptographic-standards-and-guidelines/example-values
    # for input message lengths 0 bit and 1600 bit
    (
        128,
        shake_128,
        shake_128_sponge,
        hex2bytes("7F9C2BA4E88F827D616045507605853ED73B8093F6EFBC88EB1A6EACFA66EF263CB1EEA988004B93103CFB0AEEFD2A686E01FA4A58E8A3639CA8A1E3F9AE57E235B8CC873C23DC62B8D260169AFA2F75AB916A58D974918835D25E6A435085B2BADFD6DFAAC359A5EFBB7BCC4B59D538DF9A04302E10C8BC1CBF1A0B3A5120EA17CDA7CFAD765F5623474D368CCCA8AF0007CD9F5E4C849F167A580B14AABDEFAEE7EEF47CB0FCA9767BE1FDA69419DFB927E9DF07348B196691ABAEB580B32DEF58538B8D23F87732EA63B02B4FA0F4873360E2841928CD60DD4CEE8CC0D4C922A96188D032675C8AC850933C7AFF1533B94C834ADBB69C6115BAD4692D8619F90B0CDF8A7B9C264029AC185B70B83F2801F2F4B3F70C593EA3AEEB613A7F1B1DE33FD75081F592305F2E4526EDC09631B10958F464D889F31BA010250FDA7F1368EC2967FC84EF2AE9AFF268E0B1700AFFC6820B523A3D917135F2DFF2EE06BFE72B3124721D4A26C04E53A75E30E73A7A9C4A95D91C55D495E9F51DD0B5E9D83C6D5E8CE803AA62B8D654DB53D09B8DCFF273CDFEB573FAD8BCD45578BEC2E770D01EFDE86E721A3F7C6CCE275DABE6E2143F1AF18DA7EFDDC4C7B70B5E345DB93CC936BEA323491CCB38A388F546A9FF00DD4E1300B9B2153D2041D205B443E41B45A653F2A5C4492C1ADD544512DDA2529833462B71A41A45BE97290B6F"),
        hex2bytes("131AB8D2B594946B9C81333F9BB6E0CE75C3B93104FA3469D3917457385DA037CF232EF7164A6D1EB448C8908186AD852D3F85A5CF28DA1AB6FE3438171978467F1C05D58C7EF38C284C41F6C2221A76F12AB1C04082660250802294FB87180213FDEF5B0ECB7DF50CA1F8555BE14D32E10F6EDCDE892C09424B29F597AFC270C904556BFCB47A7D40778D390923642B3CBD0579E60908D5A000C1D08B98EF933F806445BF87F8B009BA9E94F7266122ED7AC24E5E266C42A82FA1BBEFB7B8DB0066E16A85E0493F07DF4809AEC084A593748AC3DDE5A6D7AAE1E8B6E5352B2D71EFBB47D4CAEED5E6D633805D2D323E6FD81B4684B93A2677D45E7421C2C6AEA259B855A698FD7D13477A1FE53E5A4A6197DBEC5CE95F505B520BCD9570C4A8265A7E01F89C0C002C59BFEC6CD4A5C109258953EE5EE70CD577EE217AF21FA70178F0946C9BF6CA8751793479F6B537737E40B6ED28511D8A2D7E73EB75F8DAAC912FF906E0AB955B083BAC45A8E5E9B744C8506F37E9B4E749A184B30F43EB188D855F1B70D71FF3E50C537AC1B0F8974F0FE1A6AD295BA42F6AEC74D123A7ABEDDE6E2C0711CAB36BE5ACB1A5A11A4B1DB08BA6982EFCCD716929A7741CFC63AA4435E0B69A9063E880795C3DC5EF3272E11C497A91ACF699FEFEE206227A44C9FB359FD56AC0A9A75A743CFF6862F17D7259AB075216C0699511643B6439"),
    ),
    (
        256,
        shake_256,
        shake_256_sponge,
        hex2bytes("46B9DD2B0BA88D13233B3FEB743EEB243FCD52EA62B81B82B50C27646ED5762FD75DC4DDD8C0F200CB05019D67B592F6FC821C49479AB48640292EACB3B7C4BE141E96616FB13957692CC7EDD0B45AE3DC07223C8E92937BEF84BC0EAB862853349EC75546F58FB7C2775C38462C5010D846C185C15111E595522A6BCD16CF86F3D122109E3B1FDD943B6AEC468A2D621A7C06C6A957C62B54DAFC3BE87567D677231395F6147293B68CEAB7A9E0C58D864E8EFDE4E1B9A46CBE854713672F5CAAAE314ED9083DAB4B099F8E300F01B8650F1F4B1D8FCF3F3CB53FB8E9EB2EA203BDC970F50AE55428A91F7F53AC266B28419C3778A15FD248D339EDE785FB7F5A1AAA96D313EACC890936C173CDCD0FAB882C45755FEB3AED96D477FF96390BF9A66D1368B208E21F7C10D04A3DBD4E360633E5DB4B602601C14CEA737DB3DCF722632CC77851CBDDE2AAF0A33A07B373445DF490CC8FC1E4160FF118378F11F0477DE055A81A9EDA57A4A2CFB0C83929D310912F729EC6CFA36C6AC6A75837143045D791CC85EFF5B21932F23861BCF23A52B5DA67EAF7BAAE0F5FB1369DB78F3AC45F8C4AC5671D85735CDDDB09D2B1E34A1FC066FF4A162CB263D6541274AE2FCC865F618ABE27C124CD8B074CCD516301B91875824D09958F341EF274BDAB0BAE316339894304E35877B0C28A9B1FD166C796B9CC258A064A8F57E27F2A"),
        hex2bytes("CD8A920ED141AA0407A22D59288652E9D9F1A7EE0C1E7C1CA699424DA84A904D2D700CAAE7396ECE96604440577DA4F3AA22AEB8857F961C4CD8E06F0AE6610B1048A7F64E1074CD629E85AD7566048EFC4FB500B486A3309A8F26724C0ED628001A1099422468DE726F1061D99EB9E93604D5AA7467D4B1BD6484582A384317D7F47D750B8F5499512BB85A226C4243556E696F6BD072C5AA2D9B69730244B56853D16970AD817E213E470618178001C9FB56C54FEFA5FEE67D2DA524BB3B0B61EF0E9114A92CDBB6CCCB98615CFE76E3510DD88D1CC28FF99287512F24BFAFA1A76877B6F37198E3A641C68A7C42D45FA7ACC10DAE5F3CEFB7B735F12D4E589F7A456E78C0F5E4C4471FFFA5E4FA0514AE974D8C2648513B5DB494CEA847156D277AD0E141C24C7839064CD08851BC2E7CA109FD4E251C35BB0A04FB05B364FF8C4D8B59BC303E25328C09A882E952518E1A8AE0FF265D61C465896973D7490499DC639FB8502B39456791B1B6EC5BCC5D9AC36A6DF622A070D43FED781F5F149F7B62675E7D1A4D6DEC48C1C7164586EAE06A51208C0B791244D307726505C3AD4B26B6822377257AA152037560A739714A3CA79BD605547C9B78DD1F596F2D4F1791BC689A0E9B799A37339C04275733740143EF5D2B58B96A363D4E08076A1A9D7846436E4DCA5728B6F760EEF0CA92BF0BE5615E96959D767197A0BEEB"),
    ),
]
    @test shakefunc(UInt8[], length(empty_ref)) == shakefunc((), length(empty_ref)) == empty_ref
    @test shakefunc(UInt8[], Val(length(empty_ref))) == shakefunc((), Val(length(empty_ref))) == Tuple(empty_ref)
    @test squeeze(shakefunc(UInt8[]), length(empty_ref))[2] == squeeze(shakefunc(()), length(empty_ref))[2] == empty_ref
    @test squeeze(shakefunc(UInt8[]), Val(length(empty_ref)))[2] == squeeze(shakefunc(()), Val(length(empty_ref)))[2] == Tuple(empty_ref)

    @test shakefunc(fill(0xa3, 200), length(len1600_ref)) == len1600_ref
    @test shakefunc(fill(0xa3, 200), Val(length(len1600_ref))) == Tuple(len1600_ref)
    @test squeeze(shakefunc(fill(0xa3, 200)), length(len1600_ref))[2] == len1600_ref
    @test squeeze(shakefunc(fill(0xa3, 200)), Val(length(len1600_ref)))[2] == Tuple(len1600_ref)

    sp = spongefunc()
    for _ in 1:200
        sp = absorb(sp, tuple(0xa3))
    end
    @test squeeze(pad(sp), length(len1600_ref))[2] == len1600_ref

    sp = spongefunc()
    for _ in 1:200
        sp = absorb(sp, [0xa3])
    end
    @test squeeze(pad(sp), length(len1600_ref))[2] == len1600_ref

    # now verify with random input - where mis-indexing would matter - that chunked
    # input/output behaves identical to en-bloc input/output
    input = rand(UInt8, 10_000)
    sp = spongefunc()
    k = 0
    while k < length(input)
        l = min(rand(0:300), length(input)-k)
        chunk = input[k+1:k+l]
        if l > 32 || rand(Bool) # randomly choose between vector and tuple input for short chunks
            sp = absorb(sp, chunk)
        else
            sp = absorb(sp, Tuple(chunk))
        end
        k += l
    end
    sp = pad(sp)
    output = Vector{UInt8}(undef, 5_000)
    k = 0
    while k < length(output)
        l = min(rand(0:300), length(output)-k)
        if l > 32 || rand(Bool) # randomly choose between vector and tuple output for short chunks
            sp, chunk = squeeze(sp, l)
        else
            sp, chunk = squeeze(sp, Val(l))
        end
        output[k+1:k+l] .= chunk
        k += l
    end
    @test output == shakefunc(input, length(output))

    # now similarly for SIMD
    for N in 2:4
        len = 10_000
        input = [rand(UInt8, len) for _ in 1:N]
        sp = spongefunc(Val(N))
        k = 0
        while k < len
            l = min(rand(0:300), len-k)
            chunks = [inp[k+1:k+l] for inp in input]
            if l > 16 || rand(Bool) # randomly choose between vector and tuple input for short chunks
                sp = absorb(sp, chunks...)
            else
                sp = absorb(sp, Tuple.(chunks)...)
            end
            k += l
        end
        sp = pad(sp)
        len = 5_000
        output = [Vector{UInt8}(undef, len) for _ in 1:N]
        k = 0
        while k < len
            l = min(rand(0:300), len-k)
            if l > 16 || rand(Bool) # randomly choose between vector and tuple output for short chunks
                sp, chunk = squeeze(sp, l)
            else
                sp, chunk = squeeze(sp, Val(l))
            end
            for n in 1:N
                output[n][k+1:k+l] .= chunk[n]
            end
            k += l
        end
        # compare chunk-wise SIMD, en-bloc SIMD, and en-bloc non-SIMD cases
        @test output == collect(shakefunc(input..., len)) == [shakefunc(inp, len) for inp in input]
    end
end
