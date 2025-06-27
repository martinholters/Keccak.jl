using Keccak
using OffsetArrays: Origin
using SIMD: Vec
using Test: @test, @test_throws, @testset

@testset "Sponge" begin
    @testset "identity sponge" for (maybe_tup, maybe_val) in ((identity, identity), (Tuple, Val))
        sponge = Keccak.Sponge{3}(identity, Tuple(zeros(UInt16,7)), 0)
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
        sponge = Keccak.Sponge{3}(f, Tuple(zeros(UInt16,7)), 0)
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
        ref_sponges = [Keccak.Sponge{3}(f, Tuple(zeros(UInt16,7)), 0) for _ in 1:N]
        simd_sponge = Keccak.Sponge{3}(f, Tuple(zeros(Vec{N,UInt16},7)), 0)
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
    # test data from
    # https://csrc.nist.gov/projects/cryptographic-standards-and-guidelines/example-values
    # for SHAKE-128, 0-bit input
    empty_ref = hex2bytes("7F9C2BA4E88F827D616045507605853ED73B8093F6EFBC88EB1A6EACFA66EF263CB1EEA988004B93103CFB0AEEFD2A686E01FA4A58E8A3639CA8A1E3F9AE57E235B8CC873C23DC62B8D260169AFA2F75AB916A58D974918835D25E6A435085B2BADFD6DFAAC359A5EFBB7BCC4B59D538DF9A04302E10C8BC1CBF1A0B3A5120EA17CDA7CFAD765F5623474D368CCCA8AF0007CD9F5E4C849F167A580B14AABDEFAEE7EEF47CB0FCA9767BE1FDA69419DFB927E9DF07348B196691ABAEB580B32DEF58538B8D23F87732EA63B02B4FA0F4873360E2841928CD60DD4CEE8CC0D4C922A96188D032675C8AC850933C7AFF1533B94C834ADBB69C6115BAD4692D8619F90B0CDF8A7B9C264029AC185B70B83F2801F2F4B3F70C593EA3AEEB613A7F1B1DE33FD75081F592305F2E4526EDC09631B10958F464D889F31BA010250FDA7F1368EC2967FC84EF2AE9AFF268E0B1700AFFC6820B523A3D917135F2DFF2EE06BFE72B3124721D4A26C04E53A75E30E73A7A9C4A95D91C55D495E9F51DD0B5E9D83C6D5E8CE803AA62B8D654DB53D09B8DCFF273CDFEB573FAD8BCD45578BEC2E770D01EFDE86E721A3F7C6CCE275DABE6E2143F1AF18DA7EFDDC4C7B70B5E345DB93CC936BEA323491CCB38A388F546A9FF00DD4E1300B9B2153D2041D205B443E41B45A653F2A5C4492C1ADD544512DDA2529833462B71A41A45BE97290B6F")
    sponge = Keccak.absorb(Keccak.KeccakSponge{21,UInt64}(), [0x1f; zeros(UInt8, 21*8-2); 0x80])
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
    sponge = Keccak.absorb(Keccak.KeccakSponge{21,UInt64}(), [fill(0xa3, 200); 0x1f; zeros(UInt8, 134); 0x80])
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
    sponge = Keccak.KeccakSponge{21,UInt64}()
    for _ in 1:200
        sponge = Keccak.absorb(sponge, tuple(0xa3))
    end
    sponge = Keccak.absorb(sponge, [0x1f; zeros(UInt8, 134); 0x80])
    @test Keccak.squeeze(sponge, length(len1600_ref))[2] == len1600_ref

    # piece-wise input
    sponge = Keccak.KeccakSponge{21,UInt64}()
    k = 0
    while k < 200
        len = min(rand(0:10), 200-k)
        sponge = Keccak.absorb(sponge, fill(0xa3, len))
        k += len
    end
    sponge = Keccak.absorb(sponge, [0x1f; zeros(UInt8, 134); 0x80])
    @test Keccak.squeeze(sponge, length(len1600_ref))[2] == len1600_ref
end
