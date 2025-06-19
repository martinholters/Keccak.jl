using Keccak
using OffsetArrays: Origin
using SIMD: Vec
using Test: @test, @testset

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
