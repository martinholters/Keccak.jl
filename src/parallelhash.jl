const ParallelHash_N = Tuple(codeunits("ParallelHash"))

for d in [128, 256]
    R = (1600-2*d) ÷ 64
    shakefunc = Symbol("shake_$(d)")
    cshake_spongefunc = Symbol("cshake_$(d)_sponge")
    parallelhashhelper = Symbol("_parallelhash_$(d)")
    parallelhashfunc = Symbol("parallelhash_$(d)")
    parallelhashxoffunc = Symbol("parallelhash_xof_$(d)")
    @eval begin
        function $(parallelhashhelper)(
                data,
                blocksize,
                ::Val{len},
                S::AbsorbableData,
                threaded::Val,
            ) where {len}
            return $(parallelhashhelper)(data, blocksize, len, S, threaded)
        end
        function $(parallelhashhelper)(
                data::String,
                blocksize,
                len::Integer,
                S::AbsorbableData,
                threaded::Val,
            )
            return $(parallelhashhelper)(codeunits(data), blocksize, len, S, threaded)
        end
        function $(parallelhashhelper)(
                data::AbstractVector{UInt8},
                blocksize::Integer,
                len::Integer,
                S::AbsorbableData,
                ::Val{false},
            )
            sponge = $(cshake_spongefunc)(ParallelHash_N, S)
            sponge = absorb_left_encoded(sponge, blocksize)
            n = div(absorblength(data), blocksize, RoundUp)
            offset = 0
            while offset+4blocksize <= length(data)
                zi = $(shakefunc)(
                    @view(data[1+offset:offset+blocksize]),
                    @view(data[1+offset+blocksize:offset+2blocksize]),
                    @view(data[1+offset+2blocksize:offset+3blocksize]),
                    @view(data[1+offset+3blocksize:offset+4blocksize]),
                    Val($(d÷4)),
                )
                sponge = absorb(sponge, (zi[1]..., zi[2]..., zi[3]..., zi[4]...))
                offset += 4*blocksize
            end
            if offset+2blocksize <= length(data)
                zi = $(shakefunc)(
                    @view(data[1+offset:offset+blocksize]),
                    @view(data[1+offset+blocksize:offset+2blocksize]),
                    Val($(d÷4)),
                )
                sponge = absorb(sponge, (zi[1]..., zi[2]...))
                offset += 2*blocksize
            end
            if offset+blocksize <= length(data)
                zi = $(shakefunc)(
                    @view(data[1+offset:offset+blocksize]),
                    Val($(d÷4)),
                )
                sponge = absorb(sponge, zi)
                offset += blocksize
            end
            if offset < length(data)
                zi = $(shakefunc)(
                    @view(data[1+offset:end]),
                    Val($(d÷4)),
                )
                sponge = absorb(sponge, zi)
            end
            sponge = absorb_right_encoded(sponge, n)
            sponge = absorb_right_encoded(sponge, 8*len)
            sponge = pad(sponge)
            return sponge
        end

        function $(parallelhashhelper)(
                data::AbstractVector{UInt8},
                blocksize::Integer,
                len::Integer,
                S::AbsorbableData,
                ::Val{true},
            )
            sponge = $(cshake_spongefunc)(ParallelHash_N, S)
            sponge = absorb_left_encoded(sponge, blocksize)
            n = div(absorblength(data), blocksize, RoundUp)
            n_threaded = (length(data) ÷ (4 * blocksize)) * 4
            zi_threaded = Vector{NTuple{$(d÷4),UInt8}}(undef, n_threaded)
            Threads.@threads for i in 0:4:n_threaded-1
                zi_threaded[i+1:i+4] .= $(shakefunc)(
                    @view(data[1+i*blocksize:(i+1)*blocksize]),
                    @view(data[1+(i+1)*blocksize:(i+2)*blocksize]),
                    @view(data[1+(i+2)*blocksize:(i+3)*blocksize]),
                    @view(data[1+(i+3)*blocksize:(i+4)*blocksize]),
                    Val($(d÷4)),
                )
            end
            sponge = absorb(sponge, reinterpret(UInt8, zi_threaded))
            offset = n_threaded*blocksize
            if offset+2blocksize <= length(data)
                zi = $(shakefunc)(
                    @view(data[1+offset:offset+blocksize]),
                    @view(data[1+offset+blocksize:offset+2blocksize]),
                    Val($(d÷4)),
                )
                sponge = absorb(sponge, (zi[1]..., zi[2]...))
                offset += 2*blocksize
            end
            if offset+blocksize <= length(data)
                zi = $(shakefunc)(
                    @view(data[1+offset:offset+blocksize]),
                    Val($(d÷4)),
                )
                sponge = absorb(sponge, zi)
                offset += blocksize
            end
            if offset < length(data)
                zi = $(shakefunc)(
                    @view(data[1+offset:end]),
                    Val($(d÷4)),
                )
                sponge = absorb(sponge, zi)
            end
            sponge = absorb_right_encoded(sponge, n)
            sponge = absorb_right_encoded(sponge, 8*len)
            sponge = pad(sponge)
            return sponge
        end

        """
            $($(parallelhashfunc))(data, blocksize, len, S=(); threaded=true)

        Computes the ParallelHash$($(d)) of length `len` with block size `blocksize` from
        `data` and the optional customization string `S`.

        The `data` can be a `AbstractVector` of `UInt8`s or a`String`. (Contrary to the
        other hashing functions, `Tuple`s are not supported). It is processed in chunks
        of size `blocksize`, using SIMD and (unless disabled by passing `threaded=false`)
        multi-threading to parallelize the computation. Note that changing `blocksize` not
        only affects performance, but also the result, so it is not purely an optimization.

        The hash is returned as a tuple of `UInt8`s if `len` is a `Val`, as a
        `Vector{UInt8}` otherwise.
        """
        function $(parallelhashfunc)(
                data::Union{AbstractVector{UInt8}, String},
                blocksize::Integer,
                len::Union{Val, Integer},
                S::AbsorbableData = ();
                threaded = true,
            )
            sponge = threaded ?
                $(parallelhashhelper)(data, blocksize, len, S, Val(true)) :
                $(parallelhashhelper)(data, blocksize, len, S, Val(false))
            return squeeze(sponge, len)[2]
        end

        """
            $($(parallelhashxoffunc))(data, blocksize, S=(); threaded=true)

        Creates a sponge for the arbitrary-length output ParallelHashXOF$($(d)) with
        block size `blocksize` from `data` and the optional customization string `S`.

        The `data` can be a `AbstractVector` of `UInt8`s or a`String`. (Contrary to the
        other hashing functions, `Tuple`s are not supported). It is processed in chunks
        of size `blocksize`, using SIMD and (unless disabled by passing `threaded=false`)
        multi-threading to parallelize the computation. Note that changing `blocksize` not
        only affects performance, but also the result, so it is not purely an optimization.

        The returned sponge is ready to squeeze output of the desired length from.
        """
        function $(parallelhashxoffunc)(
            data::AbsorbableData,
            blocksize::Integer,
            S::AbsorbableData = ();
            threaded = true,
        )
            if threaded
                return $(parallelhashhelper)(data, blocksize, 0, S, Val(true))
            else
                return $(parallelhashhelper)(data, blocksize, 0, S, Val(false))
            end
        end

        """
            $($(parallelhashxoffunc))(data, blocksize, len, S=(); threaded=true)

        Computes the arbitrary-length output ParallelHashXOF$($(d)) of length `len` with
        block size `blocksize` from `data` and the optional customization string `S`.

        The `data` can be a `AbstractVector` of `UInt8`s or a`String`. (Contrary to the
        other hashing functions, `Tuple`s are not supported). It is processed in chunks
        of size `blocksize`, using SIMD and (unless disabled by passing `threaded=false`)
        multi-threading to parallelize the computation. Note that changing `blocksize` not
        only affects performance, but also the result, so it is not purely an optimization.

        The hash is returned as a tuple of `UInt8`s if `len` is a `Val`, as a
        `Vector{UInt8}` otherwise.

        The difference to `$($(parallelhashfunc))` is that the `len` is not itself hashed.
        """
        function $(parallelhashxoffunc)(
            data::AbsorbableData,
            blocksize::Integer,
            len::Union{Val,Integer},
            S::AbsorbableData = ();
            threaded = true,
        )
            return squeeze($(parallelhashxoffunc)(data, blocksize, S; threaded), len)[2]
        end
    end
end
