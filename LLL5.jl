# LLL-Inplace.jl: Optimized LLL algorithm using direct MPFR/GMP ccalls.
# This code bypasses Julia's standard arithmetic to minimize GC pressure.

using Base.MPFR: libmpfr, MPFRRoundNearest, MPFRRoundingMode, CdoubleMax
using Base.GMP: libgmp, CulongMax

const bighalf = BigFloat(1)/2
const mbighalf = -BigFloat(1)/2
const Clong_0 = Clong(0)
const Culong_0 = Culong(0)
const Culong_1 = Culong(1)
const mpfrRN = MPFRRoundNearest

@inline function mpfr_set(rop::BigFloat , op::BigFloat, rnd::MPFRRoundingMode)::Int32
    ccall((:mpfr_set, libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, MPFRRoundingMode), rop, op, rnd)
end
@inline function mpfr_set_z(rop::BigFloat , op::BigInt, rnd::MPFRRoundingMode)::Int32
    ccall((:mpfr_set_z, libmpfr), Int32, (Ref{BigFloat}, Ref{BigInt}, MPFRRoundingMode), rop, op, rnd)
end
@inline function mpfr_set_ui(rop::BigFloat, op::CulongMax, rnd::MPFRRoundingMode)::Int32
    ccall((:mpfr_set_ui, libmpfr), Int32, (Ref{BigFloat}, Culong, MPFRRoundingMode), rop, op, rnd)
end
@inline function mpfr_neg(rop::BigFloat, op::BigFloat, rnd::MPFRRoundingMode)::Int32
    ccall((:mpfr_neg, libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, MPFRRoundingMode), rop, op, rnd)
end
@inline function mpfr_cmpabs_ui(op1::BigFloat, op2::CulongMax)::Int32
    return ccall((:mpfr_cmpabs_ui, libmpfr), Int32, (Ref{BigFloat}, Culong), op1, op2)
end
@inline function mpfr_greaterequal_p(op1::BigFloat, op2::BigFloat)::Int32
    return ccall((:mpfr_greaterequal_p, libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}), op1, op2)
end
@inline function mpfr_sgn(op::BigFloat)::Int32
    return ccall((:mpfr_sgn, libmpfr), Int32, (Ref{BigFloat},), op)
end
@inline function mpfr_fits_ulong_p(op::BigFloat,rnd::MPFRRoundingMode)::Int32
    return ccall((:mpfr_fits_ulong_p, libmpfr), Int32, (Ref{BigFloat}, MPFRRoundingMode), op, rnd)
end
@inline function mpfr_get_ui(op::BigFloat, rnd::MPFRRoundingMode)::UInt32
    return ccall((:mpfr_get_ui, libmpfr), UInt32, (Ref{BigFloat}, MPFRRoundingMode), op, rnd)
end
@inline function mpfr_get_z(rop::BigInt, op::BigFloat, rnd::MPFRRoundingMode)::Int32
    ccall((:mpfr_get_z, libmpfr), Int32, (Ref{BigInt}, Ref{BigFloat}, MPFRRoundingMode), rop, op, rnd)
end
@inline function mpfr_mul(rop::BigFloat, op1::BigFloat, op2::BigFloat, rnd::MPFRRoundingMode)::Int32
    ccall((:mpfr_mul, libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, MPFRRoundingMode), rop, op1, op2, rnd)
end
@inline function mpfr_mul_ui(rop::BigFloat, op1::BigFloat, op2::CulongMax, rnd::MPFRRoundingMode)::Int32
    ccall((:mpfr_mul_ui, libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Culong, MPFRRoundingMode), rop, op1, op2, rnd)
end
@inline function mpfr_mul_z(rop::BigFloat, op1::BigFloat, op2::BigInt, rnd::MPFRRoundingMode)::Int32
    ccall((:mpfr_mul_z, libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigInt}, MPFRRoundingMode), rop, op1, op2, rnd)
end
@inline function mpfr_mul_2ui(rop::BigFloat, op1::BigFloat, op2::CulongMax, rnd::MPFRRoundingMode)::Int32
    ccall((:mpfr_mul_2ui, libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Culong, MPFRRoundingMode), rop, op1, op2, rnd)
end
@inline function mpfr_fma(rop::BigFloat, op1::BigFloat, op2::BigFloat, op3::BigFloat, rnd::MPFRRoundingMode)::Int32
    ccall((:mpfr_fma, libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, MPFRRoundingMode), rop, op1, op2, op3, rnd)
end
@inline function mpfr_fms(rop::BigFloat, op1::BigFloat, op2::BigFloat, op3::BigFloat, rnd::MPFRRoundingMode)::Int32
    ccall((:mpfr_fms, libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, MPFRRoundingMode), rop, op1, op2, op3, rnd)
end
@inline function mpfr_sqr(rop::BigFloat, op::BigFloat, rnd::MPFRRoundingMode)::Int32
    ccall((:mpfr_sqr, libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, MPFRRoundingMode), rop, op, rnd)
end
@inline function mpfr_add(rop::BigFloat, op1::BigFloat, op2::BigFloat, rnd::MPFRRoundingMode)::Int32
    ccall((:mpfr_add, libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, MPFRRoundingMode), rop, op1, op2, rnd)
end
@inline function mpfr_sub(rop::BigFloat, op1::BigFloat, op2::BigFloat, rnd::MPFRRoundingMode)::Int32
    ccall((:mpfr_sub, libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, MPFRRoundingMode), rop, op1, op2, rnd)
end
@inline function mpfr_div(rop::BigFloat, op1::BigFloat, op2::BigFloat, rnd::MPFRRoundingMode)::Int32
    ccall((:mpfr_div, libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, MPFRRoundingMode), rop, op1, op2, rnd)
end
@inline function mpz_set_ui(rop::BigInt, op::CulongMax)
    ccall((:__gmpz_set_ui, libgmp), Cvoid, (Ref{BigInt}, Culong), rop, op)
end
@inline function mpz_swap(rop1::BigInt, rop2::BigInt)
    ccall((:__gmpz_swap, libgmp), Cvoid, (Ref{BigInt}, Ref{BigInt}), rop1, rop2)
end
@inline function mpz_addmul_ui!(rop::BigInt, op1::BigInt, op2::CulongMax)
    ccall((:__gmpz_addmul_ui, libgmp), Cvoid, (Ref{BigInt}, Ref{BigInt}, Culong), rop, op1, op2)
end
@inline function mpz_addmul!(rop::BigInt, op1::BigInt, op2::BigInt)
    ccall((:__gmpz_addmul, libgmp), Cvoid, (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), rop, op1, op2)
end
@inline function mpz_submul_ui!(rop::BigInt, op1::BigInt, op2::CulongMax)
    ccall((:__gmpz_submul_ui, libgmp), Cvoid, (Ref{BigInt}, Ref{BigInt}, Culong), rop, op1, op2)
end
@inline function mpz_submul!(rop::BigInt, op1::BigInt, op2::BigInt)
    ccall((:__gmpz_submul, libgmp), Cvoid, (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), rop, op1, op2)
end

function dot!(rop::BigFloat, bi::AbstractVector{BigFloat}, bj::AbstractVector{BigFloat})
    mpfr_set_ui(rop, Culong_0, mpfrRN)
    @inbounds for k in eachindex(bi,bj)
        mpfr_fma(rop, bi[k], bj[k], rop, mpfrRN)
    end
end
function norm2!(rop::BigFloat, b::AbstractVector{BigFloat})
    mpfr_set_ui(rop, Culong_0, mpfrRN)
    @inbounds for k in eachindex(b)
        mpfr_fma(rop, b[k], b[k], rop, mpfrRN)
    end
end

@inline function GSO3!(mu::Matrix{BigFloat},Bs::Vector{BigFloat}, b_int::Matrix{BigInt}, b_real::Vector{BigFloat})

    n = size(b_int, 1)
    b =[BigFloat(0) for i in 1:n, j in 1:n+1]     # GSOのため一度だけMatrix{BigFloat}に戻す

    BFbuf = BigFloat(0)

    @inbounds for i in 1:n, j in 1:n
        mpfr_set_z(b[i,j], b_int[i,j], mpfrRN)
    end
    @inbounds for i in 1:n
        mpfr_set(b[i,n+1], b_real[i], mpfrRN)
    end

    bs = [BigFloat(0) for i in 1:n, j in 1:n+1]

    bsi = [BigFloat(0) for _ in 1:n+1]

    @views begin
        for i in 1:n+1
            mpfr_set(bs[1,i], b[1,i], mpfrRN)
        end
        norm2!(Bs[1], bs[1,:])
    end
    
    mu[1,1] = BigFloat(1)

    @views for i in 2:n
        for j in 1:n+1
            mpfr_set(bsi[j], b[i,j], mpfrRN)
        end
        
        for j in 1:i-1

            dot!(BFbuf, b[i,:], bs[j,:])
            mpfr_div(mu[i,j], BFbuf, Bs[j], mpfrRN)

            for l in 1:n+1
                mpfr_neg(BFbuf, mu[i,j], mpfrRN)
                mpfr_fma(bsi[l], BFbuf, bs[j,l], bsi[l], mpfrRN)
            end

        end
        for j in 1:n+1
            mpfr_set(bs[i,j], bsi[j], mpfrRN)
        end
        norm2!(Bs[i],bsi)
        mpfr_set_ui(mu[i,i], Culong_1, mpfrRN)
    end
end

function LLL5(x::Vector{BigFloat},C::BigInt,delta::BigFloat=BigFloat(3)/4,maxiter=10^6)
    
    #@assert delta > BigFloat(1)/4 && delta < one(BigFloat) "LLL requires 1/4 < delta < 1"

    n = lastindex(x)
    k = 2
    iter = 1
    itlim = 1
    q_Big = BigInt()
    mukj_abs = BigFloat()
    buf_BF = BigFloat()
    mukj_double = BigFloat()
    t = BigFloat()
    mu_ = BigFloat()
    B = BigFloat()

    b_int = [BigInt() for i in 1:n, j in 1:n]
    b_real = [BigFloat() for i in 1:n]  # b[i,n+1]

    for i in 1:n
        mpz_set_ui(b_int[i,i],Culong_1)
        mpfr_mul_z(b_real[i], x[i], C, mpfrRN)
    end

    mu = [BigFloat(0) for i in 1:n, j in 1:n]
    Bs = [BigFloat(0) for i in 1:n]

    GSO3!(mu, Bs, b_int, b_real)

    while k <= n && iter < maxiter
        iter += 1
        itlim = max(itlim,k)

        for j in k-1:-1:1
            
            mukj = mu[k,j]

            mpfr_mul_2ui(mukj_double,mukj,Culong_1,mpfrRN)

            if mpfr_cmpabs_ui(mukj_double,Culong_1) > Clong_0
                # |mukj| > 0.5

                if mpfr_sgn(mukj) > Clong_0
                    # mukj is positive

                    if Clong_0 != mpfr_fits_ulong_p(mukj,mpfrRN)
                        # if mukj fits Culong
                        # q is Culong

                        q = mpfr_get_ui(mukj,mpfrRN)
                        
                        @views begin
                            rowk = b_int[k, :]
                            rowj = b_int[j, :]
                            @inbounds for l in 1:itlim
                                mpz_submul_ui!(rowk[l],rowj[l],q)   # rowk[l] = rowk[l] - rowj[l] * q
                            end
                        end

                        mpfr_mul_ui(buf_BF, b_real[j], q, mpfrRN)     # buf_BF = b_real[j] * q
                        mpfr_sub(b_real[k], b_real[k], buf_BF, mpfrRN)  # b_real[k] = b_real[k] - buf_BF
                        
                        @inbounds for l in 1:j
                            mpfr_mul_ui(buf_BF, mu[j,l], q, mpfrRN)
                            mpfr_sub(mu[k,l], mu[k,l], buf_BF, mpfrRN)
                        end
                    
                    else
                        # if mukj doesn't fit Culong
                        # q is BigInt

                        mpfr_get_z(q_Big,mukj,mpfrRN)

                        @views begin
                            rowk = b_int[k, :]
                            rowj = b_int[j, :]
                            @inbounds for l in 1:itlim
                                mpz_submul!(rowk[l], rowj[l], q_Big)
                            end
                        end

                        mpfr_mul_z(buf_BF, b_real[j], q_Big, mpfrRN)
                        mpfr_sub(b_real[k], b_real[k], buf_BF, mpfrRN)

                        @inbounds for l in 1:j
                            mpfr_mul_z(buf_BF, mu[j,l], q_Big, mpfrRN)
                            mpfr_sub(mu[k,l], mu[k,l], buf_BF, mpfrRN)
                        end
                    end

                else
                    # mukj < -0.5
                    # mukj is negative
                    mpfr_neg(mukj_abs,mukj,mpfrRN)

                    if Clong_0 != mpfr_fits_ulong_p(mukj_abs,mpfrRN)
                        # if mukj fits Culong
                        # q is Culong

                        q = mpfr_get_ui(mukj_abs,mpfrRN)
                        
                        @views begin
                            rowk = b_int[k, :]
                            rowj = b_int[j, :]
                            @inbounds for l in 1:itlim
                                mpz_addmul_ui!(rowk[l],rowj[l],q)   # rowk[l] = rowk[l] + rowj[l] * q
                            end
                        end

                        mpfr_mul_ui(buf_BF, b_real[j], q, mpfrRN)     # buf_BF = b_real[j] * q
                        mpfr_add(b_real[k], b_real[k], buf_BF, mpfrRN)  # b_real[k] = b_real[k] + buf_BF
                        
                        @inbounds for l in 1:j
                            mpfr_mul_ui(buf_BF, mu[j,l], q, mpfrRN)
                            mpfr_add(mu[k,l], mu[k,l], buf_BF, mpfrRN)
                        end
                    else
                        # if mukj doesn't fit Culong
                        # q is BigInt

                        mpfr_get_z(q_Big,mukj_abs,mpfrRN)

                        @views begin
                            rowk = b_int[k, :]
                            rowj = b_int[j, :]
                            @inbounds for l in 1:itlim
                                mpz_addmul!(rowk[l], rowj[l], q_Big)
                            end
                        end

                        mpfr_mul_z(buf_BF, b_real[j], q_Big, mpfrRN)
                        mpfr_add(b_real[k], b_real[k], buf_BF, mpfrRN)

                        @inbounds for l in 1:j
                            mpfr_mul_z(buf_BF, mu[j,l], q_Big, mpfrRN)
                            mpfr_add(mu[k,l], mu[k,l], buf_BF, mpfrRN)
                        end
                    end

                end
            end
        end
        
        mpfr_sqr(buf_BF, mu[k,k-1], mpfrRN)         # buf_BF = mu[k,k-1]^2
        mpfr_sub(buf_BF, delta, buf_BF, mpfrRN)     # buf_BF = delta - mu[k,k-1]^2
        mpfr_mul(buf_BF, buf_BF, Bs[k-1], mpfrRN)   # buf_BF = (delta - mu[k,k-1]^2)*Bs[k-1]

        if mpfr_greaterequal_p(Bs[k], buf_BF) != Clong_0

            k += 1
        
        else
        
            for i in 1:n
                mpz_swap(b_int[k-1,i],b_int[k,i])
            end
            b_real[k-1]  , b_real[k]  = b_real[k]  , b_real[k-1]

            mpfr_set(mu_, mu[k,k-1], mpfrRN)    # mu_ = mu[k,k-1]
            mpfr_sqr(buf_BF, mu_, mpfrRN)
            mpfr_fma(B, buf_BF, Bs[k-1], Bs[k], mpfrRN)     # B = Bs[k] + mu_* mu_ * Bs[k-1]
            
            mpfr_div(buf_BF, Bs[k-1], B, mpfrRN)
            mpfr_mul(mu[k,k-1], mu_, buf_BF, mpfrRN)    # mu[k,k-1] = mu_ * Bs[k-1] / B
            mpfr_mul(Bs[k], Bs[k], buf_BF, mpfrRN)      # Bs[k] = Bs[k] * Bs[k-1] / B
            mpfr_set(Bs[k-1], B, mpfrRN)                # Bs[k-1] = B
            
            @inbounds for j in 1:k-2
                mu[k-1,j] , mu[k,j] = mu[k,j] , mu[k-1,j]   # swap
            end
            
            @inbounds for j in k+1:n
                mpfr_set(t,mu[j,k],mpfrRN)

                mpfr_fms(buf_BF, mu_, t, mu[j,k-1], mpfrRN)
                mpfr_neg(mu[j,k], buf_BF, mpfrRN)
                mpfr_fma(mu[j,k-1], mu[k,k-1], mu[j,k], t, mpfrRN)
            end
            
            k = max(2,k-1)

        end
    end
    return b_int, b_real, iter
end

