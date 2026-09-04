! The operations bench_logical times, in their own compilation unit.
!
! This is the whole point of the separate file: with everything in one program, ifx sees
! that a copied mask is only ever counted, rewrites count(dst) as count(src), and reports a
! copy that costs a tenth of what the same copy costs when it cannot be elided. Compile
! without -flto / -ipo and the calls below are opaque, which is also what a real generated
! wrapper is to its caller.
module bench_kernels
    use, intrinsic :: iso_fortran_env, only: real64, int32, int64
    use, intrinsic :: iso_c_binding, only: c_bool
    implicit none
    private
    public :: marshal_in, marshal_out, marshal_automatic, make_default, make_cbool
    public :: read_mask, read_mask_cbool, branch_default, branch_cbool

contains

    !> intent(in) marshalling: the caller's c_bool buffer into a default logical.
    subroutine marshal_in(n, src, dst, sink)
        integer(int32), intent(in) :: n
        logical(c_bool), intent(in) :: src(n)
        logical, intent(out) :: dst(n)
        integer(int64), intent(inout) :: sink

        dst = src
        sink = sink + count(dst)
    end subroutine marshal_in

    !> intent(out) marshalling: the kernel's default logical back into c_bool.
    subroutine marshal_out(n, kern, out_c, sink)
        integer(int32), intent(in) :: n
        logical, intent(in) :: kern(n)
        logical(c_bool), intent(out) :: out_c(n)
        integer(int64), intent(inout) :: sink

        out_c = kern
        sink = sink + count(out_c)
    end subroutine marshal_out

    !> What a generated wrapper actually declares: an automatic array, sized at run time
    !| and allocated per call.
    subroutine marshal_automatic(n, src, sink)
        integer(int32), intent(in) :: n
        logical(c_bool), intent(in) :: src(n)
        integer(int64), intent(inout) :: sink
        logical :: dst(n)

        dst = src
        sink = sink + count(dst)
    end subroutine marshal_automatic

    !> Computing a mask into a default logical, and into c_bool: the second-order cost of
    !| moving the implementations to c_bool, where every comparison then needs converting.
    subroutine make_default(n, x, threshold, mask, sink)
        integer(int32), intent(in) :: n
        real(real64), intent(in) :: x(n), threshold
        logical, intent(out) :: mask(n)
        integer(int64), intent(inout) :: sink

        mask = x > threshold
        sink = sink + count(mask)
    end subroutine make_default

    subroutine make_cbool(n, x, threshold, mask, sink)
        integer(int32), intent(in) :: n
        real(real64), intent(in) :: x(n), threshold
        logical(c_bool), intent(out) :: mask(n)
        integer(int64), intent(inout) :: sink

        mask = x > threshold
        sink = sink + count(mask)
    end subroutine make_cbool

    !> Reading a mask with `count`, in both kinds. This is the pair that matters now that
    !| every mask in the framework is `logical(c_bool)`: the marshalling columns measure a
    !| copy that no longer exists, and `mk:` only covers writing one. A whole-array
    !| reduction is the vectorisable read, where a quarter of the bytes should show.
    subroutine read_mask(n, mask, sink)
        integer(int32), intent(in) :: n
        logical, intent(in) :: mask(n)
        integer(int64), intent(inout) :: sink

        sink = sink + count(mask)
    end subroutine read_mask

    subroutine read_mask_cbool(n, mask, sink)
        integer(int32), intent(in) :: n
        logical(c_bool), intent(in) :: mask(n)
        integer(int64), intent(inout) :: sink

        sink = sink + count(mask)
    end subroutine read_mask_cbool

    !> The read the implementations actually perform: a scalar element test per iteration,
    !| guarding work on the selected elements (`if (vectors_selection_mask(i_vec))`,
    !| `if (.not. neighbor_mask(i_neighbor, i_point)) cycle`). Branch behaviour, not
    !| bandwidth, may dominate here -- which is the point of measuring it separately.
    subroutine branch_default(n, mask, x, sink)
        integer(int32), intent(in) :: n
        logical, intent(in) :: mask(n)
        real(real64), intent(in) :: x(n)
        integer(int64), intent(inout) :: sink
        integer(int32) :: i
        real(real64) :: acc

        acc = 0.0_real64
        do i = 1, n
            if (mask(i)) acc = acc + x(i)
        end do
        sink = sink + int(acc, int64)
    end subroutine branch_default

    subroutine branch_cbool(n, mask, x, sink)
        integer(int32), intent(in) :: n
        logical(c_bool), intent(in) :: mask(n)
        real(real64), intent(in) :: x(n)
        integer(int64), intent(inout) :: sink
        integer(int32) :: i
        real(real64) :: acc

        acc = 0.0_real64
        do i = 1, n
            if (mask(i)) acc = acc + x(i)
        end do
        sink = sink + int(acc, int64)
    end subroutine branch_cbool

end module bench_kernels
