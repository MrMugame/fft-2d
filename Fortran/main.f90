module img_mod
  implicit none

  type Image
     integer :: width, height
     character(len=1), allocatable :: pixels(:, :)
  end type Image
contains
  type(Image) function make_image(width, height)
    integer :: width, height

    make_image%width = width;
    make_image%height = height;

    allocate(make_image%pixels(width, height))
    make_image%pixels = achar(0)
  end function make_image

  subroutine save_pgm(path, img)
    character(*), intent(in) :: path
    type(Image), intent(in) :: img
    integer :: i, j
    
    open(1, file=path)
   
    write(1, "(A/I0AI0/A)") "P5", img%width, " ", img%height, "255"
    ! I don't know if there is a better way to sovle the newline issue but this atleast works
    do i = 1, img%height
       do j = 1, img%width
          write(1, "(A)", advance="no") img%pixels(j, i)
       end do
    end do
    
    close(1)
  end subroutine save_pgm
  
  type(Image) function load_pgm(path)
    character(*) :: path
    character(len = 256) :: buf
    integer :: width, height
    character :: char
    integer :: max, i, j, status

    open(1, file=path, status="OLD")
    
    read(1, "(A)") buf

    if (buf(:2) /= "P5") stop "File type is not valid"

    do 
       read(1, "(A)") buf
       if (buf(1:1) /= "#") exit
    end do

    read(buf, *) width, height
    
    read(1, *) max
    if (max /= 255) stop "Wrong max value in file"
    
    load_pgm = make_image(width, height);
    load_pgm%width = width;
    load_pgm%height = height;

    do i = 1, height
       do j = 1, width
          read(1, "(A1)", advance="no", iostat=status) char
          load_pgm%pixels(j, i) = char
       end do
    end do

    close(1)
  end function load_pgm
end module img_mod

module processing_mod
  real(kind=8) :: PI = 4 * atan(1.0)
contains
  ! TOOD: Don't know if passing n here is necessary or if we could use len() or smth
  subroutine fft_1D(buffer, n)
    complex(kind=8), intent(inout) :: buffer(n)
    integer :: m, bit, j, k
    complex(kind=8) :: temp, w_m, w, t, u

    if (n == 0 .or. (iand(n, n - 1) /= 0)) stop "Buffer is not a power of 2"
    
    j = 0
    do k = 1, n-1
       bit = n;
       
       do
          bit = ishft(bit, -1)
          
          j = ieor(j, bit)

          if (iand(j, bit) /= 0 .or. bit == 1) exit
       end do
       
       ! Add one, because this language thinks it's lua (Or lua thinks it's this language)
       if (k < j) then
          temp = buffer(j+1)
          buffer(j+1) = buffer(k+1)
          buffer(k+1) = temp
       end if
    end do

    m = 2
    do 
       if (m > n) exit
       
       w_m = cdexp(complex(0, -2.0 * PI / m))

       do k = 0, n-1, m
          w = 1
          do j = 1, m/2
             t = w * buffer(k + j + m/2)
             u = buffer(k + j)

             buffer(k + j) = u + t
             buffer(k + j + m/2) = u - t
             w = w * w_m
          end do
       end do

       m = ishft(m, 1)
    end do
  end subroutine fft_1D
  
  ! TODO: Again dont know if necessary
  subroutine fft_2D(buffer, width, height)
    integer, intent(in) :: width, height
    complex(kind=8), intent(inout) :: buffer(width, height)
    integer :: i

    do i = 1, height
       call fft_1D(buffer(i, 1:width), width)
    end do

    buffer = transpose(buffer)
    
    do i = 1, width
       call fft_1D(buffer(i, 1:height), height)
    end do
    
    buffer = transpose(buffer)
  end subroutine fft_2D

  type(Image) function fft(img)
    use img_mod
    type(Image) :: img
    complex(kind=8), allocatable :: buffer(:, :)
    integer :: i, j, x, y
    real(kind=8) :: max, amp, c, processed

    if (img%width == 0 .or. (iand(img%width, img%width - 1) /= 0)) stop "Image width is not a power of 2"
    if (img%height == 0 .or. (iand(img%height, img%height - 1) /= 0)) stop "Image height is not a power of 2"

    allocate(buffer(img%width, img%height))

    do i = 1, img%width
       do j = 1, img%height
          buffer(i, j) = ichar(img%pixels(i, j))
       end do
    end do

    call fft_2d(buffer, img%width, img%height)

    do i = 1, img%width
       do j = 1, img%height
          amp = cdabs(buffer(i, j))
          if (amp > max) then
             max = amp
          end if
       end do
    end do

    c = 255 / dlog(1 + dabs(max))

    fft = make_image(img%width, img%height)
    
    do i = 1, img%width
       do j = 1, img%height
          amp = cdabs(buffer(i, j))
          processed = c * dlog(1 + amp)

          x = merge(i + img%width/2, i - img%width/2, i <= img%width/2)
          y = merge(j + img%height/2, j - img%height/2, j <= img%height/2)
          fft%pixels(x, y) = char(int(processed))
       end do
    end do
  end function fft
end module processing_mod

program main
  use img_mod
  use processing_mod
  implicit none
  
  type(Image) :: img, out

  img = load_pgm("../cameraman.pgm");
  
  out = fft(img)
  
  call save_pgm("out.pgm", img)
end program
