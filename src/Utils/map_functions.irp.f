subroutine map_save_to_disk(filename,map)
  use map_module
  use mmap_module
  implicit none
  character*(*), intent(in)      :: filename
  type(map_type), intent(inout)  :: map
  type(c_ptr)                    :: c_pointer(3)
  integer                        :: fd(3)
  integer*8                      :: i,k
  integer                        :: j


  if (map % consolidated) then
    stop 'map already consolidated'
  endif

  call mmap(trim(filename)//'_consolidated_idx', (/ map % map_size + 2_8 /), 8, fd(1), .False., c_pointer(1))
  call c_f_pointer(c_pointer(1),map % consolidated_idx, (/ map % map_size +2_8/))

  call mmap(trim(filename)//'_consolidated_key', (/ map % n_elements /), cache_key_kind, fd(2), .False., c_pointer(2))
  call c_f_pointer(c_pointer(2),map % consolidated_key, (/ map % n_elements /))

  call mmap(trim(filename)//'_consolidated_value', (/ map % n_elements /), integral_kind, fd(3), .False., c_pointer(3))
  call c_f_pointer(c_pointer(3),map % consolidated_value, (/ map % n_elements /))

  if (.not.associated(map%consolidated_key)) then
    stop 'cannot consolidate map : consolidated_key not associated'
  endif

  if (.not.associated(map%consolidated_value)) then
    stop 'cannot consolidate map : consolidated_value not associated'
  endif

  if (.not.associated(map%consolidated_idx)) then
    stop 'cannot consolidate map : consolidated_idx not associated'
  endif

  call map_sort(map)
  k = 1_8
  do i=0_8, map % map_size
    map % consolidated_idx (i+1) = k
    do j=1, map % map(i) % n_elements
      map % consolidated_value(k) = map % map(i) % value(j)
      map % consolidated_key  (k) = map % map(i) % key(j)
      k = k+1_8
    enddo
    deallocate(map % map(i) % value)
    deallocate(map % map(i) % key)
    map % map(i) % value => map % consolidated_value ( map % consolidated_idx (i+1) :)
    map % map(i) % key   => map % consolidated_key   ( map % consolidated_idx (i+1) :)
  enddo
  map % consolidated_idx (map % map_size + 2_8) = k
  map % consolidated = .True.


!  call munmap( (/ map % map_size + 2_8 /), 8, fd(1), c_pointer(1))
!  call mmap(trim(filename)//'_consolidated_idx', (/ map % map_size + 2_8 /), 8, fd(1), .True., c_pointer(1))
!  call c_f_pointer(c_pointer(1),map % consolidated_idx, (/ map % map_size +2_8/))
!
!  call munmap( (/ map % n_elements /), cache_key_kind, fd(2), c_pointer(2))
!  call mmap(trim(filename)//'_consolidated_key', (/ map % n_elements /), cache_key_kind, fd(2), .True., c_pointer(2))
!  call c_f_pointer(c_pointer(2),map % consolidated_key, (/ map % n_elements /))
!
!  call munmap( (/ map % n_elements /), integral_kind, fd(3), c_pointer(3))
!  call mmap(trim(filename)//'_consolidated_value', (/ map % n_elements /), integral_kind, fd(3), .True., c_pointer(3))
!  call c_f_pointer(c_pointer(3),map % consolidated_value, (/ map % n_elements /))

end

subroutine map_load_from_disk(filename,map)
  use map_module
  use mmap_module
  implicit none
  character*(*), intent(in)      :: filename
  type(map_type), intent(inout)  :: map
  type(c_ptr)                    :: c_pointer(3)
  integer                        :: fd(3)
  integer*8                      :: i,k
  integer                        :: n_elements



  if (map % consolidated) then
    stop 'map already consolidated'
  endif

  call mmap(trim(filename)//'_consolidated_idx', (/ map % map_size + 2_8 /), 8, fd(1), .True., c_pointer(1))
  call c_f_pointer(c_pointer(1),map % consolidated_idx, (/ map % map_size + 2_8/))

  map% n_elements = map % consolidated_idx (map % map_size+2_8)-1

  call mmap(trim(filename)//'_consolidated_key', (/ map % n_elements /), cache_key_kind, fd(2), .True., c_pointer(2))
  call c_f_pointer(c_pointer(2),map % consolidated_key, (/ map % n_elements /))

  call mmap(trim(filename)//'_consolidated_value', (/ map % n_elements /), integral_kind, fd(3), .True., c_pointer(3))
  call c_f_pointer(c_pointer(3),map % consolidated_value, (/ map % n_elements /))

  k = 1_8
  do i=0_8, map % map_size
    deallocate(map % map(i) % value)
    deallocate(map % map(i) % key)
    map % map(i) % value => map % consolidated_value ( map % consolidated_idx (i+1) :)
    map % map(i) % key   => map % consolidated_key   ( map % consolidated_idx (i+1) :)
    map % map(i) % sorted = .True.
    n_elements = map % consolidated_idx (i+2) - k
    k = map % consolidated_idx (i+2)
    map % map(i) % map_size = n_elements
    map % map(i) % n_elements = n_elements
  enddo
  map % n_elements = k-1
  map % sorted = .True. 
  map % consolidated = .True.

end

