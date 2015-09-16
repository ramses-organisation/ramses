module device_memory
#ifdef CLAUDIO   
   use iso_c_binding
   interface
           subroutine get_dev_mem(total,free) bind(C,name="get_dev_mem")
               use iso_c_binding
               integer(kind=c_size_t)::total,free
           end subroutine get_dev_mem
   end interface
#else
   interface
           subroutine get_dev_mem(total,free) 
               integer::total,free
           end subroutine get_dev_mem
   end interface
#endif

contains

   subroutine print_dev_mem
#ifdef CLAUDIO   
      integer(kind=c_size_t)::total, free
#else
      integer::total,free
#endif
      !call get_dev_mem(total, free)
      total=1
      free =1
      print *, "Free device memory ", real(free)/real(total)*100.0, " %"
   end subroutine print_dev_mem
   
end module device_memory
