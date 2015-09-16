#include <cuda.h> // to get memory on the device
#include <cuda_runtime.h> // to get device count

extern "C" 
{
void get_dev_mem(size_t& total, size_t& free) 
{
	cuMemGetInfo(&free, &total);
}
}