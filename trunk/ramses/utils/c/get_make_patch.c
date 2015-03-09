#include <stdio.h>
extern char _binary_make_patch_start[];

int main (void) {
	  puts (_binary_make_patch_start);
		  return 0;
}
