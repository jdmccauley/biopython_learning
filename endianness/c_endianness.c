#include <stdio.h>

/* The purpose of this program is to show endianness of the system.
 * You can only see it in multi-byte datatypes.
*/

/* Shows that endianness cannot be seen with characters.
*/
void strings() {
    char* mystr = "ABCDEF";
    
    printf("String:\t%s\n", mystr);

    printf("Size of str:\t%lu\n", sizeof(mystr));
    printf("Size of 'A':\t%lu\n", sizeof(mystr[0]));
    printf("Size of 'F':\t%lu\n", sizeof(mystr[5]));

    printf("Location of string:\t%p\n", (void *)&mystr);
    printf("Location of 'A':\t%p\n", (void *)&mystr[0]);
    printf("Location of 'F':\t%p\n", (void *)&mystr[5]);

}

/* Shows that endinaness can be seen with ints.
 * Big endian = 0x0 0x0 0x0 0x1
 * Little endian = 0x1 0x0 0x0 0x0
 * 
*/
void ints() {
    int x = 1;
    char *loc = 0;

    printf("Value of x:\t%d\n", x);
    printf("Location of x:\t%p\n",(void *)&x);

    printf("Representation of x:\n");
    for (long i = 0; i < sizeof(x); i++) {
        loc = ((char *)&x + i);
        printf("Address\t%p", loc);
        printf("Value\t0x%x\n ", *loc);
    }
    printf("\n");
}


int main(int argc, char *argv[]) {
    strings();
    ints();
    return 0;
}
