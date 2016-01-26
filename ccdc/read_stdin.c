#include <stdio.h>
#include <stdlib.h>

int main( ) {

   char str[100];
    int i,j;
    int number_of_scenes;
    int number_of_bands=6;
    unsigned char *fmask_buffer;
    int *pixel_buffer;



   //printf( "Enter a value :");
   //scanf("%s %d", str, &i);
    scanf("%d", &number_of_scenes);

    pixel_buffer = (int *) calloc ((number_of_scenes * number_of_bands), sizeof (int));
    if (pixel_buffer == NULL)
    {
        printf ("Allocating pixel_buffer memory\n");
        return (-1);
    }
    fmask_buffer = (char *) calloc (number_of_scenes, sizeof (char));
    if (fmask_buffer == NULL)
    {
        printf ("Allocating fmask_buffer memory\n");
        return (-1);
    }

    for (i = 0; i < number_of_scenes; i++)
    {
        for (j = 0; j < number_of_bands; j++)
        {
            scanf("%d", &pixel_buffer[i+j]);
            printf( "You entered: %d\n", pixel_buffer[i+j]);
        }
        scanf("%d", &fmask_buffer[i]);
        printf( "You entered: %d\n", fmask_buffer[i]);
    }
   
    printf ("\n");


   //printf( "\nYou entered: %s %d ", str, i);

   return 0;
}
