
const char * read_str_value(const char name[], const char *def_value);
int read_int_value(const char name[], int def_value);

unsigned char * read_file(const char path[], int *size);
void write_file(const char *path, unsigned char *data, int width, int height);
