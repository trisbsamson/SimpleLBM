
#ifndef VTKOUTPUT_H

void fprintB64(FILE* f, void * tab, int len);

class vtkFileOut {
	FILE * f;
	FILE * fp;
	char * name; int name_size;
	int parallel;
	int size;
public:
	vtkFileOut ();
	int Open(const char* filename);
	void WriteB64(void * tab, int len);
	void Init(char* selection, int dx, int dy, int dz, int nx, int ny, int nz, double spacing);
	void WriteField(char * name, void * data, int elem, char * tp, int components);
	inline void WriteField(char * name, double * data) { WriteField(name, (void*) data, sizeof(double), (char*)"Float64", 1); };
	void Finish();
	void Close();
};

#endif
#define VTKOUTPUT_H 1
