#define WSQ_IMPORT __declspec( dllimport )



#if defined(__cplusplus)
extern "C" {
#endif


/* saves the contents of an HBITMAP to a file.  The extension of the file name
// determines the file type.  returns 1 if successfull, 0 otherwise */
extern
WSQ_IMPORT SaveBMPToFile( HBITMAP hBitmap,      /* bitmap to be saved */
                                 const char *filename, int filetype); /* name of output file */


/* Creates an HBITMAP from an image file.  The extension of the file name
// determines the file type.  returns an HBITMAP if successfull, NULL
// otherwise */
extern
HBITMAP WSQ_IMPORT CreateBMPFromFile( const char *filename);


extern
int WSQ_IMPORT RegisterWSQ();



extern
void WSQ_IMPORT WriteWSQ_bitrate(double bitrate);

extern
double WSQ_IMPORT ReadWSQ_bitrate();

extern
void WSQ_IMPORT WriteWSQ_ppi(int ppi);

extern
int WSQ_IMPORT ReadWSQ_ppi();

extern
void WSQ_IMPORT WriteWSQ_comment(char *comment);

extern
char* WSQ_IMPORT ReadWSQ_comment();



extern
void WSQ_IMPORT WriteTIFFcompression(int tiff_compression);

extern
void WSQ_IMPORT WriteTIFFpredictor(int tiff_predictor);

extern
void WSQ_IMPORT SetShowFilePropertiesDialog(int file_properties_dialog);


extern
void WSQ_IMPORT ShowFileConverter();


#if defined(__cplusplus)
 }
#endif

