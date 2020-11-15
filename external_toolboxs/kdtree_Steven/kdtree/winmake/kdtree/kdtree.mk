#######################################################################
# Makefile for MatlabCPP
#
#######################################################################

# MATLAB directory -- this may need to change depending on where you have MATLAB installed
MATDIR = C:\\Program Files\\MATLAB\\R2008a

INCDIR = /I "." /I "../../src" -I "$(MATDIR)/extern/include" -I"../Libs/"
CPP = cl
CPPFLAGS = /c /Zp8 /GR /W3 /EHs /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_DEPRECATE \
		 /D_SECURE_SCL=0 /DMATLAB_MEX_FILE /nologo /DWIN64
#CPPFLAGS = /c /Zp8 /MD /GR /W3 /EHs /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_DEPRECATE \#
#	#/D_SECURE_SCL=0 /DMATLAB_MEX_FILE /nologo /D "CPP_ACCEPT_EXPORTS" /D_USERDLL /D_WINDLL \
#	/DUNICODE /D_UNICODE /DWIN64
LINK = link
LINKFLAGS = /dll /export:mexFunction /MAP /MACHINE:X64 \
	/LIBPATH:"$(MATDIR)\extern\lib\win64\microsoft" \
	/LIBPATH:"../Libs/$(OUTDIR)" \
	libmex.lib libmx.lib  libmat.lib 

INSTDIR = ./../../@kdtree/
# DEBUGBUILD = 1

!IF DEFINED(DEBUGBUILD)
OUTDIR = Debug/
CPPFLAGS = $(CPPFLAGS)  /Fo"$(OUTDIR)" /DDEBUG
LINKFLAGS = $(LINKFLAGS) /INCREMENTAL /DEBUG

!ELSE
OUTDIR = Release/
CPPFLAGS = $(CPPFLAGS)  /O2 /Oy- /Fo"$(OUTDIR)" /DNDEBUG 
!ENDIF

.SILENT :

# Rules for making the targets
TARGETS = $(OUTDIR)kdtree.mexw64 \
	$(OUTDIR)kdtree_closestpoint.mexw64 \
	$(OUTDIR)kdtree_range.mexw64

all: $(TARGETS)
	@copy $(OUTDIR:/=\)*.mexw64 $(INSTDIR:/=\)
	@echo Files Built Successfully

clean: 
	@echo Cleaning output filder
	@del $(OUTDIR:/=\)*.mexw64
	@del $(OUTDIR:/=\)*.lib
	@del $(OUTDIR:/=\)*.exp
	@del $(OUTDIR:/=\)*.obj
	@del $(OUTDIR:/=\)*.manifest
	@del $(OUTDIR:/=\)*.map
	@del $(INSTDIR:/=\)*.mexw64
	
rebuild: clean all

.SUFFIXES : mexw64
.SILENT :

# The below two lines don't seem to work -- I'll do it manually
#{$(OUTDIR)}.mexw64{$(OUTDIR)}.obj:
#	$(LINK) $(OUTDIR)$< $(LINKFLAGS) /OUT:$(OUTDIR)$<.mexw64

{./../../src/}.c{$(OUTDIR)}.obj:
    $(CPP) $(CPPFLAGS) $(INCDIR) $<

{./../../src/}.cpp{$(OUTDIR)}.obj:
    $(CPP) $(CPPFLAGS) $(INCDIR) $<


$(OUTDIR)kdtree.mexw64 : $(OUTDIR)kdtree.obj $(OUTDIR)kdtree_create.obj
	$(LINK) $(OUTDIR)kdtree.obj $(OUTDIR)kdtree_create.obj \
	$(LINKFLAGS) /PDB:"$(OUTDIR)kdtree.pdb" \
	/OUT:"$(OUTDIR)kdtree.mexw64"

$(OUTDIR)kdtree_closestpoint.mexw64 : $(OUTDIR)kdtree.obj $(OUTDIR)kdtree_closestpoint.obj
	$(LINK) $(OUTDIR)kdtree.obj $(OUTDIR)kdtree_closestpoint.obj \
	$(LINKFLAGS) /PDB:"$(OUTDIR)kdtree_closestpoint.pdb" \
	/OUT:"$(OUTDIR)kdtree_closestpoint.mexw64"

$(OUTDIR)kdtree_range.mexw64 : $(OUTDIR)kdtree.obj $(OUTDIR)kdtree_range.obj
	$(LINK) $(OUTDIR)kdtree.obj $(OUTDIR)kdtree_range.obj \
	$(LINKFLAGS) /PDB:"$(OUTDIR)kdtree_range.pdb" \
	/OUT:"$(OUTDIR)kdtree_range.mexw64"
