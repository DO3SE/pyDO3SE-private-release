# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build
BUILDDIRPDF      = _build_pdf

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile production

copyfiles:
	cp -f $(BUILDDIRPDF)/* .

convertWordFiles:
	echo "convertWordFiles not implemented"
#	DOCX_FILES=True python convert_word_htm_files.py

buildSource:
# Build source code docstring to rst ignoring tests and snapshots
	rm -f source/*.rst
	SPHINX_APIDOC_OPTIONS=members,undoc-members sphinx-apidoc -f -a -d 1  -M -e -o source ../do3se_phenology ../*[s_][nt][ae][ps][st]*

buildPdf: convertWordFiles
	echo "buildPdf not implemented"
	# @$(SPHINXBUILD) -b pdf "$(SOURCEDIR)" "$(BUILDDIRPDF)" $(O)
	# cp -f ./$(BUILDDIRPDF)/*.pdf .


production: buildPdf convertWordFiles buildSource
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

sourceOnly:buildSource
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)


# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
