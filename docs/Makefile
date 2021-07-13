# Makefile for Sphinx documentation

#Set some variables
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
BUILDDIR      = _build

# Internal variables.
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees  $(SPHINXOPTS) .

.PHONY: help clean html linkcheck doctest api

help:
	@echo " "
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       make only the HTML files from the existing rst sources"
	@echo "  linkcheck  check all external links for integrity"
	@echo "  doctest    run all doctests embedded in the documentation (if enabled)"
	@echo "  clean      remove _build folder"
	@echo "  server     create a server"
	@echo " "

clean:
	rm -rf $(BUILDDIR)
	@echo
	@echo "Remove $(BUILDDIR)"

html:
	@echo
	@echo "Building HTML files"
	@echo
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

doctest:

	$(SPHINXBUILD) -b doctest $(ALLSPHINXOPTS) $(BUILDDIR)/doctest
	@echo "Testing of doctests in the sources finished, look at the " \
	      "results in $(BUILDDIR)/doctest/output.txt."


linkcheck:
	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) $(BUILDDIR)/linkcheck
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in $(BUILDDIR)/linkcheck/output.txt."

server:
	cd $(BUILDDIR)/html && python -m http.server 8001
